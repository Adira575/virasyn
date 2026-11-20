"""VIRASYN baseline RBP–receptor pair model.

The model is intentionally a transparent baseline. It must be trained on
experimentally labelled positive and tested-negative pairs. Unknown pairs
must not be encoded as negatives. Evaluation is group-held-out to reduce
sequence/lineage leakage; reported metrics are not a substitute for external
validation on an independent dataset.
"""
from __future__ import annotations

from io import StringIO
import math
import re

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    confusion_matrix, f1_score, precision_score, recall_score, roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit

AA = "ACDEFGHIKLMNPQRSTVWY"
REQUIRED_COLUMNS = ["rbp_sequence", "receptor_sequence", "label", "split_group"]


def training_template() -> str:
    return (
        "pair_id,rbp_sequence,receptor_sequence,label,split_group,evidence,citation\n"
        "pair_001,ACDEFGHIKLMNPQRSTVWY,ACDEFGHIKLMNPQRSTVWY,1,rbp_cluster_A,experimental_binding,doi:...\n"
        "pair_002,ACDEFGHIKLMNPQRSTVWY,YYYYYYYYYYYYYYYYYYYY,0,rbp_cluster_B,experimentally_tested_nonbinding,doi:...\n"
    )


def _clean_sequence(value: object, column: str) -> str:
    sequence = re.sub(r"\s+", "", str(value or "")).upper().replace("*", "")
    if not sequence:
        raise ValueError(f"{column} contains an empty protein sequence.")
    invalid = sorted(set(sequence) - set(AA))
    if invalid:
        raise ValueError(f"{column} contains unsupported amino-acid symbols: {''.join(invalid)}")
    if len(sequence) < 5:
        raise ValueError(f"{column} sequences must contain at least 5 amino acids.")
    return sequence


def read_training_pairs(csv_text: str) -> pd.DataFrame:
    try:
        frame = pd.read_csv(StringIO(csv_text))
    except Exception as exc:
        raise ValueError(f"Training-pair CSV could not be parsed: {exc}") from exc
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError("Training CSV is missing required columns: " + ", ".join(missing))
    if len(frame) < 10:
        raise ValueError("At least 10 labelled pairs are required for a preliminary grouped evaluation.")
    frame = frame.copy()
    frame["rbp_sequence"] = frame["rbp_sequence"].map(lambda x: _clean_sequence(x, "rbp_sequence"))
    frame["receptor_sequence"] = frame["receptor_sequence"].map(lambda x: _clean_sequence(x, "receptor_sequence"))
    labels = pd.to_numeric(frame["label"], errors="coerce")
    if labels.isna().any() or not set(labels.unique()).issubset({0, 1}):
        raise ValueError("label must contain only 1 (experimentally supported binding) or 0 (experimentally tested non-binding).")
    frame["label"] = labels.astype(int)
    if frame["label"].nunique() != 2:
        raise ValueError("Training data must contain both positive and negative labelled pairs.")
    frame["split_group"] = frame["split_group"].astype(str).str.strip()
    if (frame["split_group"] == "").any():
        raise ValueError("Every row needs a split_group such as RBP cluster, phage lineage or a stricter independent group.")
    if frame["split_group"].nunique() < 4:
        raise ValueError("At least 4 distinct split_group values are required for grouped holdout evaluation.")
    if frame["label"].value_counts().min() < 3:
        raise ValueError("At least 3 examples of each label are required for preliminary training.")
    return frame.reset_index(drop=True)


def _composition(sequence: str) -> list[float]:
    length = float(len(sequence))
    return [sequence.count(aa) / length for aa in AA]


def pair_features(rbp_sequence: str, receptor_sequence: str) -> list[float]:
    rbp = _clean_sequence(rbp_sequence, "rbp_sequence")
    receptor = _clean_sequence(receptor_sequence, "receptor_sequence")
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    score = float(aligner.score(rbp, receptor))
    self_rbp = max(float(aligner.score(rbp, rbp)), 1.0)
    self_receptor = max(float(aligner.score(receptor, receptor)), 1.0)
    normalized_score = max(0.0, min(1.0, score / math.sqrt(self_rbp * self_receptor)))
    alignment = aligner.align(rbp, receptor)[0]
    aligned_pairs = sum(int(end - start) for start, end in alignment.aligned[0])
    identity = 0.0
    if aligned_pairs:
        identity = sum(
            rbp[a0 + offset] == receptor[b0 + offset]
            for (a0, a1), (b0, b1) in zip(alignment.aligned[0], alignment.aligned[1])
            for offset in range(min(a1 - a0, b1 - b0))
        ) / aligned_pairs
    features = (
        _composition(rbp) + _composition(receptor) +
        [math.log1p(len(rbp)), math.log1p(len(receptor)),
         min(len(rbp), len(receptor)) / max(len(rbp), len(receptor)),
         normalized_score, identity]
    )
    return features


def _feature_matrix(frame: pd.DataFrame) -> np.ndarray:
    return np.asarray([
        pair_features(row.rbp_sequence, row.receptor_sequence)
        for row in frame.itertuples(index=False)
    ], dtype=float)


def train_group_heldout_model(csv_text: str, test_size: float = 0.25, random_state: int = 42) -> dict:
    """Fit a baseline RF model and evaluate on unseen split_group values."""
    frame = read_training_pairs(csv_text)
    X = _feature_matrix(frame)
    y = frame["label"].to_numpy()
    groups = frame["split_group"].to_numpy()

    splitter = GroupShuffleSplit(n_splits=100, test_size=test_size, random_state=random_state)
    selected = None
    for train_idx, test_idx in splitter.split(X, y, groups):
        if len(set(y[train_idx])) == 2 and len(set(y[test_idx])) == 2:
            selected = (train_idx, test_idx)
            break
    if selected is None:
        raise ValueError(
            "Could not create a group-held-out split containing both labels in train and test. "
            "Add more independent groups and examples of each class; do not fall back to a random row split."
        )
    train_idx, test_idx = selected
    model = RandomForestClassifier(
        n_estimators=400, max_features="sqrt", min_samples_leaf=2,
        class_weight="balanced_subsample", random_state=random_state,
        n_jobs=-1,
    )
    model.fit(X[train_idx], y[train_idx])
    probabilities = model.predict_proba(X[test_idx])[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    y_test = y[test_idx]
    metrics = {
        "n_pairs": int(len(frame)),
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "n_groups_total": int(pd.Series(groups).nunique()),
        "n_groups_train": int(pd.Series(groups[train_idx]).nunique()),
        "n_groups_test": int(pd.Series(groups[test_idx]).nunique()),
        "train_groups": sorted(set(groups[train_idx])),
        "test_groups": sorted(set(groups[test_idx])),
        "accuracy": float(accuracy_score(y_test, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, predictions)),
        "precision": float(precision_score(y_test, predictions, zero_division=0)),
        "recall": float(recall_score(y_test, predictions, zero_division=0)),
        "f1": float(f1_score(y_test, predictions, zero_division=0)),
        "roc_auc": float(roc_auc_score(y_test, probabilities)),
        "average_precision": float(average_precision_score(y_test, probabilities)),
        "confusion_matrix": confusion_matrix(y_test, predictions, labels=[0, 1]).tolist(),
        "evaluation": "single group-held-out split; preliminary internal evaluation only",
        "calibrated": False,
        "external_validation": False,
    }
    prediction_table = frame.iloc[test_idx][[c for c in ["pair_id", "split_group", "evidence", "citation"] if c in frame.columns]].copy()
    prediction_table["true_label"] = y_test
    prediction_table["predicted_label"] = predictions
    prediction_table["binding_score_uncalibrated"] = probabilities
    return {"model": model, "metrics": metrics, "test_predictions": prediction_table}


def predict_pairs(model, csv_text: str) -> pd.DataFrame:
    try:
        frame = pd.read_csv(StringIO(csv_text))
    except Exception as exc:
        raise ValueError(f"Prediction-pair CSV could not be parsed: {exc}") from exc
    missing = [column for column in ["rbp_sequence", "receptor_sequence"] if column not in frame.columns]
    if missing:
        raise ValueError("Prediction CSV is missing columns: " + ", ".join(missing))
    frame = frame.copy()
    frame["rbp_sequence"] = frame["rbp_sequence"].map(lambda x: _clean_sequence(x, "rbp_sequence"))
    frame["receptor_sequence"] = frame["receptor_sequence"].map(lambda x: _clean_sequence(x, "receptor_sequence"))
    probabilities = model.predict_proba(_feature_matrix(frame))[:, 1]
    frame["binding_score_uncalibrated"] = probabilities
    frame["model_output_label_at_0_5_threshold"] = (probabilities >= 0.5).astype(int)
    frame["interpretation"] = "uncalibrated model score; experimental validation required"
    return frame
