import pytest

from virasyn.rbp_model import (
    pair_features, read_training_pairs, train_group_heldout_model,
)

def _dataset():
    rows = ["pair_id,rbp_sequence,receptor_sequence,label,split_group,evidence,citation"]
    positives = [
        ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
        ("VVVVACDEFGHIKLMNPQRS", "VVVVACDEFGHIKLMNPQRS"),
        ("ACDEACDEFGHIKLMNPQRS", "ACDEACDEFGHIKLMNPQRS"),
        ("MNPQRSTVWYACDEFGHIKL", "MNPQRSTVWYACDEFGHIKL"),
        ("FGHIKLMNPQRSTVWYACDE", "FGHIKLMNPQRSTVWYACDE"),
    ]
    negatives = [
        ("YYYYYYYYYYYYYYYYYYYY", "ACDEFGHIKLMNPQRSTVWY"),
        ("WWWWWWWWWWWWWWWWWWWW", "VVVVACDEFGHIKLMNPQRS"),
        ("KKKKKKKKKKKKKKKKKKKK", "ACDEACDEFGHIKLMNPQRS"),
        ("DDDDDDDDDDDDDDDDDDDD", "MNPQRSTVWYACDEFGHIKL"),
        ("EEEEEEEEEEEEEEEEEEEE", "FGHIKLMNPQRSTVWYACDE"),
    ]
    for i, (rbp, receptor) in enumerate(positives):
        rows.append(f"p{i},{rbp},{receptor},1,group_{i},experimental,doi:p{i}")
    for i, (rbp, receptor) in enumerate(negatives):
        rows.append(f"n{i},{rbp},{receptor},0,group_{i+5},tested_nonbinding,doi:n{i}")
    return "\n".join(rows) + "\n"

def test_pair_features_are_finite_and_fixed_length():
    values = pair_features("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY")
    assert len(values) == 45
    assert all(value == value for value in values)

def test_training_data_requires_both_classes_and_groups():
    frame = read_training_pairs(_dataset())
    assert len(frame) == 10
    with pytest.raises(ValueError, match="both positive and negative"):
        read_training_pairs(_dataset().replace(",0,", ",1,"))

def test_group_heldout_model_never_overlaps_groups():
    result = train_group_heldout_model(_dataset(), test_size=0.3, random_state=7)
    metrics = result["metrics"]
    assert set(metrics["train_groups"]).isdisjoint(metrics["test_groups"])
    assert metrics["n_test"] > 0
    assert 0 <= metrics["roc_auc"] <= 1
