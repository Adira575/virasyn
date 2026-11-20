from pathlib import Path

from streamlit.testing.v1 import AppTest


APP_PATH = Path(__file__).resolve().parents[1] / "app.py"
TAB_LABELS = [
    "🔎 Research database",
    "🧬 Sequence & structure analysis",
    "🧪 PCR / qPCR analysis",
    "🧰 Bioinformatics workbench",
    "🧬 AlphaFold",
    "🦠 PhARIS RBP intelligence",
    "🧫 Phage intelligence",
    "📚 Database directory",
    "ℹ️ About VIRASYN",
]


def test_all_nine_horizontal_workspaces_render():
    app = AppTest.from_file(APP_PATH, default_timeout=60).run()

    assert not app.exception, "\\n".join(
        str(exception.message) for exception in app.exception
    )

    source = APP_PATH.read_text(encoding="utf-8")
    tabs_start = source.index("database_tab, analysis_tab, qpcr_tab, workbench_tab, alphafold_tab, pharis_tab, phage_tab, resources_tab, about_tab = st.tabs(")
    tabs_end = source.index("\n)", tabs_start)
    tabs_declaration = source[tabs_start:tabs_end]
    positions = [tabs_declaration.index(label) for label in TAB_LABELS]
    assert positions == sorted(positions)

    visible_subheaders = [element.value for element in app.subheader]
    expected_headings = [
        "Search viral research literature and molecular data",
        "Sequence & structure analysis",
        "Real-time PCR (qPCR / RT-qPCR) analysis",
        "VIRASYN Bioinformatics Workbench",
        "AlphaFold protein structure lookup",
        "PhARIS-inspired RBP intelligence",
        "VIRASYN Phage Intelligence · Bacteria",
        "A wider directory of public virology and molecular-biology resources",
        "A unified discovery and analysis workspace for virology",
    ]
    for heading in expected_headings:
        assert heading in visible_subheaders

    assert any("virasyn-rbp-model-2026-11-20-r7" in item.value for item in app.caption)


def test_sequence_analysis_upload_controls_are_rendered():
    app = AppTest.from_file(APP_PATH, default_timeout=60).run()

    assert not app.exception, "\\n".join(
        str(exception.message) for exception in app.exception
    )

    labels = [element.label for element in app.file_uploader]
    assert "Genome/protein FASTA and PDB/mmCIF structures" in labels
    assert "Optional BLASTP reference proteins (FASTA)" in labels
    assert "Optional HMMER profile library (.hmm)" in labels
