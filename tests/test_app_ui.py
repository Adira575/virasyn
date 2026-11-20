from streamlit.testing.v1 import AppTest


def test_sequence_analysis_workspace_renders_by_default():
    app = AppTest.from_file("app.py", default_timeout=60).run()

    assert not app.exception, "\\n".join(
        str(exception.message) for exception in app.exception
    )

    selected = [
        element.value
        for element in app.selectbox
        if getattr(element, "key", None) == "virasyn_workspace_section"
    ]
    assert selected == ["🧬 Sequence & structure analysis"]

    visible_subheaders = [element.value for element in app.subheader]
    assert "Sequence & structure analysis" in visible_subheaders
    assert "Upload genomes, proteins and structures" in visible_subheaders

    visible_captions = [element.value for element in app.caption]
    assert any("virasyn-workspace-ui-2026-11-20-r1" in value for value in visible_captions)


def test_workspace_selector_switches_sections():
    app = AppTest.from_file("app.py", default_timeout=60).run()
    selector = next(
        element
        for element in app.selectbox
        if getattr(element, "key", None) == "virasyn_workspace_section"
    )

    app.selectbox(key="virasyn_workspace_section").set_value(
        "🧰 Bioinformatics workbench"
    ).run()

    visible_subheaders = [element.value for element in app.subheader]
    assert "VIRASYN Bioinformatics Workbench" in visible_subheaders
    assert "Upload genomes, proteins and structures" not in visible_subheaders
