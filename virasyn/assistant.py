from __future__ import annotations

import os

import streamlit as st
from streamlit_float import float_init


SYSTEM_PROMPT = """You are VIRASYN AI, a scientific help assistant inside the VIRASYN viral research platform.
Help users navigate the platform, understand database records, sequence analysis, ORF tables,
protein properties, BLASTP/HMMER outputs and structural-biology terminology.
Be careful and evidence-led. Never invent citations, accession IDs, database records, experimental results,
host-range predictions, receptor interactions or confidence scores. Distinguish similarity from function,
attachment from productive infection, and a computational candidate from experimentally validated biology.
The platform is a research aid, not a diagnostic or clinical decision system.
When a user asks about their visible results, rely only on the context supplied to you; say when it is insufficient.
Keep answers clear, practical and concise. If unsure, explain the uncertainty.
"""


def _get_secret(name: str) -> str:
    value = os.getenv(name, "")
    if value:
        return value
    try:
        return str(st.secrets.get(name, ""))
    except Exception:
        return ""


def _fallback_answer(question: str, context: str) -> str:
    q = question.lower()
    if any(word in q for word in ["roadmap", "future", "next development", "planned feature"]):
        return (
            "The next planned VIRASYN development areas are individual virus profile pages, saved research "
            "collections, automated genome annotation, and scientifically validated virus-host prediction "
            "models. These are roadmap items, not features in the current release. Host-prediction models "
            "will need virus-group-appropriate datasets and independent validation; a single universal model "
            "would not be scientifically appropriate."
        )
    if any(word in q for word in ["what can", "features", "how to use", "help"]):
        return (
            "I can help you use VIRASYN. The Research Database searches literature, nucleotide and protein "
            "records, taxonomy, structures and protein-family resources. Sequence & Structure Analysis accepts "
            "FASTA and PDB/mmCIF files, screens six-frame ORFs, calculates protein properties and can run "
            "BLASTP/HMMER when the server tools and reference files are available. The Database Directory "
            "links to additional specialist resources."
        )
    if any(word in q for word in ["orf", "reading frame", "six frame"]):
        return (
            "Six-frame ORF screening scans the three forward and three reverse-complement reading frames for "
            "start-to-stop coding candidates. It is a preliminary screen, not a full gene caller. Complex "
            "animal and plant viral genomes may require splice-aware or virus-specific annotation."
        )
    if any(word in q for word in ["blast", "similarity", "identity", "coverage"]):
        return (
            "BLASTP compares an input protein against the reference FASTA you provide. Interpret identity, "
            "query coverage, E-value and bit score together. A significant sequence match supports homology; "
            "it does not independently prove receptor binding, host range or biological function."
        )
    if any(word in q for word in ["hmmer", "hmm", "domain", "pfam", "interpro"]):
        return (
            "HMMER searches your uploaded profile library for conserved protein-family/domain patterns. "
            "The quality and scope of the supplied HMM library determine what can be detected. InterPro "
            "records in the research database are annotations for reference entries; they are not a substitute "
            "for running a profile scan on your own sequence."
        )
    if any(word in q for word in ["host", "receptor", "infect", "tropism"]):
        return (
            "VIRASYN does not currently provide a validated universal host-range or receptor-prediction model. "
            "A predicted attachment protein, sequence similarity or receptor-binding hypothesis is not proof "
            "of productive infection. Use organism-specific experimental evidence and report uncertainty."
        )
    if any(word in q for word in ["pdb", "structure", "alphafold", "3d"]):
        return (
            "Upload PDB or mmCIF coordinate files in the analysis tab. VIRASYN reports chain/residue metadata, "
            "and PDB uploads can be shown in the interactive 3D viewer. Automated structure prediction and "
            "structure alignment are not part of the current release."
        )
    return (
        "Quick-help mode is active, but generative AI has not been connected yet. The platform owner can enable "
        "the full assistant by adding OPENAI_API_KEY to Streamlit Community Cloud → App settings → Secrets. "
        "That key is kept server-side.\n\n"
        "I can still explain ORF screening, BLASTP, HMMER, the research databases and the current analysis "
        "limitations.\n\n"
        + (f"Current page context: {context[:900]}" if context else "")
    )


def _generate_answer(question: str, history: list[dict], context: str) -> str:
    api_key = _get_secret("OPENAI_API_KEY") or _get_secret("VIRASYN_OPENAI_API_KEY")
    if not api_key:
        return _fallback_answer(question, context)
    try:
        from openai import OpenAI
        model = _get_secret("VIRASYN_AI_MODEL") or "gpt-6.2-astra"
        client = OpenAI(api_key=api_key)
        recent = history[-10:]
        input_messages = [{"role": item["role"], "content": item["content"]} for item in recent]
        if not input_messages or input_messages[-1].get("role") != "user" or input_messages[-1].get("content") != question:
            input_messages.append({"role": "user", "content": question})
        if context:
            instructions = SYSTEM_PROMPT + "\n\nCurrent VIRASYN page/search context (unverified user-visible data):\n" + context[:5000]
        else:
            instructions = SYSTEM_PROMPT
        response = client.responses.create(
            model=model,
            instructions=instructions,
            input=input_messages,
        )
        answer = getattr(response, "output_text", "") or ""
        return answer.strip() or "I couldn't generate an answer for that question. Please try again."
    except Exception as exc:
        return (
            "The AI service could not complete this request. Check that the API key and configured model "
            "are valid in Streamlit Secrets. Technical detail: " + str(exc)[:500]
        )


def render_floating_assistant() -> None:
    """Render a persistent bottom-right assistant bubble and chat panel."""
    try:
        float_init()
    except Exception:
        return

    if "virasyn_ai_open" not in st.session_state:
        st.session_state["virasyn_ai_open"] = False
    if "virasyn_ai_messages" not in st.session_state:
        st.session_state["virasyn_ai_messages"] = [
            {
                "role": "assistant",
                "content": "Hello! I'm VIRASYN AI. I can help you search research databases, understand viral sequence results and navigate the platform.",
            }
        ]

    toggle = st.container()
    with toggle:
        label = "✕ Close assistant" if st.session_state["virasyn_ai_open"] else "✦ Ask VIRASYN AI"
        if st.button(label, key="virasyn_ai_toggle", type="primary"):
            st.session_state["virasyn_ai_open"] = not st.session_state["virasyn_ai_open"]
            st.rerun()
    toggle.float(
        "position:fixed;right:1.15rem;bottom:1.1rem;z-index:999999;"
        "border-radius:999px;box-shadow:0 8px 28px rgba(0,0,0,.22);"
    )

    if st.session_state["virasyn_ai_open"]:
        # A keyed container lets the card and all of its descendants receive an
        # opaque surface, rather than only positioning a transparent wrapper.
        st.markdown(
            """
            <style>
            .st-key-virasyn_ai_panel {
                width: min(390px, calc(100vw - 2rem)) !important;
                max-width: calc(100vw - 2rem) !important;
                box-sizing: border-box !important;
                padding: .85rem .9rem !important;
                background: #ffffff !important;
                color: #172333 !important;
                border: 1px solid #d7dee8 !important;
                border-radius: 18px !important;
                box-shadow: 0 14px 42px rgba(15, 23, 42, .28) !important;
                isolation: isolate !important;
                opacity: 1 !important;
            }
            .st-key-virasyn_ai_panel > div,
            .st-key-virasyn_ai_panel [data-testid="stVerticalBlock"] {
                background: #ffffff !important;
                color: #172333 !important;
            }
            .st-key-virasyn_ai_panel [data-testid="stChatMessage"] {
                background: #f3f6fa !important;
                color: #172333 !important;
                border-radius: 12px !important;
                padding: .55rem .7rem !important;
            }
            .st-key-virasyn_ai_panel input {
                background: #ffffff !important;
                color: #172333 !important;
            }
            .st-key-virasyn_ai_panel [data-testid="stForm"] {
                background: #ffffff !important;
                border: 0 !important;
                padding: .2rem 0 0 !important;
            }
            @media (max-width: 480px) {
                .st-key-virasyn_ai_panel {
                    width: calc(100vw - 1rem) !important;
                    max-width: calc(100vw - 1rem) !important;
                    padding: .7rem !important;
                }
            }
            </style>
            """,
            unsafe_allow_html=True,
        )
        panel = st.container(key="virasyn_ai_panel", border=False)
        with panel:
            st.markdown("### ✦ VIRASYN AI")
            st.caption("Research guide · Database help · Bioinformatics explanations")
            if not (_get_secret("OPENAI_API_KEY") or _get_secret("VIRASYN_OPENAI_API_KEY")):
                st.info("Quick-help mode. Full generative answers can be enabled by the site owner with a server-side OpenAI API key.")
            chat_area = st.container(height=245)
            with chat_area:
                for message in st.session_state["virasyn_ai_messages"][-12:]:
                    with st.chat_message(message["role"]):
                        st.markdown(message["content"])
            count = len(st.session_state["virasyn_ai_messages"])
            with st.form(f"virasyn_ai_form_{count}", clear_on_submit=True):
                prompt = st.text_input(
                    "Ask a question",
                    placeholder="e.g. What does query coverage mean?",
                    label_visibility="collapsed",
                )
                sent = st.form_submit_button("Send", type="primary", use_container_width=True)
            if sent and prompt.strip():
                messages = st.session_state["virasyn_ai_messages"]
                messages.append({"role": "user", "content": prompt.strip()})
                answer = _generate_answer(
                    prompt.strip(),
                    messages[:-1],
                    st.session_state.get("virasyn_ai_context", ""),
                )
                messages.append({"role": "assistant", "content": answer})
                st.session_state["virasyn_ai_messages"] = messages[-20:]
                st.rerun()
            st.caption("AI answers may be incomplete. Verify scientific claims against original records and papers.")
        panel.float(
            "position:fixed!important;right:1rem!important;bottom:5.15rem!important;"
            "width:min(390px,calc(100vw - 2rem))!important;"
            "max-width:calc(100vw - 2rem)!important;"
            "max-height:calc(100vh - 6.5rem)!important;overflow-y:auto!important;"
            "overscroll-behavior:contain;z-index:999998!important;"
            "background-color:#ffffff!important;color:#172333!important;"
            "border:1px solid #d7dee8!important;border-radius:18px!important;"
            "box-shadow:0 14px 42px rgba(15,23,42,.28)!important;"
            "isolation:isolate!important;opacity:1!important;box-sizing:border-box!important;"
        )
