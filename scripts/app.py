#!/usr/bin/env python3
"""MST Pipeline — combined web UI.

Combines QC & V4 Extraction, DADA2 Analysis, and SourceTracker
into a single multi-tab Streamlit application.

Run with:
    streamlit run scripts/app.py
"""

import streamlit as st

# Must be the very first Streamlit call.
st.set_page_config(page_title="MST Pipeline", layout="wide")

# Prevent sub-apps from calling set_page_config again (only one call allowed).
st.set_page_config = lambda *a, **kw: None

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import app_qc_ui
import app_dada2_ui
import app_sourcetracker_ui
import app_pathogen_ui

tab_qc, tab_dada2, tab_st, tab_pathogen = st.tabs(
    ["QC & Extraction", "DADA2", "SourceTracker", "Pathogen Detection"]
)

with tab_qc:
    app_qc_ui.main()

with tab_dada2:
    app_dada2_ui.main()

with tab_st:
    app_sourcetracker_ui.main()

with tab_pathogen:
    app_pathogen_ui.main()
