#!/usr/bin/env python3
"""MST Full Pipeline — end-to-end single-page UI.

Upload raw FASTQ.gz files → QC → DADA2 → SourceTracker2 → Pathogen Detection.
Results (source contribution + pathogen abundance) are shown when done.

Run with:
    streamlit run scripts/app2.py
"""

import glob, io, os, shutil, subprocess, sys, tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

st.set_page_config(page_title="MST Full Pipeline", layout="wide")

_HERE   = os.path.dirname(os.path.abspath(__file__))
_DB_DIR = os.path.join(_HERE, "..", "DB")

sys.path.insert(0, _HERE)

from app_sourcetracker_ui import make_figure as _st_figure, _fig_to_png
from app_pathogen_ui import (
    _load_csv_gz, _is_sequence_index, _get_silva_path,
    _classify, _detect, _make_figure as _path_figure,
    PATHOGENS, _PRIORITY_ORDER, _PRIORITY_LABELS,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _run(cmd, spinner_label):
    """Run cmd under a spinner; return (ok, log_text)."""
    with st.spinner(spinner_label):
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return result.returncode == 0, result.stdout


def _log_expander(label, log_text, ok):
    icon = "✅" if ok else "❌"
    with st.expander(f"{icon} {label} — {'OK' if ok else 'FAILED'}", expanded=not ok):
        st.code(log_text, language="")


# ── Session state init ────────────────────────────────────────────────────────

for key in ("done", "st_dir", "asv_table", "st_prop", "path_ra"):
    if key not in st.session_state:
        st.session_state[key] = None


# ── UI layout ─────────────────────────────────────────────────────────────────

st.title("MST Full Pipeline")
st.markdown(
    "Upload raw paired-end FASTQ.gz files. The pipeline runs "
    "**QC → DADA2 → SourceTracker2 → Pathogen Detection** automatically "
    "and shows both result graphs at the end."
)

uploaded = st.file_uploader(
    "FASTQ.gz files (all R1/R2 pairs)",
    type=["gz"],
    accept_multiple_files=True,
)

with st.expander("Options", expanded=True):
    _ALL_GROUPS   = ["bat", "chicken", "cow", "duck", "groundwater",
                     "horse", "human", "pig", "river", "seawater"]
    _DEFAULT_ON   = {"bat", "chicken", "cow", "duck", "human", "pig", "river"}
    st.markdown("**Source groups**")
    _gcols = st.columns(len(_ALL_GROUPS))
    _sel   = [g for i, g in enumerate(_ALL_GROUPS)
              if _gcols[i].checkbox(g, value=(g in _DEFAULT_ON), key=f"grp2_{g}")]
    groups = ",".join(_sel)
    if not _sel:
        st.warning("Select at least one source group.")

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        mode      = st.selectbox("Feature mode", ["asv", "otu"])
    with c2:
        src_depth = st.number_input("Source rarefaction depth", value=1000, min_value=100)
        snk_depth = st.number_input("Sink rarefaction depth (0=auto)", value=0, min_value=0)
    with c3:
        restarts  = st.number_input("MCMC restarts",      value=10,  min_value=1)
        burnin    = st.number_input("Burnin iterations",  value=100, min_value=1)
        draws     = st.number_input("Draws per restart",  value=1,   min_value=1)
    with c4:
        threads   = st.number_input("Threads",            value=4,   min_value=1, max_value=256)
        pat_id    = st.slider("Pathogen SILVA identity (%)", 90, 100, 97) / 100

run_btn = st.button("▶  Run Pipeline", type="primary", disabled=not uploaded)

if not uploaded:
    st.info("Upload FASTQ.gz files to enable the pipeline.")

# ── Pipeline execution ────────────────────────────────────────────────────────

if run_btn and uploaded:
    st.session_state.done = False

    # Fresh working directory
    if st.session_state.get("_work_dir") and os.path.exists(st.session_state["_work_dir"]):
        shutil.rmtree(st.session_state["_work_dir"])
    work_dir = tempfile.mkdtemp(prefix="mst_run_")
    st.session_state["_work_dir"] = work_dir

    sink_dir  = os.path.join(work_dir, "sink")
    out_dir   = os.path.join(work_dir, "output")
    qc_dir    = os.path.join(out_dir, "01_qc")
    dada2_dir = os.path.join(out_dir, "02_dada2")
    st_dir    = os.path.join(out_dir, "03_sourcetracker")
    for d in (sink_dir, qc_dir, dada2_dir, st_dir):
        os.makedirs(d)

    # Save uploaded files
    for f in uploaded:
        with open(os.path.join(sink_dir, f.name), "wb") as fh:
            fh.write(f.read())

    fastq_files = sorted(glob.glob(os.path.join(sink_dir, "*.fastq.gz")))

    st.divider()
    st.subheader("Pipeline Progress")

    # ── Step 1: QC ────────────────────────────────────────────────────────────
    ok, log = _run(
        [sys.executable, os.path.join(_HERE, "app_qc.py")]
        + fastq_files
        + ["-o", qc_dir, "-t", str(int(threads))],
        "Step 1/3 — QC & V4 Extraction …",
    )
    _log_expander("Step 1/3 — QC & V4 Extraction", log, ok)
    if not ok:
        st.error("QC step failed.")
        st.stop()

    trimmed = sorted(glob.glob(os.path.join(qc_dir, "*_trimmed_R*.fastq.gz")))
    if not trimmed:
        st.error("No trimmed files found after QC.")
        st.stop()

    # ── Step 2: DADA2 ─────────────────────────────────────────────────────────
    ok, log = _run(
        [sys.executable, os.path.join(_HERE, "app_dada2.py")]
        + trimmed
        + ["-o", dada2_dir, "--threads", str(int(threads))],
        "Step 2/3 — DADA2 ASV Table …",
    )
    _log_expander("Step 2/3 — DADA2 ASV Table", log, ok)
    if not ok:
        st.error("DADA2 step failed.")
        st.stop()

    # Pick the right ASV table
    merged = os.path.join(dada2_dir, "ASV_table_merged.csv.gz")
    if os.path.exists(merged):
        asv_table = merged
    else:
        candidates = sorted(
            t for t in glob.glob(os.path.join(dada2_dir, "ASV_table_*.csv.gz"))
            if "merged" not in t
        )
        if not candidates:
            st.error("No ASV table produced by DADA2.")
            st.stop()
        asv_table = candidates[0]

    # ── Step 3: SourceTracker ─────────────────────────────────────────────────
    st_cmd = [
        sys.executable, os.path.join(_HERE, "run_sourcetracker.py"),
        asv_table,
        os.path.join(_DB_DIR, "db_table.csv.gz"),
        os.path.join(_DB_DIR, "db.fasta"),
        os.path.join(_DB_DIR, "MST.design"),
        st_dir,
        "--mode",      mode,
        "--src-depth", str(int(src_depth)),
        "--snk-depth", str(int(snk_depth)),
        "--restarts",  str(int(restarts)),
        "--burnin",    str(int(burnin)),
        "--draws",     str(int(draws)),
        "--threads",   str(int(threads)),
    ]
    if groups.strip():
        st_cmd += ["--groups"] + [g.strip() for g in groups.split(",") if g.strip()]

    ok, log = _run(st_cmd, "Step 3/3 — SourceTracker2 …")
    _log_expander("Step 3/3 — SourceTracker2", log, ok)
    if not ok:
        st.error("SourceTracker2 failed.")
        st.stop()

    # ── Step 4: Pathogen Detection (inline) ───────────────────────────────────
    with st.spinner("Step 4/4 — Pathogen Detection …"):
        asv_df = _load_csv_gz(asv_table)
        if _is_sequence_index(asv_df.index) and os.path.exists(
            os.path.join(_HERE, "..", "config", "silva_nr99_v138.1_train_set_v4.fa.gz")
        ):
            silva = _get_silva_path()
            seq_to_genus = _classify(list(asv_df.index), silva, pat_id, int(threads))
            _, ra_df = _detect(asv_df, seq_to_genus)
        else:
            ra_df = pd.DataFrame()

    st.success("✅ All steps complete")

    # Store results in session state
    st.session_state.done      = True
    st.session_state.st_dir    = st_dir
    st.session_state.asv_table = asv_table
    st.session_state.st_prop   = pd.read_csv(
        os.path.join(st_dir, "proportions.csv"), index_col=0
    )
    st.session_state.path_ra   = ra_df if not ra_df.empty else None


# ── Results ───────────────────────────────────────────────────────────────────

if st.session_state.done:
    st.divider()
    st.header("Results")

    col_st, col_path = st.columns(2)

    # ── Source contributions ──────────────────────────────────────────────────
    with col_st:
        st.subheader("Source Contributions")
        mp = st.session_state.st_prop
        fig = _st_figure(mp, "Source Contribution (SourceTracker2)")
        png = _fig_to_png(fig)
        st.image(png)
        plt.close(fig)

        pct = (mp * 100).round(2)
        pct.index.name = "Sample"
        st.dataframe(pct, use_container_width=True)
        st.download_button(
            "Download contributions CSV",
            pct.to_csv().encode(),
            file_name="contributions_pct.csv",
            mime="text/csv",
        )

    # ── Pathogen detection ────────────────────────────────────────────────────
    with col_path:
        st.subheader("Pathogen Detection")
        ra_df = st.session_state.path_ra

        if ra_df is None or ra_df.empty:
            st.info("No pathogenic bacteria detected (or SILVA database unavailable).")
        else:
            fig = _path_figure(ra_df)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
            buf.seek(0)
            st.image(buf)
            plt.close(fig)
            st.download_button(
                "Download pathogen figure (PNG)", buf,
                file_name="pathogen_abundance.png", mime="image/png",
            )

            ra_df.index.name = "Sample"
            st.dataframe(
                ra_df.round(4).style.background_gradient(cmap="Reds", axis=None),
                use_container_width=True,
            )
            st.download_button(
                "Download pathogen RA CSV",
                ra_df.to_csv().encode(),
                file_name="pathogen_ra.csv",
                mime="text/csv",
            )
