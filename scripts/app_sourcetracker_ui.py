#!/usr/bin/env python3
"""Streamlit SourceTracker UI — microbial source tracking.

Accepts:
  • Sink table  : CSV.gz — rows = ASVs (sequences or IDs), cols = samples.
  • Source table: CSV.gz — rows = ASV IDs, cols = samples.
  • Source FASTA: sequences for the source ASV IDs.
  • Source design: TSV/TXT  Sample<tab>Group  (one per line).

Feature alignment options:
  • ASV  — 100 % exact match
  • OTU  — 99 % identity via vsearch

Gibbs sampling runs inside the ST conda environment via subprocess.
"""

import io
import os
import re
import subprocess
import tempfile
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import streamlit as st

# Paths resolved relative to this script's directory
_HERE         = os.path.dirname(os.path.abspath(__file__))
_GIBBS_SCRIPT = os.path.join(_HERE, "_run_gibbs.py")
_DB_DIR       = os.path.join(_HERE, "..", "DB")
_SOURCE_TABLE = os.path.join(_DB_DIR, "db_table.csv.gz")
_SOURCE_FASTA = os.path.join(_DB_DIR, "db.fasta")
_SOURCE_DESIGN= os.path.join(_DB_DIR, "MST.design")

# ── Colour palette ────────────────────────────────────────────────────────────
_FALLBACK_COLORS = [
    "#E07B7B", "#A0785A", "#F0C040", "#5DB0A0", "#9B6BB5",
    "#C8A06E", "#4A90D9", "#90C8F0", "#2060A0", "#60C8C0",
    "#E8A060", "#A0D090", "#D070A0", "#70B0D0",
]
GROUP_COLORS = {
    "pig":         "#E07B7B",
    "cow":         "#A0785A",
    "chicken":     "#F0C040",
    "duck":        "#5DB0A0",
    "bat":         "#9B6BB5",
    "horse":       "#C8A06E",
    "human":       "#4A90D9",
    "groundwater": "#90C8F0",
    "seawater":    "#2060A0",
    "river":       "#60C8C0",
    "Unknown":     "#C8C8C8",
}


# ════════════════════════════════════════════════════════════════════════════
#  Data loaders
# ════════════════════════════════════════════════════════════════════════════

def _load_csv_gz(file_obj) -> pd.DataFrame:
    if hasattr(file_obj, "seek"):
        file_obj.seek(0)
        data = file_obj.read()
    else:
        with open(file_obj, "rb") as f:
            data = f.read()
    return pd.read_csv(io.BytesIO(data), compression="gzip", index_col=0)


def _load_fasta(file_obj) -> dict:
    if hasattr(file_obj, "seek"):
        file_obj.seek(0)
        text = file_obj.read().decode("utf-8", errors="replace")
    else:
        with open(file_obj) as f:
            text = f.read()
    seqs, cur_id, cur_seq = {}, None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if cur_id:
                seqs[cur_id] = "".join(cur_seq)
            cur_id, cur_seq = line[1:].split()[0], []
        elif cur_id:
            cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs


def _load_design(file_obj) -> dict:
    if hasattr(file_obj, "seek"):
        file_obj.seek(0)
        text = file_obj.read().decode("utf-8", errors="replace")
    else:
        with open(file_obj) as f:
            text = f.read()
    design = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.lower().startswith("sample"):
            continue
        parts = s.split("\t")
        if len(parts) >= 2:
            design[parts[0]] = parts[1]
    return design


# ════════════════════════════════════════════════════════════════════════════
#  Feature alignment  (sink sequences → source ASV IDs)
# ════════════════════════════════════════════════════════════════════════════

def _is_sequence_index(index: pd.Index) -> bool:
    return bool(re.match(r"^[ACGTacgt]{50,}$", str(index[0])))


def _vsearch_map(query_fa: str, db_fa: str, identity: float,
                 threads: int = 16) -> dict:
    uc = query_fa + ".uc"
    subprocess.run(
        ["vsearch", "--usearch_global", query_fa, "--db", db_fa,
         "--id", str(identity), "--uc", uc,
         "--threads", str(threads), "--strand", "plus",
         "--maxaccepts", "1", "--maxrejects", "32", "--quiet"],
        capture_output=True,
    )
    hits = {}
    if os.path.exists(uc):
        with open(uc) as f:
            for line in f:
                p = line.split("\t")
                if p[0] == "H":
                    hits[p[8]] = p[9].strip()
    return hits


def align_features(sink_df, source_df, db_fasta, mode="asv", threads=16):
    """
    Map sink features (sequences or IDs) to source ASV IDs.

    Returns
    -------
    sink_aligned   : DataFrame  rows=sink_samples,   cols=common_asv_ids
    source_aligned : DataFrame  rows=source_samples, cols=common_asv_ids
    n_matched      : (matched_features, total_sink_features)
    """
    sink_seqbased = _is_sequence_index(sink_df.index)

    with tempfile.TemporaryDirectory() as tmp:
        if sink_seqbased:
            # Write query FASTA with short IDs
            sq_map = {}          # sq_id → original sequence string
            q_fa = os.path.join(tmp, "query.fasta")
            with open(q_fa, "w") as f:
                for i, seq in enumerate(sink_df.index):
                    sid = f"sq{i}"
                    sq_map[sid] = str(seq)
                    f.write(f">{sid}\n{seq}\n")

            # Write db FASTA (only ASVs present in source_df)
            db_fa = os.path.join(tmp, "db.fasta")
            with open(db_fa, "w") as f:
                for asv_id, seq in db_fasta.items():
                    if asv_id in source_df.index:
                        f.write(f">{asv_id}\n{seq}\n")

            identity = 1.0 if mode == "asv" else 0.99
            raw_hits = _vsearch_map(q_fa, db_fa, identity, threads)

            # feature_map: original_seq → source_asv_id
            feature_map = {sq_map[sq_id]: asv_id
                           for sq_id, asv_id in raw_hits.items()}
        else:
            # Both tables already use ASV IDs — direct intersection
            common_ids = sink_df.index.intersection(source_df.index)
            feature_map = {fid: fid for fid in common_ids}

    if not feature_map:
        raise ValueError("No features matched between sink and source tables.")

    # ── Build aligned DataFrames ───────────────────────────────────────────
    matched_source_ids = sorted(set(feature_map.values()))

    # source_aligned: rows=samples, cols=asv_ids
    source_aligned = source_df.loc[
        source_df.index.isin(matched_source_ids)
    ].T  # transpose: samples × features

    # sink_aligned: remap sink feature names → source ASV IDs, then transpose
    new_index = [feature_map.get(str(f)) for f in sink_df.index]
    sink_remapped = sink_df.copy()
    sink_remapped.index = new_index
    sink_remapped = sink_remapped[sink_remapped.index.notna()]
    sink_remapped = sink_remapped.groupby(sink_remapped.index).sum()
    sink_aligned = sink_remapped.T  # samples × features

    common_cols = source_aligned.columns.intersection(sink_aligned.columns)
    source_aligned = source_aligned[common_cols].fillna(0).astype(int)
    sink_aligned   = sink_aligned[common_cols].fillna(0).astype(int)

    return sink_aligned, source_aligned, (len(common_cols), len(sink_df.index))


def collapse_by_group(source_aligned: pd.DataFrame, design: dict) -> pd.DataFrame:
    """Sum source counts within each group → one row per group."""
    cols = source_aligned.columns
    groups = defaultdict(lambda: np.zeros(len(cols), dtype=int))
    for sample in source_aligned.index:
        grp = design.get(sample)
        if grp:
            groups[grp] += source_aligned.loc[sample].values.astype(int)
    if not groups:
        raise ValueError("No source samples matched design entries.")
    return pd.DataFrame(
        {g: v for g, v in groups.items()}, index=cols
    ).T  # groups × features


# ════════════════════════════════════════════════════════════════════════════
#  Gibbs via ST conda env
# ════════════════════════════════════════════════════════════════════════════

def run_gibbs_subprocess(source_collapsed, sink_aligned,
                         src_depth, snk_depth,
                         restarts, burnin, draws,
                         status_fn=None):
    """
    Save data, run _run_gibbs.py via ST conda env, return (mp_df, mps_df).
    status_fn: optional callable(str) for progress messages.
    """
    with tempfile.TemporaryDirectory() as tmp:
        src_csv  = os.path.join(tmp, "sources.csv")
        snk_csv  = os.path.join(tmp, "sinks.csv")
        out_dir  = os.path.join(tmp, "results")

        source_collapsed.to_csv(src_csv)
        sink_aligned.to_csv(snk_csv)

        cmd = [
            "conda", "run", "--no-capture-output", "-n", "ST",
            "python3", _GIBBS_SCRIPT,
            src_csv, snk_csv, out_dir,
            str(src_depth), str(snk_depth),
            str(restarts), str(burnin), str(draws),
        ]

        if status_fn:
            status_fn("Running Gibbs sampling (ST conda env)…")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.stdout:
            for line in result.stdout.strip().splitlines():
                if status_fn:
                    status_fn(line)

        if result.returncode != 0:
            raise RuntimeError(
                f"Gibbs sampler failed:\n{result.stdout}\n{result.stderr}"
            )

        mp  = pd.read_csv(os.path.join(out_dir, "proportions.csv"), index_col=0)
        mps = pd.read_csv(os.path.join(out_dir, "stds.csv"),        index_col=0)

    return mp, mps


# ════════════════════════════════════════════════════════════════════════════
#  Figure
# ════════════════════════════════════════════════════════════════════════════

def make_figure(mp: pd.DataFrame, title: str = "Source Contribution") -> plt.Figure:
    groups = list(mp.columns)
    colors = [GROUP_COLORS.get(g, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
              for i, g in enumerate(groups)]

    samples = list(mp.index)
    x = np.arange(len(samples))
    width = 0.65

    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 1.5), 5.5))
    bottoms = np.zeros(len(samples))

    for gi, (grp, color) in enumerate(zip(groups, colors)):
        vals = mp[grp].values.clip(0)   # clip tiny negatives from MCMC noise
        ax.bar(x, vals, width, bottom=bottoms, color=color,
               label=grp, edgecolor="white", linewidth=0.5)
        for si, (val, bot) in enumerate(zip(vals, bottoms)):
            if val > 0.05:
                ax.text(si, bot + val / 2, f"{val:.0%}",
                        ha="center", va="center", fontsize=8,
                        color="white", fontweight="bold")
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=30, ha="right", fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Proportion", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_patches = [mpatches.Patch(color=c, label=g)
                      for g, c in zip(groups, colors)]
    ax.legend(handles=legend_patches, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), fontsize=9, frameon=False)

    plt.tight_layout()
    return fig


def _fig_to_png(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return buf.read()


def _df_to_csv_bytes(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_csv(buf)
    buf.seek(0)
    return buf.read()


# ════════════════════════════════════════════════════════════════════════════
#  Streamlit UI
# ════════════════════════════════════════════════════════════════════════════

def main():
    st.set_page_config(page_title="SourceTracker", layout="wide")
    st.title("Microbial Source Tracking")
    st.caption("SourceTracker2 Gibbs sampling · feature alignment via vsearch")

    if "st_results" not in st.session_state:
        st.session_state.st_results = None

    # ── 1. Source DB info (fixed) ─────────────────────────────────────────────
    for path, label in [(_SOURCE_TABLE, "db_table.csv.gz"),
                        (_SOURCE_FASTA, "db.fasta"),
                        (_SOURCE_DESIGN, "MST.design")]:
        if not os.path.exists(path):
            st.error(f"Required source file not found: `{path}`")
            return

    # ── 2. Sink table upload ──────────────────────────────────────────────────
    st.subheader("1  Upload sink (environmental) table")
    sink_file = st.file_uploader(
        "Sink table (CSV.gz) — rows = features (sequences or ASV IDs), cols = samples",
        type=["gz"], key="sink",
    )

    if not sink_file:
        st.info(
            "Upload your environmental ASV table (CSV.gz).  \n"
            "Expected format: first column = sequences or ASV IDs, remaining columns = samples.  \n\n"
            "Source reference database (fixed):  \n"
            f"• `db_table.csv.gz` · `db.fasta` · `MST.design`"
        )
        if st.session_state.st_results:
            _show_results(st.session_state.st_results)
        return

    # ── 3. Load all files ─────────────────────────────────────────────────────
    try:
        with st.spinner("Loading files…"):
            sink_df   = _load_csv_gz(sink_file)
            source_df = _load_csv_gz(_SOURCE_TABLE)
            db_fasta  = _load_fasta(_SOURCE_FASTA)
            design    = _load_design(_SOURCE_DESIGN)
    except Exception as e:
        st.error(f"Error loading files: {e}")
        return

    # ── 4. Summary ────────────────────────────────────────────────────────────
    st.subheader("2  Data summary")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Sink samples",   sink_df.shape[1])
    m2.metric("Sink features",  sink_df.shape[0])
    m3.metric("Source samples", source_df.shape[1])
    m4.metric("Source ASVs",    source_df.shape[0])

    groups = sorted(set(design.values()))
    st.markdown(f"**Source groups ({len(groups)}):** {', '.join(groups)}")

    # ── 4. Source group selection ─────────────────────────────────────────────
    st.subheader("3  Select source groups")
    st.caption("Uncheck groups you want to exclude from the analysis.")

    n_cols = min(len(groups), 5)
    check_cols = st.columns(n_cols)
    selected_groups = []
    for i, grp in enumerate(groups):
        color = GROUP_COLORS.get(grp, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
        label = f":{grp}"   # plain label; colour swatch added via markdown below
        checked = check_cols[i % n_cols].checkbox(
            grp, value=True, key=f"grp_{grp}"
        )
        if checked:
            selected_groups.append(grp)

    if not selected_groups:
        st.warning("Select at least one source group.")
        return

    # Filter design to selected groups only
    design_filtered = {s: g for s, g in design.items() if g in selected_groups}

    # ── 5. Parameters ─────────────────────────────────────────────────────────
    st.subheader("4  Parameters")
    pa, pb, pc = st.columns(3)
    with pa:
        mode = st.radio(
            "Feature matching",
            ["ASV  (100% exact)", "OTU  (99% via vsearch)"],
            help="ASV: exact sequence matches only.\nOTU: 99% identity clusters bridge cross-run noise."
        )
        mode_key = "asv" if mode.startswith("ASV") else "otu"

    with pb:
        sink_depths = sink_df.sum(axis=0)
        # Use 50% of the original min depth as default — feature alignment
        # reduces per-sample counts to only matched features, so the
        # post-alignment depth will be lower than the original table total.
        # _run_gibbs.py will further auto-cap if still too high.
        auto_snk = max(500, int(sink_depths.min() * 0.5))
        snk_depth = st.number_input(
            "Sink rarefaction depth", value=auto_snk,
            min_value=100, step=100,
            help=(
                "Applied after feature alignment, so effective counts per sample "
                "are lower than in the original table. The runner will auto-cap "
                "this to the available minimum if it is still too high."
            )
        )
        src_depth = st.number_input(
            "Source rarefaction depth", value=1000,
            min_value=100, step=100,
            help="Applied after collapsing source samples by group."
        )

    with pc:
        restarts = st.number_input("MCMC restarts",        value=10, min_value=1, max_value=100)
        burnin   = st.number_input("Burnin iterations",    value=100, min_value=10, max_value=1000)
        draws    = st.number_input("Draws per restart",    value=1,  min_value=1, max_value=10)

    threads = min(16, os.cpu_count() or 4)
    st.caption(f"Running with {len(selected_groups)} source group(s): {', '.join(selected_groups)}")

    if not st.button("Run SourceTracker", type="primary"):
        if st.session_state.st_results:
            _show_results(st.session_state.st_results)
        return

    # ── 5. Run pipeline ───────────────────────────────────────────────────────
    st.session_state.st_results = None
    progress = st.progress(0)
    status   = st.empty()

    try:
        status.text("Aligning features (vsearch)…")
        sink_al, src_al, (n_m, n_t) = align_features(
            sink_df, source_df, db_fasta, mode=mode_key, threads=threads
        )
        progress.progress(0.25)
        pct = 100 * n_m / max(n_t, 1)
        st.info(
            f"**Feature alignment:** {n_m} / {n_t} sink features matched "
            f"({pct:.1f}%) — unmatched reads will count as Unknown."
        )

        status.text("Collapsing source samples by group…")
        src_col = collapse_by_group(src_al, design_filtered)
        progress.progress(0.35)
        st.caption(
            f"Source groups after collapse: "
            + ", ".join(f"{g} ({int(src_col.loc[g].sum()):,} reads)"
                        for g in src_col.index)
        )

        status.text("Running Gibbs sampling…")
        log_area = st.empty()
        log_lines = []

        def _status(msg):
            log_lines.append(msg)
            log_area.code("\n".join(log_lines[-10:]))

        mp, mps = run_gibbs_subprocess(
            src_col, sink_al,
            src_depth=int(src_depth),
            snk_depth=int(snk_depth),
            restarts=int(restarts),
            burnin=int(burnin),
            draws=int(draws),
            status_fn=_status,
        )
        progress.progress(0.90)

        status.text("Building figure…")
        fig = make_figure(mp, title="Source Contribution (SourceTracker2 Gibbs)")
        progress.progress(1.0)
        status.text("Done.")

        st.session_state.st_results = {
            "mp":  mp, "mps": mps, "fig": fig,
            "fig_bytes": _fig_to_png(fig),
            "pct": (mp * 100).round(2),
        }

    except Exception as e:
        st.error(f"Pipeline failed: {e}")
        st.exception(e)
        return

    _show_results(st.session_state.st_results)


def _show_results(results):
    st.subheader("Results")
    st.pyplot(results["fig"], use_container_width=True)

    mp, mps, pct = results["mp"], results["mps"], results["pct"]

    # Combined mean ± std table
    display = mp.copy().applymap(lambda v: f"{v:.4f}")
    for col in mp.columns:
        if col in mps.columns:
            display[col] = (mp[col].map(lambda v: f"{v:.4f}") + " ± " +
                            mps[col].map(lambda v: f"{v:.4f}"))
    st.markdown("**Mixing proportions (mean ± std)**")
    st.dataframe(display, use_container_width=True)

    # Downloads
    st.subheader("Download")
    d1, d2, d3, d4 = st.columns(4)
    with d1:
        st.download_button("Figure (PNG)", results["fig_bytes"],
                           "sourcetracker_results.png", "image/png")
    with d2:
        st.download_button("Proportions (CSV)", _df_to_csv_bytes(mp),
                           "sourcetracker_proportions.csv", "text/csv")
    with d3:
        st.download_button("Contributions % (CSV)", _df_to_csv_bytes(pct),
                           "sourcetracker_contributions_pct.csv", "text/csv")
    with d4:
        st.download_button("Std dev (CSV)", _df_to_csv_bytes(mps),
                           "sourcetracker_stds.csv", "text/csv")

    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("sourcetracker_results.png",             results["fig_bytes"])
        zf.writestr("sourcetracker_proportions.csv",         _df_to_csv_bytes(mp))
        zf.writestr("sourcetracker_contributions_pct.csv",   _df_to_csv_bytes(pct))
        zf.writestr("sourcetracker_stds.csv",                _df_to_csv_bytes(mps))
    zip_buf.seek(0)
    st.download_button("All results (.zip)", zip_buf.read(),
                       "sourcetracker_results.zip", "application/zip")


if __name__ == "__main__":
    main()
