#!/usr/bin/env python3
"""Pathogenic bacteria detection from 16S ASV table.

Classifies ASV sequences against the SILVA database using vsearch and
reports relative abundance of pathogenic genera curated from the
WHO Priority Pathogens List (2017) and other authoritative sources.
"""

import gzip, io, os, re, subprocess, tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

_HERE  = os.path.dirname(os.path.abspath(__file__))
_SILVA = os.path.join(_HERE, "..", "config", "silva_nr99_v138.1_train_set_v4.fa.gz")

# ── Pathogen list ─────────────────────────────────────────────────────────────
# Sources: WHO Priority Pathogens 2017, CDC ESKAPE, ECDC, waterborne pathogen lists
PATHOGENS = {
    # WHO Critical Priority (carbapenem-resistant / ESBL threats)
    "Acinetobacter":  "critical",
    "Pseudomonas":    "critical",
    "Klebsiella":     "critical",
    "Escherichia":    "critical",
    "Enterobacter":   "critical",
    "Serratia":       "critical",
    "Proteus":        "critical",
    "Morganella":     "critical",
    "Providencia":    "critical",
    # WHO High Priority
    "Enterococcus":   "high",
    "Staphylococcus": "high",
    "Helicobacter":   "high",
    "Campylobacter":  "high",
    "Salmonella":     "high",
    "Neisseria":      "high",
    # WHO Medium Priority
    "Streptococcus":  "medium",
    "Haemophilus":    "medium",
    "Shigella":       "medium",
    # Other clinically / environmentally important pathogens
    "Clostridioides": "other",
    "Clostridium":    "other",
    "Vibrio":         "other",
    "Listeria":       "other",
    "Yersinia":       "other",
    "Legionella":     "other",
    "Brucella":       "other",
    "Francisella":    "other",
    "Bordetella":     "other",
    "Mycobacterium":  "other",
    "Arcobacter":     "other",
    "Aliarcobacter":  "other",
}

_PRIORITY_ORDER  = ["critical", "high", "medium", "other"]
_PRIORITY_LABELS = {
    "critical": "WHO Critical",
    "high":     "WHO High",
    "medium":   "WHO Medium",
    "other":    "Other Notable",
}
# Each priority group gets its own colour family; genera get distinct shades
_CMAPS = {
    "critical": plt.cm.Reds,
    "high":     plt.cm.Oranges,
    "medium":   plt.cm.YlOrBr,
    "other":    plt.cm.Blues,
}


# ── Data helpers ──────────────────────────────────────────────────────────────

def _load_csv_gz(file_obj):
    if hasattr(file_obj, "seek"):
        file_obj.seek(0)
        data = file_obj.read()
    else:
        with open(file_obj, "rb") as f:
            data = f.read()
    return pd.read_csv(io.BytesIO(data), compression="gzip", index_col=0)


def _is_sequence_index(index):
    return bool(re.match(r"^[ACGTacgt]{50,}$", str(index[0])))


# ── SILVA ─────────────────────────────────────────────────────────────────────

@st.cache_resource(show_spinner="Decompressing SILVA database …")
def _get_silva_path():
    """Decompress SILVA gz once; cache the path for the app's lifetime."""
    tmp = tempfile.mktemp(suffix=".fasta", prefix="silva_mst_")
    with gzip.open(_SILVA, "rt") as gz, open(tmp, "w") as out:
        out.write(gz.read())
    return tmp


# ── Classification ────────────────────────────────────────────────────────────

def _classify(sequences, silva_path, identity, threads):
    """
    Run vsearch --usearch_global against SILVA.
    SILVA headers are 'Bacteria;Phylum;...;Genus;' — the last non-empty
    field is the genus.
    Returns {sequence: genus_str | None}.
    """
    with tempfile.TemporaryDirectory() as tmp:
        # Query FASTA
        q_fa = os.path.join(tmp, "query.fasta")
        id_to_seq = {}
        with open(q_fa, "w") as fh:
            for i, seq in enumerate(sequences):
                sid = f"sq{i}"
                id_to_seq[sid] = seq
                fh.write(f">{sid}\n{seq.upper()}\n")

        uc = os.path.join(tmp, "out.uc")
        subprocess.run(
            ["vsearch", "--usearch_global", q_fa,
             "--db", silva_path,
             "--id", str(identity),
             "--uc", uc,
             "--threads", str(threads),
             "--strand", "plus",
             "--maxaccepts", "1",
             "--maxrejects", "32",
             "--quiet"],
            capture_output=True,
        )

        genus_map = {}   # sid → genus
        if os.path.exists(uc):
            with open(uc) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 10 or parts[0] != "H":
                        continue
                    sid, tax_str = parts[8], parts[9]
                    taxa = [t.strip() for t in tax_str.split(";") if t.strip()]
                    genus_map[sid] = taxa[-1] if taxa else None

    return {id_to_seq[sid]: g for sid, g in genus_map.items()}


# ── Pathogen detection ────────────────────────────────────────────────────────

def _detect(asv_df, seq_to_genus):
    """
    Group ASV counts by pathogenic genus, compute relative abundance.
    Returns (count_df, ra_df) — both indexed samples × genera.
    """
    genera = pd.Series(
        [seq_to_genus.get(s) for s in asv_df.index], index=asv_df.index
    )
    mask = genera.isin(PATHOGENS)
    if not mask.any():
        return pd.DataFrame(), pd.DataFrame()

    grouped = asv_df.loc[mask].copy()
    grouped.index = genera[mask]
    grouped = grouped.groupby(grouped.index).sum()   # genera × samples

    total = asv_df.sum(axis=0)
    ra = grouped.div(total, axis=1) * 100

    return grouped.T, ra.T   # samples × genera


# ── Figure ────────────────────────────────────────────────────────────────────

def _genus_colors(genera):
    """Assign shades within each priority colour family."""
    by_priority = {}
    for g in genera:
        by_priority.setdefault(PATHOGENS.get(g, "other"), []).append(g)
    colors = {}
    for p, gens in by_priority.items():
        n = len(gens)
        for i, g in enumerate(gens):
            colors[g] = _CMAPS[p](0.40 + 0.50 * i / max(n - 1, 1))
    return colors


def _make_figure(ra_df):
    genera = sorted(
        ra_df.columns,
        key=lambda g: (_PRIORITY_ORDER.index(PATHOGENS.get(g, "other")), g),
    )
    ra_sorted = ra_df[genera]
    samples   = ra_sorted.index.tolist()
    colors    = _genus_colors(genera)

    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.7 + 3), 6))

    bottom = np.zeros(len(samples))
    for genus in genera:
        vals = ra_sorted[genus].values
        ax.bar(samples, vals, bottom=bottom,
               color=colors[genus], width=0.7, label=genus)
        bottom += vals

    ax.set_xlabel("Sample")
    ax.set_ylabel("Relative Abundance (% of total reads)")
    ax.set_title("Pathogenic Bacteria — Relative Abundance")
    ax.tick_params(axis="x", rotation=45)
    ax.set_ylim(0, max(bottom.max() * 1.15, 0.01))

    # Legend grouped by priority
    handles = []
    seen_priorities = {PATHOGENS.get(g, "other") for g in genera}
    for p in _PRIORITY_ORDER:
        if p not in seen_priorities:
            continue
        # Section header patch (light, narrow)
        handles.append(mpatches.Patch(
            facecolor=_CMAPS[p](0.25), edgecolor="none",
            label=f"▌ {_PRIORITY_LABELS[p]}",
        ))
        for g in genera:
            if PATHOGENS.get(g, "other") == p:
                handles.append(mpatches.Patch(color=colors[g], label=f"  {g}"))

    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, frameon=False)
    plt.tight_layout()
    return fig


# ── Streamlit UI ──────────────────────────────────────────────────────────────

def main():
    st.header("Pathogenic Bacteria Detection")
    st.markdown(
        "Upload a 16S **ASV table** (CSV.gz, sequences as first column) to identify "
        "pathogenic bacteria. ASVs are classified against **SILVA v138** using vsearch. "
        "Pathogens are curated from the "
        "[WHO Priority Pathogens List (2017)](https://www.who.int/publications/i/item/WHO-EMP-IAU-2017.12), "
        "CDC ESKAPE, and environmental pathogen sources."
    )

    with st.expander("View pathogen list"):
        rows = [
            {"Genus": g, "Priority": _PRIORITY_LABELS[p]}
            for g, p in sorted(
                PATHOGENS.items(),
                key=lambda x: (_PRIORITY_ORDER.index(x[1]), x[0]),
            )
        ]
        st.dataframe(
            pd.DataFrame(rows).set_index("Genus"),
            use_container_width=True,
        )

    st.divider()

    col1, col2, col3 = st.columns([3, 1, 1])
    with col1:
        uploaded = st.file_uploader("ASV table (CSV.gz)", type=["gz"])
    with col2:
        identity = st.slider("SILVA identity (%)", 90, 100, 97) / 100
    with col3:
        threads = st.number_input("Threads", 1, 64, 4)

    if not uploaded:
        st.info("Upload an ASV table to begin.")
        return

    # Load
    try:
        asv_df = _load_csv_gz(uploaded)
    except Exception as e:
        st.error(f"Failed to load table: {e}")
        return

    if asv_df.empty:
        st.error("Uploaded table is empty.")
        return

    if not _is_sequence_index(asv_df.index):
        st.error(
            "Row IDs do not look like DNA sequences. "
            "Please upload a DADA2 ASV table with sequences as the first column."
        )
        return

    st.success(
        f"Loaded: **{asv_df.shape[0]:,}** ASVs × **{asv_df.shape[1]}** samples"
    )

    if not os.path.exists(_SILVA):
        st.error(f"SILVA database not found: `{_SILVA}`")
        return

    silva_path = _get_silva_path()

    # Classify
    with st.spinner(
        f"Classifying {asv_df.shape[0]:,} ASVs against SILVA "
        f"(≥{identity*100:.0f}% identity) …"
    ):
        seq_to_genus = _classify(
            list(asv_df.index), silva_path, identity, int(threads)
        )

    n_classified = sum(1 for g in seq_to_genus.values() if g)
    st.write(
        f"Classified **{n_classified:,} / {asv_df.shape[0]:,}** ASVs at genus level"
    )

    # Detect
    count_df, ra_df = _detect(asv_df, seq_to_genus)

    if ra_df.empty:
        st.warning("No pathogenic bacteria detected in this dataset.")
        return

    st.write(
        f"Detected **{ra_df.shape[1]}** pathogenic genera · "
        f"max relative abundance **{ra_df.values.max():.3f}%**"
    )

    # Figure
    fig = _make_figure(ra_df)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    st.image(buf)
    plt.close(fig)
    st.download_button(
        "Download figure (PNG)", buf,
        file_name="pathogen_abundance.png", mime="image/png",
    )

    # RA table
    st.subheader("Relative Abundance (% of total reads per sample)")
    ra_display = ra_df.copy()
    ra_display.index.name = "Sample"
    st.dataframe(
        ra_display.round(4).style.background_gradient(cmap="Reds", axis=None),
        use_container_width=True,
    )
    st.download_button(
        "Download RA table (CSV)",
        ra_display.to_csv().encode(),
        file_name="pathogen_ra.csv",
        mime="text/csv",
    )

    # Detection summary
    st.subheader("Detection Summary")
    summary = []
    for genus in sorted(
        ra_df.columns,
        key=lambda g: (_PRIORITY_ORDER.index(PATHOGENS.get(g, "other")), g),
    ):
        p = PATHOGENS.get(genus, "other")
        summary.append({
            "Genus":             genus,
            "Priority":          _PRIORITY_LABELS[p],
            "Max RA (%)":        round(float(ra_df[genus].max()), 4),
            "Mean RA (%)":       round(float(ra_df[genus].mean()), 4),
            "Samples detected":  int((ra_df[genus] > 0).sum()),
        })
    st.dataframe(
        pd.DataFrame(summary).set_index("Genus"),
        use_container_width=True,
    )


if __name__ == "__main__":
    st.set_page_config(page_title="Pathogen Detection", layout="wide")
    main()
