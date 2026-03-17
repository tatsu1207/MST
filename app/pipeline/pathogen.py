"""
MST-Pipeline — WHO pathogen dictionary and detection logic.

Supports genus-level detection (short reads) and species-level refinement
when species assignments are available (e.g., PacBio long reads).
"""
import pandas as pd

# -- Pathogen dictionary (genus level) ----------------------------------------

PATHOGENS = {
    # WHO Critical Priority
    "Acinetobacter": "critical",
    "Pseudomonas": "critical",
    "Klebsiella": "critical",
    "Escherichia": "critical",
    "Enterobacter": "critical",
    "Serratia": "critical",
    "Proteus": "critical",
    "Morganella": "critical",
    "Providencia": "critical",
    # WHO High Priority
    "Enterococcus": "high",
    "Staphylococcus": "high",
    "Helicobacter": "high",
    "Campylobacter": "high",
    "Salmonella": "high",
    "Neisseria": "high",
    # WHO Medium Priority
    "Streptococcus": "medium",
    "Haemophilus": "medium",
    "Shigella": "medium",
    # Other notable pathogens
    "Clostridioides": "other",
    "Clostridium": "other",
    "Vibrio": "other",
    "Listeria": "other",
    "Yersinia": "other",
    "Legionella": "other",
    "Brucella": "other",
    "Francisella": "other",
    "Bordetella": "other",
    "Mycobacterium": "other",
    "Arcobacter": "other",
    "Aliarcobacter": "other",
    # Waterborne & fecal-associated pathogens (MST-relevant)
    "Aeromonas": "other",
    "Plesiomonas": "other",
    "Leptospira": "other",
    "Bacteroides": "other",
    "Treponema": "other",
    "Prevotella": "other",
    "Fusobacterium": "other",
    "Veillonella": "other",
    "Edwardsiella": "other",
    "Citrobacter": "other",
    "Hafnia": "other",
    "Cronobacter": "other",
}

# -- Species-level known pathogens (used when species assignment is available) -
# Key: "Genus species" (lowercase species). Priority overrides genus default
# when a species match is found.

PATHOGENIC_SPECIES = {
    # Critical
    "Acinetobacter baumannii": "critical",
    "Pseudomonas aeruginosa": "critical",
    "Klebsiella pneumoniae": "critical",
    "Escherichia coli": "critical",
    "Enterobacter cloacae": "critical",
    # High
    "Enterococcus faecium": "high",
    "Enterococcus faecalis": "high",
    "Staphylococcus aureus": "high",
    "Helicobacter pylori": "high",
    "Campylobacter jejuni": "high",
    "Campylobacter coli": "high",
    "Salmonella enterica": "high",
    "Neisseria meningitidis": "high",
    "Neisseria gonorrhoeae": "high",
    # Medium
    "Streptococcus pneumoniae": "medium",
    "Streptococcus pyogenes": "medium",
    "Streptococcus agalactiae": "medium",
    "Haemophilus influenzae": "medium",
    # Other
    "Clostridioides difficile": "other",
    "Clostridium perfringens": "other",
    "Clostridium botulinum": "other",
    "Vibrio cholerae": "other",
    "Vibrio parahaemolyticus": "other",
    "Vibrio vulnificus": "other",
    "Listeria monocytogenes": "other",
    "Yersinia pestis": "other",
    "Yersinia enterocolitica": "other",
    "Legionella pneumophila": "other",
    "Brucella melitensis": "other",
    "Brucella abortus": "other",
    "Francisella tularensis": "other",
    "Bordetella pertussis": "other",
    "Mycobacterium tuberculosis": "other",
    "Leptospira interrogans": "other",
    "Edwardsiella tarda": "other",
}

# -- Species known to be non-pathogenic within pathogenic genera ---------------
# ASVs assigned to these species are excluded from pathogen results.

NON_PATHOGENIC_SPECIES = {
    "Escherichia fergusonii",
    "Escherichia albertii",  # debated, but rarely clinical
    "Pseudomonas fluorescens",
    "Pseudomonas putida",
    "Pseudomonas stutzeri",
    "Pseudomonas syringae",  # plant pathogen, not human
    "Acinetobacter lwoffii",
    "Acinetobacter johnsonii",
    "Acinetobacter radioresistens",
    "Staphylococcus epidermidis",
    "Staphylococcus saprophyticus",
    "Staphylococcus hominis",
    "Staphylococcus warneri",
    "Staphylococcus capitis",
    "Enterococcus casseliflavus",
    "Enterococcus gallinarum",
    "Enterococcus mundtii",
    "Streptococcus thermophilus",
    "Streptococcus salivarius",
    "Streptococcus oralis",
    "Streptococcus mitis",
    "Streptococcus sanguinis",
    "Streptococcus gordonii",
    "Vibrio fischeri",
    "Vibrio natriegens",
    "Clostridium butyricum",
    "Clostridium beijerinckii",
    "Neisseria lactamica",
    "Neisseria flavescens",
    "Neisseria sicca",
    "Mycobacterium smegmatis",
    "Mycobacterium gordonae",
    "Mycobacterium vaccae",
    "Campylobacter fetus",  # primarily veterinary
    "Aeromonas media",
    "Bacteroides vulgatus",
    "Bacteroides uniformis",
    "Bacteroides ovatus",
    "Bacteroides thetaiotaomicron",
    "Prevotella copri",
    "Fusobacterium varium",
    "Veillonella parvula",
    "Veillonella dispar",
    "Veillonella atypica",
}

_PRIORITY_ORDER = ["critical", "high", "medium", "other"]
_PRIORITY_LABELS = {
    "critical": "WHO Critical",
    "high": "WHO High",
    "medium": "WHO Medium",
    "other": "Other Notable",
}


def _resolve_label_and_priority(genus, species):
    """Determine pathogen label and priority from genus + species info.

    Returns:
        (label, priority) or (None, None) if not a pathogen.
        label is "Genus species" when species is known, else "Genus".
    """
    if not isinstance(genus, str):
        return None, None

    # Normalize genus (SILVA may use "Clostridium sensu stricto 1" etc.)
    genus_base = genus.split()[0]

    if genus_base not in PATHOGENS:
        return None, None

    # If species is available, check species-level lists
    if isinstance(species, str) and species and str(species).lower() != "nan":
        species_clean = species.strip()
        full_name = f"{genus_base} {species_clean}"

        # Exclude known non-pathogens
        if full_name in NON_PATHOGENIC_SPECIES:
            return None, None

        # Use species-level priority if defined
        if full_name in PATHOGENIC_SPECIES:
            return full_name, PATHOGENIC_SPECIES[full_name]

        # Species assigned but not in either list: report as "Genus species"
        # with genus-level priority (conservative — flag it)
        return full_name, PATHOGENS[genus_base]

    # No species info — genus-level only
    return genus_base, PATHOGENS[genus_base]


def detect(asv_df: pd.DataFrame, seq_to_genus: dict,
           seq_to_species: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Group ASV counts by pathogenic taxa, compute relative abundance.

    Args:
        asv_df: features x samples DataFrame (index = sequences)
        seq_to_genus: {sequence: genus_str | None} from taxonomy
        seq_to_species: {sequence: species_str | None} from taxonomy (optional)

    Returns:
        (count_df, ra_df) — both indexed samples x taxa labels.
        Labels are "Genus species" when species data is available,
        otherwise "Genus". Returns empty DataFrames if no pathogens found.
    """
    if seq_to_species is None:
        seq_to_species = {}

    labels = {}  # sequence -> label
    priorities = {}  # label -> priority

    for seq in asv_df.index:
        genus = seq_to_genus.get(seq)
        species = seq_to_species.get(seq)
        label, priority = _resolve_label_and_priority(genus, species)
        if label:
            labels[seq] = label
            priorities[label] = priority

    if not labels:
        return pd.DataFrame(), pd.DataFrame()

    label_series = pd.Series(labels)
    mask = asv_df.index.isin(label_series.index)

    grouped = asv_df.loc[mask].copy()
    grouped.index = grouped.index.map(label_series)
    grouped = grouped.groupby(grouped.index).sum()

    total = asv_df.sum(axis=0)
    ra = grouped.div(total, axis=1) * 100

    # Store priorities for downstream use
    detect._priorities = priorities

    return grouped.T, ra.T  # samples x taxa


def get_priority(label: str) -> str:
    """Get the priority level for a pathogen label."""
    # Check species-level first, then genus
    if label in PATHOGENIC_SPECIES:
        return PATHOGENIC_SPECIES[label]
    genus = label.split()[0]
    return PATHOGENS.get(genus, "other")


def genus_colors(genera: list[str]) -> dict[str, str]:
    """Assign colors to taxa labels based on their priority group."""
    _COLOR_FAMILIES = {
        "critical": ["#ef5350", "#e53935", "#c62828", "#b71c1c", "#d32f2f",
                      "#f44336", "#e57373", "#ff8a80", "#ff5252"],
        "high": ["#ff9800", "#f57c00", "#ef6c00", "#e65100", "#fb8c00", "#ffa726"],
        "medium": ["#fdd835", "#f9a825", "#f57f17"],
        "other": ["#42a5f5", "#1e88e5", "#1565c0", "#0d47a1", "#2196f3",
                   "#64b5f6", "#90caf9", "#bbdefb", "#1976d2", "#0277bd", "#01579b"],
    }

    by_priority = {}
    for g in genera:
        p = get_priority(g)
        by_priority.setdefault(p, []).append(g)

    colors = {}
    for p, gens in by_priority.items():
        family = _COLOR_FAMILIES.get(p, _COLOR_FAMILIES["other"])
        for i, g in enumerate(gens):
            colors[g] = family[i % len(family)]
    return colors


def build_summary(ra_df: pd.DataFrame) -> list[dict]:
    """Build a detection summary table from relative abundance data."""
    summary = []
    for taxon in sorted(
        ra_df.columns,
        key=lambda g: (_PRIORITY_ORDER.index(get_priority(g)), g),
    ):
        p = get_priority(taxon)
        genus = taxon.split()[0]
        species = taxon.split(maxsplit=1)[1] if " " in taxon else ""
        summary.append({
            "Taxon": taxon,
            "Genus": genus,
            "Species": species,
            "Priority": _PRIORITY_LABELS[p],
            "Max RA (%)": round(float(ra_df[taxon].max()), 4),
            "Mean RA (%)": round(float(ra_df[taxon].mean()), 4),
            "Samples detected": int((ra_df[taxon] > 0).sum()),
        })
    return summary
