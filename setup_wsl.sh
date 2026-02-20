#!/usr/bin/env bash
# =============================================================================
# setup_wsl.sh  —  MST web pipeline conda setup for Windows WSL
# =============================================================================
#
# Creates a single conda environment named ST that contains everything:
#
#   ST
#     • Python 3.9   : streamlit, pandas, numpy, matplotlib, biopython,
#                      sourcetracker (Gibbs sampler)
#     • CLI tools    : vsearch, bbmap (bbmap.sh + reformat.sh), cutadapt
#     • R + DADA2    : r-base, bioconductor-dada2
#
#   The app calls  conda run --no-capture-output -n ST python3 _run_gibbs.py
#   so the environment must be named exactly "ST".
#
# Usage:
#   bash setup_wsl.sh
#
# Requirements:
#   • Miniconda or Anaconda installed
#     https://docs.conda.io/en/latest/miniconda.html
#   • (optional but recommended) mamba for faster dependency solving
#     conda install -n base -c conda-forge mamba
# =============================================================================

set -euo pipefail

ENV_ST="ST"

# ── Terminal colours ──────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[ OK ]${RESET}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
die()     { echo -e "${RED}[ERR ]${RESET}  $*" >&2; exit 1; }
header()  { echo -e "\n${BOLD}$*${RESET}"; }

# ── Detect conda / mamba — install mamba if missing ──────────────────────────
if command -v mamba &>/dev/null; then
    CM="mamba"
    info "Using mamba (fast solver)."
elif command -v conda &>/dev/null; then
    info "mamba not found — installing into base environment…"
    conda install -y -n base -c conda-forge mamba
    # Reload shell functions so the newly installed mamba is on PATH
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    if command -v mamba &>/dev/null; then
        CM="mamba"
        success "mamba installed successfully."
    else
        warn "mamba install succeeded but is not on PATH yet — falling back to conda."
        CM="conda"
    fi
else
    die "conda not found.\nInstall Miniconda: https://docs.conda.io/en/latest/miniconda.html"
fi

# ── Helper: remove env if it already exists ───────────────────────────────────
_env_exists() {
    conda env list | awk '{print $1}' | grep -qx "$1"
}

_maybe_remove() {
    local env="$1"
    if _env_exists "$env"; then
        warn "Environment '$env' already exists."
        read -rp "  Remove and recreate? [y/N] " ans
        if [[ "$ans" =~ ^[Yy]$ ]]; then
            conda env remove -y -n "$env"
            success "Removed '$env'."
            return 0
        else
            info "Keeping existing '$env' — skipping creation."
            return 1
        fi
    fi
    return 0
}

# =============================================================================
#  STEP 1 — Build the ST environment
# =============================================================================
header "══════════════════════════════════════════════════════"
header " Step 1/2 — Building environment: ${ENV_ST}"
header "══════════════════════════════════════════════════════"

if _maybe_remove "$ENV_ST"; then

    # ── 1a. Python + core data-science + sourcetracker ───────────────────────
    # Python 3.9: required for sourcetracker compatibility.
    info "Creating '${ENV_ST}' with Python 3.9 + core packages + sourcetracker…"
    $CM create -y -n "$ENV_ST" \
        -c conda-forge -c bioconda -c defaults \
        python=3.9 \
        streamlit \
        pandas \
        numpy \
        matplotlib \
        biopython \
        sourcetracker
    success "Core Python packages + sourcetracker installed."

    # ── 1b. Bioinformatics CLI tools ──────────────────────────────────────────
    # vsearch   — OTU/ASV alignment and pathogen classification
    # bbmap     — read alignment to E. coli reference for V-region detection
    #             (provides bbmap.sh + reformat.sh; Java pulled in automatically)
    # cutadapt  — primer / adapter trimming
    info "Installing bioinformatics CLI tools (vsearch, bbmap, cutadapt)…"
    $CM install -y -n "$ENV_ST" \
        -c conda-forge -c bioconda \
        vsearch \
        bbmap \
        cutadapt
    success "CLI tools installed."

    # ── 1c. R + DADA2 (Bioconductor) ─────────────────────────────────────────
    # bioconductor-dada2 pulls in r-base automatically.
    # This step can take several minutes due to the R dependency tree.
    info "Installing R + bioconductor-dada2 (this may take several minutes)…"
    $CM install -y -n "$ENV_ST" \
        -c conda-forge -c bioconda \
        r-base \
        bioconductor-dada2
    success "R + DADA2 installed."

    # ── 1d. Pin TBB to 2020.x ────────────────────────────────────────────────
    # dada2 is compiled against the old tbb::task API removed in TBB 2021+.
    # Downgrading TBB (and the matching r-rcppparallel) prevents the
    # "undefined symbol: _ZTIN3tbb4taskE" runtime error on library load.
    info "Pinning TBB to 2020.x for dada2 compatibility…"
    $CM install -y -n "$ENV_ST" \
        -c conda-forge \
        "tbb=2020.*"
    success "TBB pinned."

fi  # end _maybe_remove ST

# =============================================================================
#  STEP 2 — Verify
# =============================================================================
header "══════════════════════════════════════════════════════"
header " Step 2/2 — Verification"
header "══════════════════════════════════════════════════════"

# ── Python packages ───────────────────────────────────────────────────────────
info "Checking Python packages in '${ENV_ST}'…"
conda run -n "$ENV_ST" python - <<'PYCHECK'
import importlib, sys
failures = []
checks = [
    ("streamlit",     "streamlit"),
    ("pandas",        "pandas"),
    ("numpy",         "numpy"),
    ("matplotlib",    "matplotlib"),
    ("biopython",     "Bio"),
    ("sourcetracker", "sourcetracker"),
]
for label, module in checks:
    try:
        m = importlib.import_module(module)
        ver = getattr(m, "__version__", "ok")
        print(f"  \u2713 {label} {ver}")
    except ImportError:
        print(f"  \u2717 {label}  MISSING")
        failures.append(label)
if failures:
    sys.exit(1)
PYCHECK
success "Python packages OK."

# ── sourcetracker Gibbs imports ───────────────────────────────────────────────
info "Checking sourcetracker Gibbs imports…"
conda run -n "$ENV_ST" python - <<'STCHECK'
import sys
try:
    from sourcetracker._sourcetracker import gibbs, subsample_dataframe
    print("  \u2713 gibbs + subsample_dataframe importable")
except ImportError as e:
    print(f"  \u2717 {e}")
    sys.exit(1)
STCHECK
success "sourcetracker Gibbs OK."

# ── CLI tools ─────────────────────────────────────────────────────────────────
info "Checking CLI tools in '${ENV_ST}'…"
conda run -n "$ENV_ST" bash <<'CLICHECK'
ok=0
check() {
    if command -v "$1" &>/dev/null; then
        echo "  ✓ $1"
    else
        echo "  ✗ $1  NOT FOUND"
        ok=1
    fi
}
check bbmap.sh
check reformat.sh
check cutadapt
check vsearch
check Rscript
exit $ok
CLICHECK
success "CLI tools OK."

# ── R dada2 package ───────────────────────────────────────────────────────────
info "Checking R dada2 package…"
conda run -n "$ENV_ST" Rscript - <<'RCHECK'
if (!requireNamespace("dada2", quietly = TRUE)) {
    cat("  \u2717 dada2  NOT FOUND\n")
    quit(status = 1)
}
cat(sprintf("  \u2713 dada2 %s\n", as.character(packageVersion("dada2"))))
RCHECK
success "dada2 R package OK."

# =============================================================================
#  Done — usage instructions
# =============================================================================
echo
echo -e "${GREEN}${BOLD}╔══════════════════════════════════════════════════════════════╗"
echo -e "║                    Setup complete!                          ║"
echo -e "╚══════════════════════════════════════════════════════════════╝${RESET}"
echo
echo -e "  ${BOLD}Run the full pipeline (QC → DADA2 → SourceTracker → Pathogen):${RESET}"
echo -e "    conda run -n ${ENV_ST} streamlit run scripts/app2.py"
echo
echo -e "  ${BOLD}Multi-tab app (all tools in one window):${RESET}"
echo -e "    conda run -n ${ENV_ST} streamlit run scripts/app.py"
echo
echo -e "  ${BOLD}Individual tools:${RESET}"
echo -e "    conda run -n ${ENV_ST} streamlit run scripts/app_qc.py"
echo -e "    conda run -n ${ENV_ST} streamlit run scripts/app_dada2.py"
echo -e "    conda run -n ${ENV_ST} streamlit run scripts/app_sourcetracker_ui.py"
echo -e "    conda run -n ${ENV_ST} streamlit run scripts/app_pathogen_ui.py"
echo
echo -e "  ${BOLD}WSL tip:${RESET} Streamlit may not open a browser automatically."
echo -e "  Navigate to ${CYAN}http://localhost:8501${RESET} in your Windows browser."
echo
echo -e "  ${BOLD}Activate interactively:${RESET}"
echo -e "    conda activate ${ENV_ST}"
echo
