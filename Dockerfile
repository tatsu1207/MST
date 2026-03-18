# ============================================================================
# MST-Pipeline — Dockerfile
# Builds both conda environments (mst + ST) into a single image.
#
# Build:   docker build -t mst-pipeline .
# Run:     docker compose up -d
# Open:    http://localhost:8050
# ============================================================================

FROM condaforge/mambaforge:latest

LABEL maintainer="MST Pipeline"
LABEL description="Microbial Source Tracking web platform for 16S rRNA amplicon sequencing"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# ── System dependencies ─────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    wget \
    procps \
    lsof \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ── Conda environment 1: mst (Python 3.11 + R/DADA2 + CLI tools) ──────────
RUN mamba create -n mst -c conda-forge python=3.11 -y && \
    mamba clean -afy

# Install R, DADA2, and bioinformatics CLI tools
RUN mamba install -n mst --override-channels -c conda-forge -c bioconda \
    r-base bioconductor-dada2 r-optparse r-jsonlite r-data.table \
    vsearch bbmap cutadapt \
    pandas numpy scipy biopython h5py \
    "tbb=2020.*" -y && \
    mamba clean -afy

# Install Python packages via pip
RUN conda run -n mst pip install --no-cache-dir \
    fastapi \
    "uvicorn[standard]" \
    python-multipart \
    dash \
    dash-bootstrap-components \
    dash-uploader \
    plotly \
    sqlalchemy \
    biom-format \
    fpdf2

# ── Conda environment 2: ST (Python 3.9 + sourcetracker) ──────────────────
RUN mamba create -n ST --override-channels -c conda-forge -c bioconda \
    python=3.9 pandas numpy sourcetracker -y && \
    mamba clean -afy

# ── Copy application code ─────────────────────────────────────────────────
COPY app/ /app/app/
COPY r_scripts/ /app/r_scripts/
COPY config/ /app/config/
COPY DB/ /app/DB/
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# ── Environment variables ─────────────────────────────────────────────────
ENV PORT=8050
ENV CONDA_BASE=/opt/conda
# Store database inside the data volume so it persists
ENV DATABASE_PATH=/app/data/mst.db

EXPOSE 8050

ENTRYPOINT ["/app/docker-entrypoint.sh"]
