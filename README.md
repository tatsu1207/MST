# MST Web — Microbial Source Tracking Web UI

Streamlit web application for 16S amplicon-based microbial source tracking and pathogenic bacteria detection.

## Features

| Tab | Description |
|-----|-------------|
| **QC & Extraction** | Primer trimming and V-region extraction (cutadapt) |
| **DADA2** | ASV inference from paired-end FASTQ files |
| **SourceTracker** | Source contribution estimation (SourceTracker2 + Gibbs sampling) |
| **Pathogen Detection** | Relative abundance of WHO priority pathogens (SILVA classification) |

`app2.py` runs the full pipeline end-to-end from raw FASTQ files in a single page.

## Windows Setup (WSL)

This pipeline runs on Linux. On Windows, use WSL (Windows Subsystem for Linux).

**1. Install WSL** (run in PowerShell or Command Prompt as Administrator):

```powershell
wsl --install
```

This installs WSL 2 and Ubuntu by default. Restart your computer when prompted.

**2. Open Ubuntu** from the Start menu and complete the first-time setup (create a Linux username and password).

**3. Install Miniconda inside WSL:**

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompts, accept the license, and allow it to initialize conda. Then restart your terminal (close and reopen Ubuntu).

**4. Clone this repository and run the setup script:**

```bash
git clone https://github.com/your-org/MST-web.git
cd MST-web
bash setup_wsl.sh
```

After setup, open the app with:

```bash
conda activate ST
streamlit run scripts/app2.py
```

Then navigate to `http://localhost:8501` in your Windows browser (Streamlit will not open a browser automatically under WSL).

---

## Requirements

### Conda environment

Run the provided setup script to create the `ST` conda environment:

```bash
bash setup_wsl.sh
```

The script installs mamba if not present, creates the `ST` environment, and verifies every dependency. It must be named exactly `ST` because the app calls `conda run -n ST` internally.

To install manually instead (Python 3.9 required for sourcetracker compatibility):

```bash
conda create -n ST python=3.9
conda activate ST
conda install -c conda-forge -c bioconda -c defaults \
    streamlit pandas numpy matplotlib biopython sourcetracker
conda install -c conda-forge -c bioconda \
    vsearch bbmap cutadapt
conda install -c conda-forge -c bioconda \
    r-base bioconductor-dada2
```

### SILVA reference database

The SILVA files are too large for git. Download from the DADA2 reference page and place in `config/`:

```bash
# DADA2-formatted SILVA v138.1 (genus-level, required for pathogen detection)
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz \
     -O config/silva_nr99_v138.1_train_set.fa.gz

# Optional: species-level assignment
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz \
     -O config/silva_nr99_v138.1_wSpecies_train_set.fa.gz
```

## Usage

### End-to-end pipeline UI (`app2.py`)

```bash
conda activate ST
streamlit run scripts/app2.py
```

Upload raw FASTQ.gz files, configure options, and run the full pipeline in one click. Displays source contribution and pathogen abundance graphs when done.

### Multi-tab UI (`app.py`)

```bash
conda activate ST
streamlit run scripts/app.py
```

Opens four tabs: QC, DADA2, SourceTracker, Pathogen Detection. Each tab can be used independently with its own file uploads.

### Individual tools

```bash
conda activate ST
streamlit run scripts/app_qc.py
streamlit run scripts/app_dada2.py
streamlit run scripts/app_sourcetracker_ui.py
streamlit run scripts/app_pathogen_ui.py
```

### WSL (Windows Subsystem for Linux)

Streamlit may not open a browser automatically under WSL. After starting the app, navigate manually to:

```
http://localhost:8501
```

in your Windows browser.

## Database files

| File | Description |
|------|-------------|
| `DB/db.fasta` | 16S ASV reference sequences |
| `DB/db_table.csv.gz` | Source ASV count table |
| `DB/db_table.norm.csv.gz` | GCN-normalized source ASV table |
| `DB/MST.design` | Sample-to-source-group mapping |

## Pathogen list

Pathogens are curated from the [WHO Priority Pathogens List (2017)](https://www.who.int/publications/i/item/WHO-EMP-IAU-2017.12), CDC ESKAPE pathogens, and environmental pathogen sources. Detection is performed by classifying ASV sequences against SILVA using vsearch and matching to known pathogenic genera.
