# MST-Pipeline: Microbial Source Tracking Web Platform

An open-source, web-based platform for end-to-end microbial source tracking and pathogen detection from 16S rRNA amplicon sequencing data.

## Features

- **Automatic sequencing detection** — identifies paired-end/single-end/long-read and variable region (V4, V3-V4, V4-V5, V1-V9)
- **Automated quality optimization** — auto-detects truncation lengths and filtering parameters
- **DADA2 ASV inference** — denoising, merging, and chimera removal
- **Cross-region harmonization** — extracts V4 sub-region from any amplicon type for unified source tracking
- **Source tracking** — Gibbs sampling via SourceTracker2 against a built-in reference database (10 source groups, 2,100 samples)
- **Pathogen detection** — screens 44 WHO-priority pathogenic genera with abundance summaries
- **BIOM file support** — upload pre-processed BIOM files directly, skipping DADA2

## Quick Start

### 1. Installation

**Linux / WSL (Windows):**

```bash
git clone https://github.com/tatsu1207/MST-Pipeline.git
cd MST-Pipeline
bash setup_ubuntu.sh    # or setup_wsl.sh on Windows WSL
```

This creates the required conda environments and downloads the SILVA v138.1 reference database.

**Windows users:** Install WSL first:
```powershell
wsl --install
```
Then open Ubuntu from the Start menu and follow the Linux instructions above.

### 2. Start the server

```bash
conda activate mst_py
uvicorn app.main:app --host 0.0.0.0 --port 8050
```

Open `http://localhost:8050` in your browser.

## Usage

### Option A: From raw FASTQ files

1. **Upload** — Drag and drop your `.fastq.gz` files into the "Upload FASTQ Files" panel, then click **Register Upload**. Both paired-end (`_R1`/`_R2`) and single-end files are accepted. The pipeline automatically detects the sequencing type and variable region.

2. **Select samples** — After registration, uploaded samples appear in a table. Check the boxes next to the samples you want to analyze.

3. **Run DADA2** — Click **Start DADA2**. The pipeline performs primer trimming, quality filtering, ASV inference, and chimera removal. Progress is displayed in real time.

4. **Run source tracking** — Once DADA2 completes, the processed dataset appears in the **Pipeline History** table. Select the dataset by checking its box. Before clicking **Run Analysis**, you can use the **Source Groups** checkboxes to include only the sources relevant to your study (e.g., uncheck seawater and bat for an inland site). Excluding irrelevant sources speeds up the analysis and reduces spurious assignments. Then click **Run Analysis** at the bottom of the page to run SourceTracker2 and pathogen detection.

5. **View results** — Source contribution charts (stacked bar, pie) and pathogen detection tables appear in the results panel.

### Combining samples from different variable regions

Samples sequenced with different primer sets (e.g., V4 and V3-V4) can be analyzed together. Upload and process each batch through DADA2 separately — the pipeline auto-detects the variable region and applies appropriate parameters for each run. At the source tracking step, select **multiple datasets** from the Pipeline History table. The pipeline automatically extracts the V4 sub-region from each dataset based on its detected variable region and merges them into a unified analysis. This allows you to combine data from different sequencing runs, primer sets, or studies in a single MST analysis.

### Option B: From a BIOM file

If you have already processed your data through DADA2 (or another pipeline), you can upload a BIOM file directly:

1. **Upload** — Drag and drop your `.biom` file into the "Upload FASTQ Files" panel, then click **Register Upload**. The BIOM file must contain ASV sequences in the observation metadata (not taxonomy).

2. **Select and run** — The uploaded BIOM file appears in the **Pipeline History** table. Select it and click **Run Analysis** to run source tracking and pathogen detection directly, skipping DADA2.

### Supported input formats

| Format | Naming convention | Notes |
|--------|-------------------|-------|
| Paired-end FASTQ | `sample_R1.fastq.gz` / `sample_R2.fastq.gz` | Also supports `_R1_001` and `_1`/`_2` conventions |
| Single-end FASTQ | `sample.fastq.gz` | |
| Long-read FASTQ | `sample.fastq.gz` | PacBio HiFi auto-detected by read length |
| BIOM (HDF5) | `filename.biom` | Must include ASV sequences in observation metadata |

### Supported variable regions

The pipeline automatically detects the variable region and applies appropriate processing:

| Region | Forward primer | Reverse primer | Expected amplicon |
|--------|---------------|----------------|-------------------|
| V4 | 515F | 806R | ~253 bp |
| V3-V4 | 341F | 806R | ~460 bp |
| V4-V5 | 515F | 926R | ~400 bp |
| V1-V9 | 27F | 1492R | ~1,500 bp |

For non-V4 amplicons, the V4 sub-region is automatically extracted before source tracking.

## Source Reference Database

The built-in source reference database contains 55,958 ASVs from 2,100 samples across 10 source groups:

| Source group | Samples |
|---|---|
| Bat | 515 |
| Human | 347 |
| Pig | 300 |
| Cow | 267 |
| Horse | 242 |
| Groundwater | 168 |
| Chicken | 121 |
| Duck | 74 |
| Seawater | 43 |
| River | 27 |

## Database files

| File | Description |
|------|-------------|
| `DB/db.fasta` | 16S ASV reference sequences |
| `DB/db_table.csv.gz` | Source ASV count table |
| `DB/MST.design` | Sample-to-source-group mapping |

## Requirements

- Linux (Ubuntu 20.04+) or Windows WSL
- Miniconda or Miniforge
- 8 GB RAM minimum (16 GB recommended)
- The setup script creates two conda environments:
  - `mst_py` — Python 3.11, FastAPI, Plotly Dash, DADA2 (R), bioinformatics tools
  - `mst_st2` — Python 3.9, SourceTracker2 (separate due to dependency constraints)

## Citation

If you use MST-Pipeline in your research, please cite:

> Unno, T. (2026). MST-Pipeline: An Integrated Web-Based Platform for Microbial Source Tracking and Pathogen Detection from 16S rRNA Amplicon Sequencing Data.

## License

MIT
