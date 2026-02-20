#!/usr/bin/env python3
"""Streamlit web UI for DADA2 paired-end analysis.

Accepts trimmed files from app_qc_ui.py (e.g. Sample_v4_250bp_trimmed_R1.fastq.gz).
Auto-detects region/length from filename, groups samples accordingly, applies
per-group DADA2 pipelines, and merges results into a single ASV table.
"""

import gzip
import io
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from collections import Counter

import pandas as pd
import streamlit as st

# ── Truncation defaults (trunc_r1, trunc_r2) keyed by (region_tag, read_length) ──
TRUNC_DEFAULTS = {
    ("v34", 250): (0,   220),
    ("v34", 300): (120, 180),
    ("v4",  250): (180, 120),
    ("v4",  300): (180, 120),
    ("v45", 300): (200, 100),
}


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _autodetect_threads():
    cpu = os.cpu_count() or 1
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    ram_gb = int(line.split()[1]) / (1024 ** 2)
                    break
            else:
                ram_gb = 4.0
    except Exception:
        ram_gb = 4.0
    return min(cpu, 8)


def get_average_length(file_obj, limit=5000):
    """Compute average read length from an uploaded file object."""
    total_len = 0
    count = 0
    try:
        file_obj.seek(0)
        with gzip.open(file_obj, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total_len += len(line.strip())
                    count += 1
                    if count >= limit:
                        break
    except Exception:
        return 0.0
    return total_len / count if count else 0.0


def _detect_region_and_length(sample_name):
    """Detect (tag, read_length) from a trimmed sample name.

    Expects names like: Sample_v4_250bp_trimmed
    """
    name_lower = sample_name.lower()
    for tag in ("v34", "v45", "v4"):
        if f"_{tag}_" in name_lower:
            m = re.search(rf"_{tag}_(\d+)bp", name_lower)
            if m:
                return tag, int(m.group(1))
            return tag, None
    return None, None


def _extract_original_name(sample_name):
    """Strip region/length suffix from a trimmed sample name.

    'MySample_v4_250bp_trimmed' -> 'MySample'
    """
    for tag in ("v34", "v45", "v4"):
        m = re.search(rf"_{tag}_\d+bp", sample_name, re.IGNORECASE)
        if m:
            return sample_name[: m.start()]
    return sample_name


def _pair_files(uploaded_files):
    """Group uploaded files into samples using a greedy regex.

    Supports: _R1/_R2 and _1/_2 naming conventions.
    Returns (samples_dict, unmatched_list).
    """
    samples = {}
    unmatched = []
    for f in uploaded_files:
        match = re.search(r"(.+)_(R?[12])(?:[_.]|\.fastq).*\.gz$", f.name)
        if match:
            s_name = match.group(1)
            raw = match.group(2)
            r_type = raw if raw.startswith("R") else "R" + raw
            if s_name not in samples:
                samples[s_name] = {"R1": None, "R2": None}
            samples[s_name][r_type] = f
        else:
            unmatched.append(f.name)
    return samples, unmatched


# ═══════════════════════════════════════════════════════════════════════════
#  R script generators
# ═══════════════════════════════════════════════════════════════════════════

def _write_filter_script(temp_dir, r_fnFs, r_fnRs, r_sample_names,
                          trunc_r1, trunc_r2, threads):
    script = f"""
library(dada2)
path <- "{temp_dir}"

fnFs <- c({', '.join(r_fnFs)})
fnRs <- c({', '.join(r_fnRs)})
sample.names <- c({', '.join(r_sample_names)})

names(fnFs) <- sample.names
names(fnRs) <- sample.names

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c({trunc_r1},{trunc_r2}),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread={threads})

if (sum(out[, 'reads.out']) == 0) {{
    stop("No reads passed filtering. Adjust truncation lengths.")
}}

saveRDS(filtFs, file.path(path, "filtFs.rds"))
saveRDS(filtRs, file.path(path, "filtRs.rds"))
"""
    path = os.path.join(temp_dir, "1_filter.R")
    with open(path, "w") as f:
        f.write(script)
    return path


def _write_learn_errors_script(temp_dir, threads):
    script = f"""
library(dada2)
path <- "{temp_dir}"

filtFs <- readRDS(file.path(path, "filtFs.rds"))
filtRs <- readRDS(file.path(path, "filtRs.rds"))

errF <- learnErrors(filtFs, multithread={threads})
errR <- learnErrors(filtRs, multithread={threads})

saveRDS(errF, file.path(path, "errF.rds"))
saveRDS(errR, file.path(path, "errR.rds"))
"""
    path = os.path.join(temp_dir, "2_learn_errors.R")
    with open(path, "w") as f:
        f.write(script)
    return path


def _write_dada_script(temp_dir, threads):
    script = f"""
library(dada2)
path <- "{temp_dir}"

filtFs <- readRDS(file.path(path, "filtFs.rds"))
filtRs <- readRDS(file.path(path, "filtRs.rds"))
errF   <- readRDS(file.path(path, "errF.rds"))
errR   <- readRDS(file.path(path, "errR.rds"))

dadaFs <- dada(filtFs, err=errF, multithread={threads})
dadaRs <- dada(filtRs, err=errR, multithread={threads})

saveRDS(dadaFs, file.path(path, "dadaFs.rds"))
saveRDS(dadaRs, file.path(path, "dadaRs.rds"))
"""
    path = os.path.join(temp_dir, "3_dada.R")
    with open(path, "w") as f:
        f.write(script)
    return path


def _write_merge_table_script(temp_dir):
    script = f"""
library(dada2)
path <- "{temp_dir}"

filtFs <- readRDS(file.path(path, "filtFs.rds"))
filtRs <- readRDS(file.path(path, "filtRs.rds"))
dadaFs <- readRDS(file.path(path, "dadaFs.rds"))
dadaRs <- readRDS(file.path(path, "dadaRs.rds"))

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
if (is.data.frame(mergers)) {{
    sname <- names(filtFs)
    if (is.null(sname)) sname <- "sample"
    mergers <- setNames(list(mergers), sname[1])
}}
seqtab <- makeSequenceTable(mergers)

total_reads      <- sum(seqtab)
relative_abundance <- colSums(seqtab) / total_reads
num_samples      <- nrow(seqtab)
prevalence       <- colSums(seqtab > 0) / num_samples
to_remove        <- relative_abundance < 0.001 & prevalence < 0.1
seqtab_filtered  <- seqtab[, !to_remove, drop=FALSE]

sample_names <- rownames(seqtab_filtered)
seq_names    <- colnames(seqtab_filtered)
final_table  <- data.frame(Sequence=seq_names, check.names=FALSE)
for (i in seq_along(sample_names)) {{
    final_table[[sample_names[i]]] <- as.integer(seqtab_filtered[i, ])
}}

results_path <- file.path(path, "dada2_results")
if (!dir.exists(results_path)) dir.create(results_path)

write.csv(final_table, gzfile(file.path(results_path, "ASV_table.csv.gz")), row.names=FALSE)
"""
    path = os.path.join(temp_dir, "4_merge_table.R")
    with open(path, "w") as f:
        f.write(script)
    return path


def _run_r(script_path):
    result = subprocess.run(
        ["Rscript", script_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, script_path,
            output=result.stdout, stderr=result.stderr,
        )


# ═══════════════════════════════════════════════════════════════════════════
#  Table merge helper
# ═══════════════════════════════════════════════════════════════════════════

def _merge_asv_tables(table_bytes_list):
    """Merge a list of (bytes) ASV tables into one combined DataFrame."""
    dfs = [pd.read_csv(io.BytesIO(b), compression="gzip") for b in table_bytes_list]
    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on="Sequence", how="outer")
    merged = merged.fillna(0)
    for col in merged.columns:
        if col != "Sequence":
            merged[col] = merged[col].astype(int)
    return merged


def _df_to_gz_bytes(df):
    buf = io.BytesIO()
    df.to_csv(buf, index=False, compression="gzip")
    buf.seek(0)
    return buf.read()


# ═══════════════════════════════════════════════════════════════════════════
#  Main app
# ═══════════════════════════════════════════════════════════════════════════

def main():
    st.set_page_config(page_title="DADA2 Analysis", layout="wide")
    st.title("DADA2 Analysis for Paired-End FASTQ Files")

    # ── Session state ────────────────────────────────────────────────────
    if "results" not in st.session_state:
        st.session_state.results = None  # dict of downloadable bytes

    # ── File upload ──────────────────────────────────────────────────────
    uploaded_files = st.file_uploader(
        "Upload trimmed paired-end FASTQ files (*.fastq.gz)",
        type=["gz"],
        accept_multiple_files=True,
    )

    if not uploaded_files:
        st.info(
            "Upload trimmed paired-end FASTQ files produced by app_qc_ui.py.  \n"
            "Expected naming: `Sample_v4_250bp_trimmed_R1.fastq.gz` (supports `_R1`/`_R2` and `_1`/`_2`)."
        )
        return

    # ── Pair files ───────────────────────────────────────────────────────
    samples, unmatched = _pair_files(uploaded_files)
    if unmatched:
        st.warning(f"Could not parse: {', '.join(unmatched)}")

    valid = {s: f for s, f in samples.items() if f["R1"] and f["R2"]}
    if not valid:
        st.error("No valid paired samples found. Files must match: `*_R1.fastq.gz` / `*_R2.fastq.gz`")
        return

    # ── Detect region/length groups ──────────────────────────────────────
    region_groups = {}   # (tag, read_len) -> [s_name, ...]
    undetected = []
    for s_name in sorted(valid):
        tag, read_len = _detect_region_and_length(s_name)
        if tag is None or read_len is None or (tag, read_len) not in TRUNC_DEFAULTS:
            undetected.append(s_name)
            continue
        key = (tag, read_len)
        region_groups.setdefault(key, []).append(s_name)

    if undetected:
        st.warning(
            f"Could not detect region/length for: {', '.join(undetected)}  \n"
            "Expected filenames like `Sample_v4_250bp_trimmed_R1.fastq.gz`."
        )
    if not region_groups:
        st.error("No samples with recognisable region tags found.")
        return

    # ── Avg read lengths table ───────────────────────────────────────────
    st.subheader("Detected Samples")
    with st.spinner("Calculating average read lengths..."):
        rows = []
        for (tag, read_len), names in sorted(region_groups.items()):
            for s_name in names:
                r1_avg = get_average_length(valid[s_name]["R1"])
                r2_avg = get_average_length(valid[s_name]["R2"])
                rows.append({
                    "Sample": _extract_original_name(s_name),
                    "Region": tag,
                    "Read length": f"{read_len} bp",
                    "R1 avg len": f"{r1_avg:.1f}",
                    "R2 avg len": f"{r2_avg:.1f}",
                })
    st.dataframe(pd.DataFrame(rows), use_container_width=True)

    # ── Parameters ───────────────────────────────────────────────────────
    st.subheader("Parameters")
    _auto_threads = _autodetect_threads()
    threads = st.number_input(
        "CPU threads", value=_auto_threads, min_value=1, max_value=os.cpu_count() or 1
    )

    # Per-group truncation overrides
    trunc_overrides = {}
    if len(region_groups) > 0:
        st.markdown("**Truncation lengths** (auto-filled from region defaults):")
        cols = st.columns(min(len(region_groups), 3))
        for i, (key, names) in enumerate(sorted(region_groups.items())):
            tag, read_len = key
            default_r1, default_r2 = TRUNC_DEFAULTS[key]
            with cols[i % len(cols)]:
                st.markdown(f"**{tag} / {read_len} bp** ({len(names)} sample(s))")
                r1 = st.number_input(f"Trunc R1 [{tag}/{read_len}bp]", value=default_r1,
                                     min_value=0, key=f"tr1_{tag}_{read_len}")
                r2 = st.number_input(f"Trunc R2 [{tag}/{read_len}bp]", value=default_r2,
                                     min_value=0, key=f"tr2_{tag}_{read_len}")
            trunc_overrides[key] = (r1, r2)

    if not st.button("Run DADA2", type="primary"):
        # Show previous results if available
        if st.session_state.results:
            _show_downloads(st.session_state.results)
        return

    # ── Run DADA2 per group ──────────────────────────────────────────────
    st.session_state.results = None
    all_table_bytes = []
    zip_contents = {}   # filename -> bytes
    total_groups = len(region_groups)

    for g_idx, ((tag, read_len), sample_names) in enumerate(sorted(region_groups.items())):
        trunc_r1, trunc_r2 = trunc_overrides[(tag, read_len)]
        group_label = f"{tag} / {read_len} bp"
        st.markdown(f"---\n**Group {g_idx+1}/{total_groups}: {group_label}** "
                    f"({len(sample_names)} sample(s))")
        progress = st.progress(0)
        status = st.empty()

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # Save files to temp dir
                file_map = {}
                for s_name in sample_names:
                    r1_path = os.path.join(temp_dir, f"{s_name}_R1.fastq.gz")
                    r2_path = os.path.join(temp_dir, f"{s_name}_R2.fastq.gz")
                    valid[s_name]["R1"].seek(0)
                    valid[s_name]["R2"].seek(0)
                    with open(r1_path, "wb") as f:
                        f.write(valid[s_name]["R1"].getvalue())
                    with open(r2_path, "wb") as f:
                        f.write(valid[s_name]["R2"].getvalue())
                    file_map[s_name] = (r1_path, r2_path)

                clean_names = [_extract_original_name(s) for s in sample_names]
                r_fnFs = [f"'{file_map[s][0]}'" for s in sample_names]
                r_fnRs = [f"'{file_map[s][1]}'" for s in sample_names]
                r_sample_names = [f"'{n}'" for n in clean_names]

                steps = [
                    ("Step 1/4: Filtering and trimming...",
                     lambda: _write_filter_script(temp_dir, r_fnFs, r_fnRs, r_sample_names,
                                                  trunc_r1, trunc_r2, threads)),
                    ("Step 2/4: Learning error rates...",
                     lambda: _write_learn_errors_script(temp_dir, threads)),
                    ("Step 3/4: Applying DADA algorithm...",
                     lambda: _write_dada_script(temp_dir, threads)),
                    ("Step 4/4: Merging pairs and building table...",
                     lambda: _write_merge_table_script(temp_dir)),
                ]

                for step_i, (label, gen_fn) in enumerate(steps):
                    status.text(label)
                    _run_r(gen_fn())
                    progress.progress((step_i + 1) / len(steps))

                status.text(f"{group_label} — done.")

                # Collect results
                asv_path = os.path.join(temp_dir, "dada2_results", "ASV_table.csv.gz")
                with open(asv_path, "rb") as f:
                    table_bytes = f.read()
                all_table_bytes.append(table_bytes)
                zip_contents[f"ASV_table_{tag}_{read_len}bp.csv.gz"] = table_bytes

                for err_file in ("errF.rds", "errR.rds"):
                    src = os.path.join(temp_dir, err_file)
                    if os.path.exists(src):
                        stem = err_file.replace(".rds", "")
                        with open(src, "rb") as f:
                            zip_contents[f"{stem}_{tag}_{read_len}bp.rds"] = f.read()

        except subprocess.CalledProcessError as e:
            st.error(f"R script failed for {group_label}:")
            st.code(e.stderr, language="r")
            return
        except Exception as e:
            st.error(f"Unexpected error for {group_label}:")
            st.exception(e)
            return

    # ── Merge tables if multiple groups ──────────────────────────────────
    if len(all_table_bytes) > 1:
        merged_df = _merge_asv_tables(all_table_bytes)
        merged_bytes = _df_to_gz_bytes(merged_df)
        zip_contents["ASV_table_merged.csv.gz"] = merged_bytes
        st.success(
            f"All {total_groups} groups processed. "
            f"Merged table: {len(merged_df)} ASVs, "
            f"{len(merged_df.columns) - 1} sample(s)."
        )
    else:
        merged_bytes = all_table_bytes[0]
        st.success("Analysis complete.")

    # Build zip
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fname, data in zip_contents.items():
            zf.writestr(fname, data)
    zip_buf.seek(0)

    st.session_state.results = {
        "merged_bytes": merged_bytes,
        "zip_bytes": zip_buf.read(),
        "has_multiple": len(all_table_bytes) > 1,
    }
    _show_downloads(st.session_state.results)


def _show_downloads(results):
    st.subheader("Download Results")
    col1, col2 = st.columns(2)
    with col1:
        label = "Download Merged ASV Table" if results["has_multiple"] else "Download ASV Table"
        st.download_button(
            label=label,
            data=results["merged_bytes"],
            file_name="ASV_table_merged.csv.gz" if results["has_multiple"] else "ASV_table.csv.gz",
            mime="application/gzip",
        )
    with col2:
        st.download_button(
            label="Download All Results (.zip)",
            data=results["zip_bytes"],
            file_name="dada2_results.zip",
            mime="application/zip",
        )


if __name__ == "__main__":
    main()
