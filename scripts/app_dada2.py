try:
    import streamlit as st
except ImportError:
    st = None
import os
import subprocess
import pandas as pd
from tempfile import NamedTemporaryFile, TemporaryDirectory
import gzip
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

def get_average_length(file):
    total_len = 0
    count = 0
    read_limit = 5000  # Process only the first 5000 reads
    with gzip.open(file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_len += len(record)
            count += 1
            if count >= read_limit:
                break
    return total_len / count if count > 0 else 0

def main():
    st.title("DADA2 Analysis for Paired-End FASTQ Files")

    # Initialize session state
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
        st.session_state.table_data = None
        st.session_state.errF_data = None
        st.session_state.errR_data = None

    st.header("Upload your FASTQ files")
    uploaded_files = st.file_uploader("Upload paired-end FASTQ files (.fastq.gz)", accept_multiple_files=True)

    if uploaded_files:
        st.header("Average Read Lengths")
        with st.spinner("Calculating average read lengths..."):
            grouped_files = {}
            for file in uploaded_files:
                import re as _re
                _m = _re.search(r'_(R[12])\.fastq(\.gz)?$', file.name)
                sample_name = file.name[:_m.start()] if _m else file.name.split('_R')[0]
                r_type = _m.group(1) if _m else ('R1' if file.name.endswith('_1.fastq.gz') else 'R2')
                if sample_name not in grouped_files:
                    grouped_files[sample_name] = {'R1': None, 'R2': None}
                grouped_files[sample_name][r_type] = file
            
            avg_lengths_data = []
            for sample_name, files in grouped_files.items():
                r1_avg_len = "N/A"
                r2_avg_len = "N/A"
                if files['R1']:
                    r1_avg_len = f"{get_average_length(files['R1']):.2f}"
                if files['R2']:
                    r2_avg_len = f"{get_average_length(files['R2']):.2f}"
                avg_lengths_data.append({"Sample Name": sample_name, "R1 Average Length": r1_avg_len, "R2 Average Length": r2_avg_len})
            
            st.table(pd.DataFrame(avg_lengths_data))

        st.header("Define Analysis Parameters")
        trunc_len_r1 = st.number_input("Truncation length for R1", value=240, min_value=1)
        trunc_len_r2 = st.number_input("Truncation length for R2", value=160, min_value=1)
        threads = st.number_input("Number of CPU threads to use", value=2, min_value=1)

        if st.button("Run DADA2 Analysis"):
            st.session_state.analysis_complete = False # Reset on new run
            st.info("Starting DADA2 analysis...")
            with TemporaryDirectory() as temp_dir:
                
                # --- 1. SAVE UPLOADED FILES ---
                samples = {}
                for sample_name, files in grouped_files.items():
                    samples[sample_name] = {'R1': None, 'R2': None}
                    if files['R1']:
                        file_path_r1 = os.path.join(temp_dir, files['R1'].name)
                        files['R1'].seek(0)
                        with open(file_path_r1, "wb") as f:
                            f.write(files['R1'].getbuffer())
                        samples[sample_name]['R1'] = file_path_r1
                    if files['R2']:
                        file_path_r2 = os.path.join(temp_dir, files['R2'].name)
                        files['R2'].seek(0)
                        with open(file_path_r2, "wb") as f:
                            f.write(files['R2'].getbuffer())
                        samples[sample_name]['R2'] = file_path_r2

                try:
                    # --- 2. VALIDATE SAMPLES ---
                    st.info(f"{len(grouped_files)} samples were detected based on file names.")
                    valid_samples = {s: f for s, f in grouped_files.items() if f['R1'] and f['R2']}
                    st.info(f"{len(valid_samples)} samples have both R1 and R2 files and will be processed.")
                    
                    if not valid_samples:
                        st.error("No valid paired-end files found. Please upload R1 and R2 files for at least one sample.")
                        return

                    r_fnFs = [f"'{samples[s]['R1']}'" for s in valid_samples]
                    r_fnRs = [f"'{samples[s]['R2']}'" for s in valid_samples]
                    r_sample_names = [f"'{s}'" for s in valid_samples]

                    # --- 3. RUN DADA2 PIPELINE ---
                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    # Step 1: Filtering and Trimming
                    status_text.text("Step 1/4: Filtering and Trimming...")
                    script_path = generate_filter_script(temp_dir, r_fnFs, r_fnRs, r_sample_names, trunc_len_r1, trunc_len_r2, threads)
                    run_r_script(script_path)
                    progress_bar.progress(25)

                    # Step 2: Learning Errors
                    status_text.text("Step 2/4: Learning Errors...")
                    script_path = generate_learn_errors_script(temp_dir, threads)
                    run_r_script(script_path)
                    progress_bar.progress(50)

                    # Step 3: Dereplication and DADA
                    status_text.text("Step 3/4: Applying DADA Algorithm...")
                    script_path = generate_dada_script(temp_dir, threads)
                    run_r_script(script_path)
                    progress_bar.progress(75)

                    # Step 4: Merging Pairs and Creating Table
                    status_text.text("Step 4/4: Merging Pairs and Finalizing...")
                    script_path = generate_merge_and_table_script(temp_dir)
                    run_r_script(script_path)
                    progress_bar.progress(100)
                    status_text.text("Analysis Complete!")

                    # --- 4. STORE RESULTS IN SESSION STATE ---
                    st.success("DADA2 analysis complete!")
                    results_path = os.path.join(temp_dir, "dada2_results")
                    asv_table_path = os.path.join(results_path, "ASV_table.csv.gz")

                    with open(asv_table_path, "rb") as f:
                        st.session_state.table_data = f.read()
                    with open(os.path.join(temp_dir, "errF.rds"), "rb") as f:
                        st.session_state.errF_data = f.read()
                    with open(os.path.join(temp_dir, "errR.rds"), "rb") as f:
                        st.session_state.errR_data = f.read()

                    st.session_state.analysis_complete = True

                except subprocess.CalledProcessError as e:
                    st.error(f"An error occurred during an R script execution:")
                    st.code(e.stderr, language="bash")
                    st.session_state.analysis_complete = False
                except Exception as e:
                    st.error(f"An unexpected error occurred:")
                    st.exception(e)
                    st.session_state.analysis_complete = False

    # Always show download buttons if analysis is complete
    if st.session_state.analysis_complete:
        st.header("Download Results")
        st.download_button(
            label="Download ASV Abundance Table",
            data=st.session_state.table_data,
            file_name="ASV_table.csv.gz",
            mime="application/gzip"
        )
        st.download_button(
            label="Download Forward Error Rates (errF.rds)",
            data=st.session_state.errF_data,
            file_name="errF.rds",
            mime="application/octet-stream"
        )
        st.download_button(
            label="Download Reverse Error Rates (errR.rds)",
            data=st.session_state.errR_data,
            file_name="errR.rds",
            mime="application/octet-stream"
        )


def run_r_script(script_path):
    """Helper function to execute an R script and handle errors."""
    try:
        subprocess.run(["Rscript", script_path], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        # Re-raise the exception to be caught by the main try-except block
        raise e

def generate_filter_script(temp_dir, r_fnFs, r_fnRs, r_sample_names, trunc_len_r1, trunc_len_r2, threads):
    script = f"""
    library(dada2)
    path <- "{temp_dir}"
    
    fnFs <- c({', '.join(r_fnFs)})
    fnRs <- c({', '.join(r_fnRs)})
    sample.names <- c({', '.join(r_sample_names)})

    if(length(fnFs) != length(fnRs) || length(fnFs) != length(sample.names)) {{
        stop("Mismatched number of forward, reverse, or sample names.")
    }}
    
    names(fnFs) <- sample.names
    names(fnRs) <- sample.names
    
    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c({trunc_len_r1},{trunc_len_r2}),
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread={threads})
    
    if(sum(out[, 'reads.out']) == 0) {{
        stop("No reads passed the filtering step. Adjust truncation lengths or other filtering parameters.")
    }}
    
    saveRDS(filtFs, file.path(path, "filtFs.rds"))
    saveRDS(filtRs, file.path(path, "filtRs.rds"))
    """
    script_path = os.path.join(temp_dir, "1_filter.R")
    with open(script_path, "w") as f:
        f.write(script)
    return script_path

def generate_learn_errors_script(temp_dir, threads):
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
    script_path = os.path.join(temp_dir, "2_learn_errors.R")
    with open(script_path, "w") as f:
        f.write(script)
    return script_path

def generate_dada_script(temp_dir, threads):
    script = f"""
    library(dada2)
    path <- "{temp_dir}"
    
    filtFs <- readRDS(file.path(path, "filtFs.rds"))
    filtRs <- readRDS(file.path(path, "filtRs.rds"))
    errF <- readRDS(file.path(path, "errF.rds"))
    errR <- readRDS(file.path(path, "errR.rds"))
    
    dadaFs <- dada(filtFs, err=errF, multithread={threads})
    dadaRs <- dada(filtRs, err=errR, multithread={threads})
    
    saveRDS(dadaFs, file.path(path, "dadaFs.rds"))
    saveRDS(dadaRs, file.path(path, "dadaRs.rds"))
    """
    script_path = os.path.join(temp_dir, "3_dada.R")
    with open(script_path, "w") as f:
        f.write(script)
    return script_path

def generate_merge_and_table_script(temp_dir):
    script = f"""
    library(dada2)
    path <- "{temp_dir}"

    filtFs <- readRDS(file.path(path, "filtFs.rds"))
    filtRs <- readRDS(file.path(path, "filtRs.rds"))
    dadaFs <- readRDS(file.path(path, "dadaFs.rds"))
    dadaRs <- readRDS(file.path(path, "dadaRs.rds"))

    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
    # Wrap single-sample result in a named list for makeSequenceTable
    if (is.data.frame(mergers)) {{
        sname <- names(filtFs)
        if (is.null(sname)) sname <- "sample"
        mergers <- setNames(list(mergers), sname[1])
    }}
    seqtab <- makeSequenceTable(mergers)

    # --- Filter out rare ASVs ---
    # 1. Calculate Relative Abundance
    total_reads <- sum(seqtab)
    relative_abundance <- colSums(seqtab) / total_reads

    # 2. Calculate Prevalence
    num_samples <- nrow(seqtab)
    prevalence <- colSums(seqtab > 0) / num_samples

    # 3. Identify ASVs to remove (rare AND not prevalent)
    to_remove <- relative_abundance < 0.001 & prevalence < 0.1
    
    # 4. Filter the sequence table (drop=FALSE keeps matrix when 1 sample)
    seqtab_filtered <- seqtab[, !to_remove, drop=FALSE]

    # Build data frame: ASVs as rows, samples as columns
    sample_names <- rownames(seqtab_filtered)
    seq_names <- colnames(seqtab_filtered)
    final_table <- data.frame(Sequence=seq_names, check.names=FALSE)
    for (i in seq_along(sample_names)) {{
        final_table[[sample_names[i]]] <- as.integer(seqtab_filtered[i, ])
    }}
    
    results_path <- file.path(path, "dada2_results")
    if (!dir.exists(results_path)) {{
        dir.create(results_path)
    }}
    
    # Save the final table as a gzipped CSV
    write.csv(final_table, gzfile(file.path(results_path, "ASV_table.csv.gz")), row.names=FALSE)
    """
    script_path = os.path.join(temp_dir, "4_merge_table.R")
    with open(script_path, "w") as f:
        f.write(script)
    return script_path

# (trunc_r1, trunc_r2) keyed by (region_tag, read_length)
TRUNC_DEFAULTS = {
    ("v34", 250): (0, 220),
    ("v34", 300): (120, 180),
    ("v4",  250): (180, 120),
    ("v4",  300): (180, 120),
    ("v45", 300): (200, 100),
}


def _detect_region_and_length(sample_name):
    """Detect region tag and read length from sample name.

    Expects filenames like: Sample_v45_300bp_trimmed_R1.fastq.gz
    Returns (tag, read_length) or (None, None).
    """
    import re
    name_lower = sample_name.lower()
    for tag in ("v34", "v45", "v4"):
        if f"_{tag}_" in name_lower:
            m = re.search(rf"_{tag}_(\d+)bp", name_lower)
            if m:
                return tag, int(m.group(1))
            return tag, None
    return None, None


def _extract_sample_name(full_name):
    """Extract original sample name from trimmed filename.

    'V45_v45_300bp_trimmed' -> 'V45'
    'v4_250_v4_250bp_trimmed' -> 'v4_250'
    """
    import re
    # Strip everything from the region tag onward
    for tag in ("v34", "v45", "v4"):
        m = re.search(rf"_{tag}_\d+bp", full_name, re.IGNORECASE)
        if m:
            return full_name[:m.start()]
    return full_name


def merge_asv_tables(table_paths, output_path):
    """Merge multiple ASV tables into one combined table.

    Each table has a 'Sequence' column and one or more sample columns.
    Tables are concatenated with missing samples filled as 0.
    """
    dfs = []
    for path in table_paths:
        df = pd.read_csv(path, compression="gzip")
        dfs.append(df)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on="Sequence", how="outer")
    merged = merged.fillna(0)

    # Convert count columns back to int
    for col in merged.columns:
        if col != "Sequence":
            merged[col] = merged[col].astype(int)

    merged.to_csv(output_path, index=False, compression="gzip")
    print(f"  Merged {len(table_paths)} tables: "
          f"{len(merged)} ASVs, {len(merged.columns) - 1} sample(s)")


def cli_main():
    """Command-line interface for testing DADA2 pipeline."""
    import argparse
    import shutil

    parser = argparse.ArgumentParser(
        description="DADA2 analysis for paired-end FASTQ files (CLI mode)")
    parser.add_argument("files", nargs="+", help="Paired-end FASTQ files (*_R1.fastq.gz *_R2.fastq.gz)")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory (default: .)")
    parser.add_argument("--trunc-r1", type=int, default=None,
                        help="Override truncation length R1 (default: auto from region)")
    parser.add_argument("--trunc-r2", type=int, default=None,
                        help="Override truncation length R2 (default: auto from region)")
    _auto_threads = min(os.cpu_count() or 1, 8)
    parser.add_argument("--threads", type=int, default=_auto_threads,
                        help=f"CPU threads (default: auto={_auto_threads})")
    parser.add_argument("--dry-run", action="store_true", help="Show what would run without executing")
    args = parser.parse_args()

    # Group files into samples
    grouped = {}
    for fpath in args.files:
        if not os.path.exists(fpath):
            print(f"ERROR: {fpath} not found")
            return 1
        basename = os.path.basename(fpath)
        import re
        m = re.search(r'_(R[12])\.fastq(\.gz)?$', basename)
        sample_name = basename[:m.start()] if m else basename.split('_R')[0]
        r_type = m.group(1) if m else ('R1' if basename.endswith('_1.fastq.gz') else 'R2')
        if sample_name not in grouped:
            grouped[sample_name] = {'R1': None, 'R2': None}
        grouped[sample_name][r_type] = fpath

    valid = {s: f for s, f in grouped.items() if f['R1'] and f['R2']}
    if not valid:
        print("ERROR: No valid paired-end samples found.")
        return 1

    # Detect regions and read lengths from filenames
    sample_info = {}  # s_name -> (tag, read_length, r1_avg, r2_avg)
    region_groups = {}  # (tag, read_length) -> [s_name, ...]
    for s_name in sorted(valid):
        tag, read_len = _detect_region_and_length(s_name)
        if tag is None:
            print(f"ERROR: Cannot detect region from sample name '{s_name}'.")
            print("  Expected filename like: Sample_v4_300bp_trimmed_R1.fastq.gz")
            known_tags = sorted(set(t for t, _ in TRUNC_DEFAULTS))
            print(f"  Known tags: {', '.join(known_tags)}")
            return 1
        if read_len is None:
            print(f"ERROR: Cannot detect read length from sample name '{s_name}'.")
            print("  Expected filename like: Sample_v4_300bp_trimmed_R1.fastq.gz")
            return 1
        key = (tag, read_len)
        if key not in TRUNC_DEFAULTS:
            print(f"ERROR: No truncation defaults for {tag} @ {read_len}bp")
            return 1
        files = valid[s_name]
        r1_avg = get_average_length(files['R1'])
        r2_avg = get_average_length(files['R2'])
        sample_info[s_name] = (tag, read_len, r1_avg, r2_avg)
        if key not in region_groups:
            region_groups[key] = []
        region_groups[key].append(s_name)

    print(f"=== {len(valid)} sample(s) detected ===")
    for s_name in sorted(valid):
        tag, read_len, r1_avg, r2_avg = sample_info[s_name]
        print(f"  {s_name}: R1 avg={r1_avg:.1f} bp, R2 avg={r2_avg:.1f} bp  [{tag} {read_len}bp]")

    # Process each region/length group separately
    for (tag, read_len), sample_names in sorted(region_groups.items()):
        trunc_r1 = args.trunc_r1 if args.trunc_r1 is not None else TRUNC_DEFAULTS[(tag, read_len)][0]
        trunc_r2 = args.trunc_r2 if args.trunc_r2 is not None else TRUNC_DEFAULTS[(tag, read_len)][1]

        print(f"\n{'='*60}")
        print(f"Region: {tag} @ {read_len}bp ({len(sample_names)} sample(s))")
        print(f"Truncation: R1={trunc_r1}, R2={trunc_r2}, threads={args.threads}")

        if args.dry_run:
            print("(dry run — skipping DADA2 execution)")
            continue

        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temp_dir:
            samples = {}
            for s_name in sample_names:
                files = valid[s_name]
                samples[s_name] = {'R1': None, 'R2': None}
                for rt in ['R1', 'R2']:
                    if files[rt]:
                        dest = os.path.join(temp_dir, os.path.basename(files[rt]))
                        shutil.copy2(files[rt], dest)
                        samples[s_name][rt] = dest

            r_fnFs = [f"'{samples[s]['R1']}'" for s in sample_names]
            r_fnRs = [f"'{samples[s]['R2']}'" for s in sample_names]
            clean_names = [_extract_sample_name(s) for s in sample_names]
            r_sample_names = [f"'{n}'" for n in clean_names]

            steps = [
                ("1/4: Filtering and Trimming",
                 lambda: generate_filter_script(temp_dir, r_fnFs, r_fnRs, r_sample_names,
                                                trunc_r1, trunc_r2, args.threads)),
                ("2/4: Learning Errors",
                 lambda: generate_learn_errors_script(temp_dir, args.threads)),
                ("3/4: Applying DADA Algorithm",
                 lambda: generate_dada_script(temp_dir, args.threads)),
                ("4/4: Merging Pairs and Finalizing",
                 lambda: generate_merge_and_table_script(temp_dir)),
            ]

            for label, gen_script in steps:
                print(f"\n  Step {label}...")
                script_path = gen_script()
                try:
                    run_r_script(script_path)
                    print("    OK")
                except subprocess.CalledProcessError as e:
                    print(f"    FAILED:\n{e.stderr}")
                    return 1

            # Copy results out
            src = os.path.join(temp_dir, "dada2_results", "ASV_table.csv.gz")
            if os.path.exists(src):
                os.makedirs(args.outdir, exist_ok=True)
                dest = os.path.join(args.outdir, f"ASV_table_{tag}_{read_len}bp.csv.gz")
                shutil.copy2(src, dest)
                print(f"\n  Output: {dest}")
            else:
                print("ERROR: ASV_table.csv.gz not found in results")
                return 1

            os.makedirs(args.outdir, exist_ok=True)
            for err_file in ("errF.rds", "errR.rds"):
                src_err = os.path.join(temp_dir, err_file)
                if os.path.exists(src_err):
                    stem = err_file.replace(".rds", "")
                    dest_err = os.path.join(args.outdir, f"{stem}_{tag}_{read_len}bp.rds")
                    shutil.copy2(src_err, dest_err)
                    print(f"  Error rates: {dest_err}")

    if args.dry_run:
        print("\n(dry run complete)")
        return 0

    # Merge all individual ASV tables into one
    asv_files = sorted(
        os.path.join(args.outdir, f)
        for f in os.listdir(args.outdir)
        if f.startswith("ASV_table_") and f.endswith(".csv.gz")
        and "merged" not in f
    )
    if len(asv_files) > 1:
        merged_path = os.path.join(args.outdir, "ASV_table_merged.csv.gz")
        merge_asv_tables(asv_files, merged_path)
        print(f"\nMerged table: {merged_path}")

    print("\n=== All regions processed successfully ===")
    return 0


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        sys.exit(cli_main() or 0)
    else:
        main()
