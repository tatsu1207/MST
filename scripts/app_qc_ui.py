#!/usr/bin/env python3
"""Streamlit web UI for paired-end FASTQ QC, variable region ID, and V4 extraction.

Auto-detects V-region, read length, and primer status from the data, then
applies the correct cutadapt command from the 12-case extraction matrix.
"""

import gzip
import io
import multiprocessing
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

ECOLI_REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config", "ecoli.fas")

VREGIONS = [
    ("V1", 69, 99),
    ("V2", 137, 242),
    ("V3", 433, 497),
    ("V4", 576, 682),
    ("V5", 822, 879),
    ("V6", 986, 1043),
    ("V7", 1117, 1173),
    ("V8", 1243, 1294),
    ("V9", 1435, 1465),
]

# Primer constants (IUPAC)
PRIMERS = {
    "515F":    "GTGYCAGCMGCCGCGGTAA",
    "806R":    "GACTACNVGGGTWTCTAATCC",
    "515F_RC": "TTACCGCGGCKGCTGRCAC",
    "806R_RC": "GGATTAGAWACCCBNGTAGTC",
}

VALID_V4_REGIONS = {"V4", "V3-V4", "V4-V5"}
VREGION_TAG = {"V4": "v4", "V3-V4": "v34", "V4-V5": "v45"}

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]",
    "S": "[GC]", "W": "[AT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

# ── 12-case matrix description (for UI) ─────────────────────────────────────
CASE_TABLE = """
| Case | Region | Length | Primers | cutadapt flags |
|------|--------|--------|---------|----------------|
| 1 | V4 | 250 bp | Included | `-g 515F -G 806R` |
| 2 | V4 | 250 bp | Trimmed | *(no action)* |
| 3 | V4 | 300 bp | Included | `-g 515F -G 806R -a 806R_RC -A 515F_RC` |
| 4 | V4 | 300 bp | Trimmed | `-G 806R -a 806R_RC -A 515F_RC` |
| 5 | V3-V4 | 250 bp | Included | `-g 515F --overlap 16 -G 806R` |
| 6 | V3-V4 | 250 bp | Trimmed | `-g 515F --overlap 16` |
| 7 | V3-V4 | 300 bp | Included | `-g 515F --overlap 16 -G 806R -A 515F_RC` |
| 8 | V3-V4 | 300 bp | Trimmed | `-g 515F --overlap 16 -A 515F_RC` |
| 9 | V4-V5 | 250 bp | Included | *Error — insufficient overlap* |
| 10 | V4-V5 | 250 bp | Trimmed | *Error — insufficient overlap* |
| 11 | V4-V5 | 300 bp | Included | `-g 515F -G 806R --overlap 16 -a 806R_RC --minimum-length 150:40` |
| 12 | V4-V5 | 300 bp | Trimmed | `-G 806R --overlap 16 -a 806R_RC --minimum-length 150:40` |
"""


# ═══════════════════════════════════════════════════════════════════════════
#  QC helpers
# ═══════════════════════════════════════════════════════════════════════════

def get_comprehensive_stats(file_path, n_bases=20, limit=50000):
    if not file_path or not os.path.exists(file_path):
        return "N/A", 0, 0, 0, 0, 0, 0

    bp_counts = Counter()
    pos_stats = []
    total_len = 0
    read_count = 0
    total_reads = 0

    try:
        with gzip.open(file_path, "rt") as f:
            for i, line in enumerate(f):
                line_type = i % 4
                if line_type == 1:
                    total_reads += 1
                    if read_count < limit:
                        seq = line.strip()
                        s_len = len(seq)
                        total_len += s_len
                        read_count += 1
                        if s_len >= n_bases:
                            bp_counts[seq[:n_bases]] += 1
                elif line_type == 3 and read_count <= limit:
                    q_str = line.strip()
                    for pos, char in enumerate(q_str):
                        q_score = ord(char) - 33
                        if len(pos_stats) <= pos:
                            pos_stats.append([0, 0])
                        pos_stats[pos][0] += q_score
                        pos_stats[pos][1] += 1
    except Exception:
        return "ERR", 0, 0, 0, 0, 0, 0

    if read_count == 0:
        return "N/A", 0, 0, 0, 0, 0, 0

    avg_len = total_len / read_count
    most_common = bp_counts.most_common(1)[0][0] if bp_counts else "N/A"
    averages = [s[0] / s[1] if s[1] > 0 else 0 for s in pos_stats]

    def find_sustained_drop(avg_list, threshold, window=5):
        for idx in range(len(avg_list) - window):
            if sum(avg_list[idx : idx + window]) / window < threshold:
                return idx
        return len(avg_list)

    q30_drop = int(min(avg_len, find_sustained_drop(averages, 30)))
    q25_drop = int(min(avg_len, find_sustained_drop(averages, 25)))
    q20_drop = int(min(avg_len, find_sustained_drop(averages, 20)))

    return most_common, avg_len, read_count, q30_drop, q25_drop, q20_drop, total_reads


def _parse_sam(sam_file):
    starts, ends = [], []
    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue
            cols = line.strip().split("\t")
            flag = int(cols[1])
            if flag & 4:
                continue
            pos = int(cols[3])
            cigar = cols[5]
            alen = sum(int(l) for l, op in re.findall(r"(\d+)([MDNX=])", cigar))
            starts.append(pos)
            ends.append(pos + alen - 1)
    return starts, ends


def _bbtools_error(stderr_bytes):
    """Extract the first meaningful error line from BBTools stderr bytes."""
    text = stderr_bytes.decode(errors="replace")
    for line in text.splitlines():
        line = line.strip()
        if line and not line.startswith("at ") and "Exception" in line or (
            line and not line.startswith("at ") and line.startswith(("Error", "java.", "Could not", "File not"))
        ):
            return line[:300]
    # fallback: first non-empty line
    for line in text.splitlines():
        if line.strip():
            return line.strip()[:300]
    return text[:300].strip()


def identify_vregion(r1_path, r2_path, ref_path, n_reads=1000, timeout=60, threads=1):
    if not r1_path or not r2_path or not ref_path:
        return "N/A"
    if not os.path.exists(ref_path):
        return "N/A (ref missing)"
    if not shutil.which("bbmap.sh") or not shutil.which("reformat.sh"):
        return "N/A (bbmap not found)"

    tmpdir = tempfile.mkdtemp()
    try:
        ref = subprocess.run(
            [
                "reformat.sh",
                f"in={r1_path}", f"in2={r2_path}",
                f"out={tmpdir}/sub_R1.fastq.gz", f"out2={tmpdir}/sub_R2.fastq.gz",
                f"samplereadstarget={n_reads}", "-Xmx1g",
            ],
            capture_output=True,
            timeout=timeout,
        )
        # If reformat.sh fails (e.g. file has fewer reads than samplereadstarget),
        # fall back to running bbmap directly on the originals with reads=n_reads.
        if ref.returncode != 0:
            bbmap_inputs = [(r1_path, "r1"), (r2_path, "r2")]
            use_reads_limit = True
        else:
            bbmap_inputs = [
                (f"{tmpdir}/sub_R1.fastq.gz", "r1"),
                (f"{tmpdir}/sub_R2.fastq.gz", "r2"),
            ]
            use_reads_limit = False

        for infile, label in bbmap_inputs:
            extra = [f"reads={n_reads}"] if use_reads_limit else []
            bm = subprocess.run(
                [
                    "bbmap.sh", f"ref={ref_path}",
                    f"in={infile}", f"out={tmpdir}/{label}.sam",
                    "nodisk", f"threads={threads}", "ambiguous=best",
                    "minid=0.50", "maxindel=100", "bwr=0.25",
                    "semiperfectmode=f", "strictmaxindel=f", "-Xmx1g",
                ] + extra,
                capture_output=True,
                timeout=timeout,
            )
            if bm.returncode != 0 or not os.path.exists(f"{tmpdir}/{label}.sam"):
                msg = _bbtools_error(bm.stderr)
                return f"N/A (bbmap failed on {label}: {msg})"

        r1_starts, r1_ends = _parse_sam(f"{tmpdir}/r1.sam")
        r2_starts, r2_ends = _parse_sam(f"{tmpdir}/r2.sam")

        n_r1, n_r2 = len(r1_starts), len(r2_starts)
        if n_r1 == 0 and n_r2 == 0:
            return "N/A (no reads mapped to reference)"

        if n_r1 > 0:
            r1_starts.sort()
            amp_start = r1_starts[n_r1 // 2]
        else:
            r2_starts.sort()
            amp_start = r2_starts[n_r2 // 2]

        if n_r2 > 0:
            r2_ends.sort()
            amp_end = r2_ends[n_r2 // 2]
        else:
            r1_ends.sort()
            amp_end = r1_ends[n_r1 // 2]

        covered = []
        for name, v_start, v_end in VREGIONS:
            ov_start = max(amp_start, v_start)
            ov_end = min(amp_end, v_end)
            if ov_end > ov_start:
                cov = (ov_end - ov_start) / (v_end - v_start)
                if cov >= 0.5:
                    covered.append(name)

        if not covered:
            return "Unknown"
        if len(covered) == 1:
            return covered[0]
        return f"{covered[0]}-{covered[-1]}"
    except Exception as e:
        return f"N/A ({e})"
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
#  Auto-detection helpers
# ═══════════════════════════════════════════════════════════════════════════

def _iupac_to_regex(seq):
    return "".join(IUPAC.get(c.upper(), c) for c in seq)


def _prefix_matches_primer(prefix, primer_seq, min_match=15):
    """Check if a read prefix starts with a primer (IUPAC-aware)."""
    if not isinstance(prefix, str) or len(prefix) < min_match:
        return False
    pattern = _iupac_to_regex(primer_seq[:min_match])
    return bool(re.match(pattern, prefix[:min_match], re.IGNORECASE))


def fast_identify_vregion(r1_prefix, r2_prefix):
    """Identify V-region from read prefixes without running bbmap.

    Works when amplification primers are still present in the reads:
      - R1 starts with 515F  + R2 starts with 806R  → V4
      - R1 starts with 515F  + R2 lacks 806R         → V4-V5
      - R1 lacks 515F        + R2 starts with 806R   → V3-V4

    Returns the region string, or None when primers are absent/ambiguous
    (caller should fall back to identify_vregion / bbmap).
    """
    r1_515F = _prefix_matches_primer(r1_prefix, PRIMERS["515F"])
    r2_806R = _prefix_matches_primer(r2_prefix, PRIMERS["806R"])
    if r1_515F and r2_806R:
        return "V4"
    if r1_515F:
        return "V4-V5"
    if r2_806R:
        return "V3-V4"
    return None


def detect_read_length(avg_len):
    """Classify average read length into 250 or 300 bp."""
    return 250 if avg_len <= 275 else 300


def detect_primer_status(r1_prefix, vregion, r2_prefix=None):
    """Detect whether amplification primers are still in the reads.

    For V4/V4-V5: check if R1 starts with 515F.
    For V3-V4: check if R2 starts with 806R (R1 uses a V3 primer, not 515F).
    """
    if vregion in ("V4", "V4-V5"):
        if _prefix_matches_primer(r1_prefix, PRIMERS["515F"]):
            return "Included"
        return "Trimmed"
    # V3-V4: check R2 for 806R
    if r2_prefix and _prefix_matches_primer(r2_prefix, PRIMERS["806R"]):
        return "Included"
    return "Trimmed"


def detect_r2_pad(r2_path, limit=1000):
    """Detect dark-cycle pad at the 5' end of R2 reads.

    Returns the number of leading bases to trim (0 if no pad detected).
    Looks for positions where >50% of reads have N within the first 5 bp.
    """
    if not r2_path or not os.path.exists(r2_path):
        return 0
    n_at_pos = {}
    count = 0
    try:
        with gzip.open(r2_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    seq = line.strip()
                    for pos in range(min(5, len(seq))):
                        if seq[pos] == "N":
                            n_at_pos[pos] = n_at_pos.get(pos, 0) + 1
                    count += 1
                    if count >= limit:
                        break
    except Exception:
        return 0
    if count == 0:
        return 0
    pad = 0
    for pos in range(5):
        if n_at_pos.get(pos, 0) / count > 0.5:
            pad = pos + 1
    return pad


# ═══════════════════════════════════════════════════════════════════════════
#  12-case cutadapt logic
# ═══════════════════════════════════════════════════════════════════════════

def get_cutadapt_args(region, length, primer_status):
    """Return (adapter_args_list, case_number, error_msg).

    adapter_args is a list of cutadapt flag/value pairs, [] for no-action,
    or None for error cases.  Includes --minimum-length per case.
    """
    included = primer_status == "Included"
    minlen_v4  = ["--minimum-length", "150:150"]
    minlen_v34 = ["--minimum-length", "40:150"]

    if region == "V4":
        if length == 250:
            if included:
                return ["-g", PRIMERS["515F"], "-G", PRIMERS["806R"]
                        ] + minlen_v4, 1, None
            return [], 2, None
        # 300 bp
        if included:
            return [
                "-g", PRIMERS["515F"], "-G", PRIMERS["806R"],
                "-a", PRIMERS["806R_RC"], "-A", PRIMERS["515F_RC"],
            ] + minlen_v4, 3, None
        return ["-a", PRIMERS["806R_RC"], "-A", PRIMERS["515F_RC"],
                ] + minlen_v4, 4, None

    if region == "V3-V4":
        if length == 250:
            if included:
                return ["-g", PRIMERS["515F"], "--overlap", "16",
                        "-G", PRIMERS["806R"]] + minlen_v34, 5, None
            return ["-g", PRIMERS["515F"], "--overlap", "16",
                    ] + minlen_v34, 6, None
        # 300 bp
        if included:
            return ["-g", PRIMERS["515F"], "--overlap", "16",
                    "-G", PRIMERS["806R"], "-A", PRIMERS["515F_RC"],
                    ] + minlen_v34, 7, None
        return ["-g", PRIMERS["515F"], "--overlap", "16",
                "-A", PRIMERS["515F_RC"]] + minlen_v34, 8, None

    if region == "V4-V5":
        if length == 250:
            case = 9 if included else 10
            return None, case, (
                "V4-V5 with 250 bp reads: insufficient overlap after V4 "
                "extraction for paired-end merging. Use 300 bp reads."
            )
        # 300 bp
        if included:
            return [
                "-g", PRIMERS["515F"],
                "-G", PRIMERS["806R"],
                "-a", PRIMERS["806R_RC"],
                "--times", "2",
                "--overlap", "16",
                "--discard-untrimmed",
                "--minimum-length", "200:40",
            ], 11, None
        return [
                "-a", PRIMERS["806R_RC"],
                "-A", PRIMERS["806R"],
                "--overlap", "16",
                "--discard-untrimmed",
                "--minimum-length", "200:40",
        ], 12, None

    return None, 0, f"Unknown region: {region}"


def run_cutadapt(args, r1_in, r2_in, r1_out, r2_out, timeout=300, threads=1):
    """Execute cutadapt with the given adapter arguments."""
    cmd = ["cutadapt", "-j", str(threads)] + args + [
        "-o", r1_out, "-p", r2_out,
        r1_in, r2_in,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return result.returncode == 0, result.stderr
    except subprocess.TimeoutExpired:
        return False, f"cutadapt timed out after {timeout}s"


def fastq_avg_len(file_path, limit=50000):
    """Fast average read length from a fastq.gz file."""
    if not file_path or not os.path.exists(file_path):
        return 0.0
    total_len = 0
    count = 0
    try:
        with gzip.open(file_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total_len += len(line.strip())
                    count += 1
                    if count >= limit:
                        break
    except Exception:
        return 0.0
    return total_len / count if count else 0.0


def _format_args(args):
    """Format cutadapt args list into a readable string."""
    if not args:
        return "(no action)"
    parts = []
    i = 0
    while i < len(args):
        flag = args[i]
        if flag.startswith("-") and i + 1 < len(args) and not args[i + 1].startswith("-"):
            parts.append(f"{flag} {_primer_name(args[i + 1])}")
            i += 2
        else:
            parts.append(flag)
            i += 1
    return " ".join(parts)


def _primer_name(seq):
    """Reverse-lookup a primer sequence to its short name."""
    for name, s in PRIMERS.items():
        if seq == s:
            return name
        if seq == f"X{s}":
            return f"X{name}"
    return seq[:12] + "..."


# ═══════════════════════════════════════════════════════════════════════════
#  HTML table builder
# ═══════════════════════════════════════════════════════════════════════════

def build_html_table(rows):
    html = """<style>
    .qc-table { border-collapse: collapse; width: 100%; font-size: 14px; }
    .qc-table th { background-color: #2c3e50; color: white; padding: 8px 12px;
                    text-align: center; border: 1px solid #ddd; }
    .qc-table td { padding: 6px 12px; text-align: center; border: 1px solid #ddd; }
    .qc-table tr:nth-child(even) { background-color: #f8f9fa; }
    .qc-table tr:hover { background-color: #e8f4f8; }
    .qc-table .sample-name { text-align: left; font-weight: bold; }
    .qc-table .vregion { background-color: #eaf7e6; font-weight: bold; }
    .qc-table .vregion-bad { background-color: #fde8e8; font-weight: bold; color: #c0392b; }
    </style>"""

    html += '<table class="qc-table">'
    html += "<thead><tr>"
    html += '<th rowspan="2">Sample</th>'
    html += '<th colspan="6">R1 (Forward)</th>'
    html += '<th colspan="6">R2 (Reverse)</th>'
    html += '<th rowspan="2">V-Region</th>'
    html += '<th rowspan="2">Length</th>'
    html += '<th rowspan="2">Primers</th>'
    html += '<th rowspan="2">Case</th>'
    html += "</tr><tr>"
    for _ in range(2):
        html += "<th>Reads</th><th>Prefix</th><th>Avg Len</th>"
        html += "<th>Q30</th><th>Q25</th><th>Q20</th>"
    html += "</tr></thead><tbody>"

    for row in rows:
        html += "<tr>"
        html += f'<td class="sample-name">{row["sample"]}</td>'
        for direction in ["r1", "r2"]:
            html += f'<td>{row[f"{direction}_reads"]:,}</td>'
            prefix = row[f"{direction}_prefix"]
            html += f"<td><code>{prefix}</code></td>"
            html += f'<td>{row[f"{direction}_avg_len"]:.1f}</td>'
            html += f'<td>{row[f"{direction}_q30"]}</td>'
            html += f'<td>{row[f"{direction}_q25"]}</td>'
            html += f'<td>{row[f"{direction}_q20"]}</td>'
        vr = row["vregion"]
        css_class = "vregion" if vr in VALID_V4_REGIONS else "vregion-bad"
        html += f'<td class="{css_class}">{vr}</td>'
        html += f'<td>{row["det_length"]} bp</td>'
        html += f'<td>{row["det_primers"]}</td>'
        html += f'<td>{row["det_case"]}</td>'
        html += "</tr>"

    html += "</tbody></table>"
    return html


# ═══════════════════════════════════════════════════════════════════════════
#  Resource auto-detection
# ═══════════════════════════════════════════════════════════════════════════

def _autodetect_resources():
    """Return sensible (jobs, threads) defaults based on CPU count and RAM.

    threads: per-tool parallelism (cutadapt -j, bbmap threads=), capped at 4.
    jobs:    parallel samples, limited by both CPU and RAM
             (each job needs ~2 GB for a bbmap JVM instance).
    """
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

    threads = min(4, cpu)
    cpu_jobs = max(1, cpu // threads)
    ram_jobs = max(1, int(ram_gb // 2))
    jobs = min(cpu_jobs, ram_jobs)
    return jobs, threads


# ═══════════════════════════════════════════════════════════════════════════
#  Parallel helpers
# ═══════════════════════════════════════════════════════════════════════════

def _process_qc_sample(s_name, paths, n_bases, read_limit, threads):
    """Run QC + auto-detection for one sample. Returns a row dict."""
    r1_stats = get_comprehensive_stats(paths["R1"], n_bases, read_limit)
    r2_stats = get_comprehensive_stats(paths["R2"], n_bases, read_limit)
    vregion = fast_identify_vregion(r1_stats[0], r2_stats[0])
    if vregion is None:
        vregion = identify_vregion(paths["R1"], paths["R2"], ECOLI_REF, threads=threads)
    det_length = detect_read_length(r1_stats[1])
    det_primers = detect_primer_status(r1_stats[0], vregion, r2_stats[0])
    args, case_num, error = get_cutadapt_args(vregion, det_length, det_primers)
    r2_pad = detect_r2_pad(paths["R2"])
    if r2_pad > 0 and "-G" not in (args or []):
        if args is None:
            args = []
        args = ["-U", str(r2_pad)] + args
    return {
        "sample": s_name,
        "r1_reads": r1_stats[6], "r1_prefix": r1_stats[0],
        "r1_avg_len": r1_stats[1], "r1_q30": r1_stats[3],
        "r1_q25": r1_stats[4], "r1_q20": r1_stats[5],
        "r2_reads": r2_stats[6], "r2_prefix": r2_stats[0],
        "r2_avg_len": r2_stats[1], "r2_q30": r2_stats[3],
        "r2_q25": r2_stats[4], "r2_q20": r2_stats[5],
        "vregion": vregion, "det_length": det_length,
        "det_primers": det_primers, "det_case": case_num,
        "cutadapt_args": args, "cutadapt_error": error,
    }


def _process_trim_sample(row, paths, tmpdir, threads):
    """Run cutadapt for one sample. Returns (s_name, tag, rlen, case_num, r1_out, r2_out, err)."""
    s_name = row["sample"]
    tag = VREGION_TAG[row["vregion"]]
    args = row["cutadapt_args"]
    case_num = row["det_case"]
    rlen = row["det_length"]
    r1_out = os.path.join(tmpdir, f"{s_name}_{tag}_{rlen}bp_trimmed_R1.fastq.gz")
    r2_out = os.path.join(tmpdir, f"{s_name}_{tag}_{rlen}bp_trimmed_R2.fastq.gz")
    if not args:
        shutil.copy2(paths["R1"], r1_out)
        shutil.copy2(paths["R2"], r2_out)
        return s_name, tag, rlen, case_num, r1_out, r2_out, None
    ok, stderr = run_cutadapt(args, paths["R1"], paths["R2"], r1_out, r2_out, threads=threads)
    if ok:
        return s_name, tag, rlen, case_num, r1_out, r2_out, None
    return s_name, tag, rlen, case_num, None, None, stderr


# ═══════════════════════════════════════════════════════════════════════════
#  Main app
# ═══════════════════════════════════════════════════════════════════════════

def main():
    import streamlit as st
    st.set_page_config(page_title="FASTQ QC & V4 Extraction", layout="wide")
    st.title("FASTQ QC & V4 Extraction")

    uploaded_files = st.file_uploader(
        "Upload paired-end FASTQ files (*.fastq.gz)",
        type=["gz"],
        accept_multiple_files=True,
    )

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        n_bases = st.number_input("Prefix length (bp)", value=20, min_value=5, max_value=50)
    with col2:
        read_limit = st.number_input("Read limit for QC stats", value=50000, min_value=1000, step=10000)
    _auto_jobs, _auto_threads = _autodetect_resources()
    with col3:
        threads = st.number_input("Threads (per tool)", value=_auto_threads, min_value=1, max_value=os.cpu_count() or 1)
    with col4:
        jobs = st.number_input("Jobs (parallel samples)", value=_auto_jobs, min_value=1, max_value=os.cpu_count() or 1)

    if not uploaded_files:
        st.info("Upload paired-end FASTQ files (e.g. Sample1_R1.fastq.gz, Sample1_R2.fastq.gz)")
        with st.expander("Supported inputs & 12-case logic matrix"):
            st.markdown(
                "Accepts **V4**, **V3-V4**, or **V4-V5** amplicons. "
                "V-region, read length (250/300), and primer status "
                "(included/trimmed) are auto-detected per sample.\n"
                + CASE_TABLE
            )
        return

    # Group files into samples
    samples = {}
    unmatched = []
    for f in uploaded_files:
        match = re.search(r"(.+)_(R?[12])(?:[_.]|\.fastq).*\.gz$", f.name)
        if match:
            s_name = match.group(1)
            r_type = match.group(2) if match.group(2).startswith("R") else "R" + match.group(2)
            if s_name not in samples:
                samples[s_name] = {"R1": None, "R2": None}
            samples[s_name][r_type] = f
        else:
            unmatched.append(f.name)

    if unmatched:
        st.warning(f"Could not parse sample name from: {', '.join(unmatched)}")
    if not samples:
        st.error("No valid paired-end files found. "
                 "Files must match: *_R1.fastq.gz / *_R2.fastq.gz")
        return

    st.write(f"**{len(samples)} sample(s) detected:** "
             f"{', '.join(sorted(samples.keys()))}")

    if not st.button("Run QC & Extract V4", type="primary"):
        return

    # Write uploaded files to temp directory
    tmpdir = tempfile.mkdtemp()
    _spawn = multiprocessing.get_context("spawn")
    try:
        sample_paths = {}
        for s_name, files in samples.items():
            sample_paths[s_name] = {"R1": None, "R2": None}
            for r_type in ["R1", "R2"]:
                if files[r_type] is not None:
                    path = os.path.join(tmpdir, f"{s_name}_{r_type}.fastq.gz")
                    with open(path, "wb") as out:
                        out.write(files[r_type].getvalue())
                    sample_paths[s_name][r_type] = path

        # ── Phase 1: QC + auto-detect ───────────────────────────────────
        sorted_names = sorted(sample_paths.keys())
        progress = st.progress(0, text="Processing samples...")
        rows_map = {}
        with ProcessPoolExecutor(max_workers=jobs, mp_context=_spawn) as executor:
            future_to_name = {
                executor.submit(
                    _process_qc_sample, s_name, sample_paths[s_name],
                    n_bases, read_limit, threads
                ): s_name
                for s_name in sorted_names
            }
            done = 0
            for future in as_completed(future_to_name):
                s_name = future_to_name[future]
                rows_map[s_name] = future.result()
                done += 1
                progress.progress(
                    done / len(sorted_names),
                    text=f"QC: {done}/{len(sorted_names)} done...",
                )
        rows = [rows_map[n] for n in sorted_names]
        progress.empty()

        # ── Show QC table ───────────────────────────────────────────────
        st.subheader("QC Results")
        html = build_html_table(rows)
        st.markdown(html, unsafe_allow_html=True)

        # ── Validate ────────────────────────────────────────────────────
        invalid_vr = [r for r in rows if r["vregion"] not in VALID_V4_REGIONS]
        error_cases = [r for r in rows if r["cutadapt_error"]]

        if invalid_vr:
            names = ", ".join(
                f'{r["sample"]} ({r["vregion"]})' for r in invalid_vr
            )
            st.error(f"Cannot extract V4 — unsupported V-region: {names}")
            return

        if error_cases:
            for r in error_cases:
                st.error(
                    f'**{r["sample"]}** (Case {r["det_case"]}): '
                    f'{r["cutadapt_error"]}'
                )
            return

        # ── Phase 2: V4 extraction ──────────────────────────────────────
        st.subheader("V4 Extraction")

        trimmed_files = []
        trim_errors = []
        trim_progress = st.progress(0, text="Trimming...")
        with ProcessPoolExecutor(max_workers=jobs, mp_context=_spawn) as executor:
            future_to_row = {
                executor.submit(
                    _process_trim_sample, row, sample_paths[row["sample"]], tmpdir, threads
                ): row
                for row in rows
            }
            done = 0
            for future in as_completed(future_to_row):
                s_name, tag, rlen, case_num, r1_out, r2_out, err = future.result()
                done += 1
                trim_progress.progress(
                    done / len(rows),
                    text=f"Trimming: {done}/{len(rows)} done...",
                )
                if err:
                    trim_errors.append((s_name, err))
                else:
                    trimmed_files.append((s_name, tag, rlen, case_num, r1_out, r2_out))
        trim_progress.empty()
        for s_name, err in trim_errors:
            st.error(f"cutadapt failed for {s_name}: {err}")
        if trim_errors:
            return

        if not trimmed_files:
            st.error("No files were trimmed successfully.")
            return

        # ── Summary with before/after avg lengths ────────────────────────
        rows_by_name = {r["sample"]: r for r in rows}

        summary_html = """<style>
        .trim-table { border-collapse: collapse; width: 100%; font-size: 14px; margin-bottom: 1em; }
        .trim-table th { background-color: #2c3e50; color: white; padding: 8px 12px;
                         text-align: center; border: 1px solid #ddd; }
        .trim-table td { padding: 6px 12px; text-align: center; border: 1px solid #ddd; }
        .trim-table tr:nth-child(even) { background-color: #f8f9fa; }
        .trim-table .sname { text-align: left; font-weight: bold; }
        </style>"""
        summary_html += '<table class="trim-table"><thead><tr>'
        summary_html += "<th>Sample</th><th>Case</th><th>Flags</th>"
        summary_html += "<th>R1 Avg Len (before)</th><th>R1 Avg Len (after)</th>"
        summary_html += "<th>R2 Avg Len (before)</th><th>R2 Avg Len (after)</th>"
        summary_html += "</tr></thead><tbody>"

        for s_name, tag, rlen, case_num, r1_out, r2_out in trimmed_files:
            row = rows_by_name[s_name]
            r1_before = row["r1_avg_len"]
            r2_before = row["r2_avg_len"]
            r1_after = fastq_avg_len(r1_out)
            r2_after = fastq_avg_len(r2_out)

            args_for_display = row["cutadapt_args"]
            flag_str = _format_args(args_for_display)

            summary_html += "<tr>"
            summary_html += f'<td class="sname">{s_name}</td>'
            summary_html += f"<td>{case_num}</td>"
            summary_html += f"<td><code>{flag_str}</code></td>"
            summary_html += f"<td>{r1_before:.1f}</td><td>{r1_after:.1f}</td>"
            summary_html += f"<td>{r2_before:.1f}</td><td>{r2_after:.1f}</td>"
            summary_html += "</tr>"

        summary_html += "</tbody></table>"
        st.markdown(summary_html, unsafe_allow_html=True)

        # ── Zip download ────────────────────────────────────────────────
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_STORED) as zf:
            for s_name, tag, rlen, _, r1_out, r2_out in trimmed_files:
                zf.write(r1_out, f"{s_name}_{tag}_{rlen}bp_trimmed_R1.fastq.gz")
                zf.write(r2_out, f"{s_name}_{tag}_{rlen}bp_trimmed_R2.fastq.gz")
        zip_buffer.seek(0)

        st.success(f"V4 extracted from {len(trimmed_files)} sample(s).")
        st.download_button(
            label="Download V4-trimmed files (.zip)",
            data=zip_buffer,
            file_name="v4_trimmed.zip",
            mime="application/zip",
        )

        with st.expander("12-Case Logic Matrix Reference"):
            st.markdown(CASE_TABLE)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
