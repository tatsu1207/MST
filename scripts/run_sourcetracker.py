#!/usr/bin/env python3
"""Command-line SourceTracker2 runner.

Usage:
    python3 run_sourcetracker.py <sink_csv_gz> <source_csv_gz> <source_fasta> <source_design> <output_dir>
                                 [--mode asv|otu] [--src-depth N] [--snk-depth N]
                                 [--restarts N] [--burnin N] [--draws N] [--threads N]
"""

import sys, os, argparse, re, tempfile, subprocess
from collections import defaultdict

import numpy as np
import pandas as pd

# ── Add script directory to path so we can import app functions ───────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from app_sourcetracker_ui import (
    _load_csv_gz, _load_fasta, _load_design,
    align_features, collapse_by_group, make_figure, _fig_to_png
)


def parse_args():
    p = argparse.ArgumentParser(description="Run SourceTracker2 Gibbs from CLI")
    p.add_argument("sink",    help="Sink table CSV.gz (rows=features, cols=samples)")
    p.add_argument("source",  help="Source table CSV.gz (rows=ASV IDs, cols=samples)")
    p.add_argument("fasta",   help="Source FASTA (ASV ID → sequence)")
    p.add_argument("design",  help="Source design TSV (Sample<tab>Group)")
    p.add_argument("outdir",  help="Output directory")
    p.add_argument("--mode",       default="asv", choices=["asv","otu"],
                   help="Feature matching: asv (100%%) or otu (99%%) [default: asv]")
    p.add_argument("--groups", nargs="+", default=None, metavar="GROUP",
                   help="Source groups to include (space-separated). "
                        "If omitted, all groups in the design file are used. "
                        "Example: --groups human pig cow")
    p.add_argument("--src-depth",  type=int, default=1000, help="Source rarefaction depth [1000]")
    p.add_argument("--snk-depth",  type=int, default=0,    help="Sink rarefaction depth (0=auto) [0]")
    p.add_argument("--restarts",   type=int, default=10,   help="MCMC restarts [10]")
    p.add_argument("--burnin",     type=int, default=100,  help="Burnin iterations [100]")
    p.add_argument("--draws",      type=int, default=1,    help="Draws per restart [1]")
    p.add_argument("--threads",    type=int, default=min(16, os.cpu_count() or 4),
                   help="vsearch threads")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("Loading files...", flush=True)
    sink_df   = _load_csv_gz(args.sink)
    source_df = _load_csv_gz(args.source)
    db_fasta  = _load_fasta(args.fasta)
    design    = _load_design(args.design)

    print(f"  Sink:   {sink_df.shape[0]} features × {sink_df.shape[1]} samples", flush=True)
    print(f"  Source: {source_df.shape[0]} ASVs × {source_df.shape[1]} samples", flush=True)
    all_groups = sorted(set(design.values()))
    print(f"  Available source groups ({len(all_groups)}): {', '.join(all_groups)}", flush=True)

    # Filter groups if --groups was specified
    if args.groups:
        unknown = [g for g in args.groups if g not in all_groups]
        if unknown:
            print(f"ERROR: Unknown group(s): {', '.join(unknown)}", file=sys.stderr)
            print(f"       Available: {', '.join(all_groups)}", file=sys.stderr)
            sys.exit(1)
        design = {s: g for s, g in design.items() if g in args.groups}
        groups = sorted(set(design.values()))
        print(f"  Selected groups ({len(groups)}): {', '.join(groups)}", flush=True)
    else:
        groups = all_groups

    print(f"\nAligning features (mode={args.mode})...", flush=True)
    sink_al, src_al, (n_m, n_t) = align_features(
        sink_df, source_df, db_fasta, mode=args.mode, threads=args.threads
    )
    pct = 100 * n_m / max(n_t, 1)
    print(f"  Matched {n_m}/{n_t} sink features ({pct:.1f}%)", flush=True)

    print("\nCollapsing source samples by group...", flush=True)
    # Count samples per group before collapsing
    sample_counts = {}
    for s, g in design.items():
        if s in src_al.index:
            sample_counts[g] = sample_counts.get(g, 0) + 1
    src_col = collapse_by_group(src_al, design)
    for g in src_col.index:
        n = sample_counts.get(g, 0)
        print(f"  {g}: {n} samples, {int(src_col.loc[g].sum()):,} reads", flush=True)

    # Auto sink depth: 50% of min sample total after alignment
    if args.snk_depth == 0:
        snk_totals = sink_al.sum(axis=1)
        args.snk_depth = max(500, int(snk_totals.min() * 0.5))
        print(f"\nAuto sink depth: {args.snk_depth}", flush=True)

    print(f"\nRunning Gibbs sampling...", flush=True)
    print(f"  restarts={args.restarts} burnin={args.burnin} draws={args.draws}", flush=True)
    print(f"  src_depth={args.src_depth} snk_depth={args.snk_depth}", flush=True)

    _GIBBS = os.path.join(_HERE, "_run_gibbs.py")
    with tempfile.TemporaryDirectory() as tmp:
        src_csv = os.path.join(tmp, "sources.csv")
        snk_csv = os.path.join(tmp, "sinks.csv")
        res_dir = os.path.join(tmp, "results")

        src_col.to_csv(src_csv)
        sink_al.to_csv(snk_csv)

        cmd = [
            "conda", "run", "--no-capture-output", "-n", "ST",
            "python3", _GIBBS,
            src_csv, snk_csv, res_dir,
            str(args.src_depth), str(args.snk_depth),
            str(args.restarts), str(args.burnin), str(args.draws),
        ]
        result = subprocess.run(cmd, capture_output=False, text=True)

        if result.returncode != 0:
            print("ERROR: Gibbs sampler failed.", file=sys.stderr)
            sys.exit(1)

        mp  = pd.read_csv(os.path.join(res_dir, "proportions.csv"), index_col=0)
        mps = pd.read_csv(os.path.join(res_dir, "stds.csv"),        index_col=0)

    # Save outputs
    mp.to_csv(os.path.join(args.outdir, "proportions.csv"))
    mps.to_csv(os.path.join(args.outdir, "stds.csv"))
    pct = (mp * 100).round(2)
    pct.to_csv(os.path.join(args.outdir, "contributions_pct.csv"))
    print(f"\nSaved proportions.csv, stds.csv, contributions_pct.csv → {args.outdir}", flush=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        fig = make_figure(mp, "Source Contribution (SourceTracker2)")
        png = _fig_to_png(fig)
        with open(os.path.join(args.outdir, "sourcetracker_results.png"), "wb") as f:
            f.write(png)
        print(f"Saved sourcetracker_results.png → {args.outdir}", flush=True)
    except Exception as e:
        print(f"Warning: could not save figure: {e}", flush=True)

    print("\nDone.", flush=True)
    print("\n--- Source Contributions (%) ---")
    print(pct.to_string())


if __name__ == "__main__":
    main()
