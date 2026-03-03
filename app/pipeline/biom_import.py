"""
MST-Pipeline — BIOM file import.

Imports a pre-processed BIOM file (HDF5 or JSON) as a Dataset,
creating per-sample csv.gz files compatible with SourceTracker and
Pathogen detection pipelines.
"""
import logging
import re
import shutil
from pathlib import Path

import pandas as pd
from biom import load_table

from app.config import DATASET_DIR, to_relative

logger = logging.getLogger(__name__)

_DNA_RE = re.compile(r"^[ACGTacgt]{50,}$")
_PREFIX_RE = re.compile(r"^[kpcofgs]__")
RANK_NAMES = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def import_biom(biom_file_path: str | Path, name: str = "BIOM Import") -> dict:
    """Import a BIOM file into the pipeline as a completed Dataset.

    Parameters
    ----------
    biom_file_path : path to .biom file (HDF5 or JSON format)
    name : dataset display name

    Returns
    -------
    dict with keys: dataset_id, sample_count, asv_count, variable_region,
                    has_taxonomy
    """
    biom_file_path = Path(biom_file_path)
    table = load_table(str(biom_file_path))

    # --- Validate ---
    sample_ids = list(table.ids(axis="sample"))
    obs_ids = list(table.ids(axis="observation"))

    if len(sample_ids) == 0:
        raise ValueError("BIOM file contains no samples.")
    if len(obs_ids) == 0:
        raise ValueError("BIOM file contains no observations (ASVs).")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("BIOM file contains duplicate sample IDs.")

    # Check for negative counts
    data = table.matrix_data
    if data.min() < 0:
        raise ValueError("BIOM file contains negative counts.")

    # --- Extract sequences ---
    sequences = _extract_sequences(table, obs_ids)

    # --- Extract taxonomy (optional) ---
    taxonomy, has_taxonomy = _extract_taxonomy(table, obs_ids, sequences)

    # --- Auto-detect variable region ---
    variable_region = _detect_region_from_sequences(sequences)

    # --- Create DB records ---
    from app.db.database import SessionLocal
    from app.db.models import Dataset, Sample

    db = SessionLocal()
    try:
        dataset = Dataset(
            upload_id=None,
            name=name,
            status="complete",
            sequencing_type="biom-import",
            variable_region=variable_region,
            sample_count=len(sample_ids),
            asv_count=len(obs_ids),
        )
        db.add(dataset)
        db.flush()
        dataset_id = dataset.id

        output_dir = DATASET_DIR / str(dataset_id)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Copy original BIOM for reference
        shutil.copy2(biom_file_path, output_dir / "original.biom")

        # Write taxonomy files if present
        if has_taxonomy:
            dada2_dir = output_dir / "dada2"
            dada2_dir.mkdir(parents=True, exist_ok=True)

            _write_taxonomy_files(dada2_dir, sequences, taxonomy)
            dataset.taxonomy_path = to_relative(dada2_dir / "taxonomy.tsv")

        # Split into per-sample csv.gz (same format as runner.py:401-418)
        samples_dir = output_dir / "samples"
        samples_dir.mkdir(parents=True, exist_ok=True)

        # Build the full ASV table as a DataFrame (sequences as index)
        asv_df = pd.DataFrame(
            table.matrix_data.toarray(),
            index=sequences,
            columns=sample_ids,
        ).astype(int)

        for sname in sample_ids:
            col = asv_df[sname]
            col_nonzero = col[col > 0]
            sample_path = samples_dir / f"{sname}.csv.gz"
            col_nonzero.to_frame(name=sname).to_csv(
                sample_path, compression="gzip"
            )

            sample = Sample(
                dataset_id=dataset_id,
                sample_name=sname,
                asv_count=len(col_nonzero),
                asv_table_path=to_relative(sample_path),
            )
            db.add(sample)

        db.commit()
        logger.info(
            f"BIOM import complete: dataset {dataset_id}, "
            f"{len(sample_ids)} samples, {len(obs_ids)} ASVs, "
            f"region={variable_region}, taxonomy={has_taxonomy}"
        )

    except Exception:
        db.rollback()
        raise
    finally:
        db.close()

    return {
        "dataset_id": dataset_id,
        "sample_count": len(sample_ids),
        "asv_count": len(obs_ids),
        "variable_region": variable_region,
        "has_taxonomy": has_taxonomy,
    }


def _extract_sequences(table, obs_ids: list[str]) -> list[str]:
    """Extract ASV sequences from a BIOM table.

    Checks three sources in order:
    1. Observation metadata "sequence" key
    2. Observation IDs that look like DNA sequences
    3. Raises ValueError if neither works
    """
    # Try observation metadata "sequence" key
    if table.metadata(axis="observation") is not None:
        first_md = table.metadata(obs_ids[0], axis="observation")
        if first_md and "sequence" in first_md:
            seqs = []
            for oid in obs_ids:
                md = table.metadata(oid, axis="observation")
                seq = md.get("sequence", "") if md else ""
                if not seq:
                    break
                seqs.append(str(seq))
            if len(seqs) == len(obs_ids):
                return seqs

    # Try observation IDs as sequences
    if all(_DNA_RE.match(oid) for oid in obs_ids):
        return list(obs_ids)

    raise ValueError(
        "Cannot extract ASV sequences from BIOM file. "
        "Observation IDs must be DNA sequences (>=50bp) or "
        "observation metadata must contain a 'sequence' key."
    )


def _extract_taxonomy(table, obs_ids: list[str], sequences: list[str]):
    """Extract taxonomy from observation metadata if present.

    Returns (taxonomy_dict, has_taxonomy) where taxonomy_dict maps
    sequence -> list of 7 ranks.
    """
    taxonomy = {}
    has_taxonomy = False

    if table.metadata(axis="observation") is None:
        return taxonomy, False

    first_md = table.metadata(obs_ids[0], axis="observation")
    if not first_md or "taxonomy" not in first_md:
        return taxonomy, False

    for oid, seq in zip(obs_ids, sequences):
        md = table.metadata(oid, axis="observation")
        if not md or "taxonomy" not in md:
            continue

        raw = md["taxonomy"]
        ranks = _parse_taxonomy(raw)
        if any(r != "NA" for r in ranks):
            has_taxonomy = True
        taxonomy[seq] = ranks

    return taxonomy, has_taxonomy


def _parse_taxonomy(raw) -> list[str]:
    """Normalize taxonomy to a list of 7 ranks (Kingdom..Species)."""
    if isinstance(raw, list):
        ranks = [str(r).strip() for r in raw]
    elif isinstance(raw, str):
        ranks = [r.strip() for r in raw.split(";")]
    else:
        return ["NA"] * 7

    # Strip prefix like k__, p__, etc.
    cleaned = []
    for r in ranks:
        r = _PREFIX_RE.sub("", r).strip()
        if not r or r.lower() in ("", "unclassified", "unidentified"):
            r = "NA"
        cleaned.append(r)

    # Pad or truncate to 7 ranks
    while len(cleaned) < 7:
        cleaned.append("NA")
    return cleaned[:7]


def _write_taxonomy_files(
    dada2_dir: Path,
    sequences: list[str],
    taxonomy: dict[str, list[str]],
):
    """Write taxonomy.tsv and rep_seqs.fasta in the same format as run_taxonomy.R."""
    # taxonomy.tsv: ASV_ID, Kingdom..Species
    rows = []
    with open(dada2_dir / "rep_seqs.fasta", "w") as fasta:
        for i, seq in enumerate(sequences, 1):
            asv_id = f"ASV_{i}"
            fasta.write(f">{asv_id}\n{seq}\n")
            ranks = taxonomy.get(seq, ["NA"] * 7)
            rows.append([asv_id] + ranks)

    tax_df = pd.DataFrame(rows, columns=["ASV_ID"] + RANK_NAMES)
    tax_df.to_csv(dada2_dir / "taxonomy.tsv", sep="\t", index=False)


def _detect_region_from_sequences(sequences: list[str]) -> str | None:
    """Auto-detect variable region from ASV sequences using primer matching + length."""
    from app.pipeline.detect import REGION_PRIMERS, _primer_matches

    if not sequences:
        return None

    avg_len = sum(len(s) for s in sequences) / len(sequences)

    # Length-based heuristic
    if avg_len >= 1000:
        return "V1-V9"

    # Try primer matching on a sample of sequences
    sample_seqs = sequences[:200]
    region_scores = {}
    for region_name, info in REGION_PRIMERS.items():
        lo, hi = info["expected_len"]
        count = 0
        for seq in sample_seqs:
            if _primer_matches(seq, info["forward"]):
                count += 1
        region_scores[region_name] = count

    best = max(region_scores, key=region_scores.get)
    if region_scores[best] > 0:
        return best

    # Fallback: length heuristic for common regions
    if 200 <= avg_len <= 300:
        return "V4"
    if 350 <= avg_len <= 500:
        return "V3-V4"
    if 300 <= avg_len <= 450:
        return "V4-V5"

    return None
