#!/usr/bin/env python3
"""
Status: terminal
Script: analysis/OCCAMS/occams_build_count_matrix.py
Description: Convert the corrected GRCh37 featureCounts table into the canonical
  gzipped gene-by-subject OCCAMS matrix and persistent sample/QC tables.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs:
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_raw_counts.txt
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_raw_counts.txt.summary
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/subject_bam_mapping.csv
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/matched_samples_summary.csv
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/OCCAMS_metadata.csv
Outputs:
  - intermediate/: not used; the reconstructable featureCounts table remains ephemeral
  - tables/: ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz;
    ref_outs/OCCAMS/tables/{subject_bam_mapping,matched_samples_summary,
    OCCAMS_RNAseq_GRCh37_featurecounts_qc,OCCAMS_count_matrix_metadata_rows}.csv
  - source_data/: ref_outs/OCCAMS/source_data/OCCAMS_metadata.csv
  - logs/: ref_outs/OCCAMS/logs/occams_count_matrix_build_summary.txt
  - reports/: updates/new_updates/summaries/occams_count_matrix_build_summary.txt
Cache/replot: Existing complete outputs are reused unless SCREF_FORCE_REBUILD=TRUE.
  Replot-only is not applicable because this script makes no figures.
Run: Called by analysis/OCCAMS/occams_featurecounts_grch37.sh after featureCounts.
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
import gzip
import os
import shutil
import re
from collections import defaultdict
from pathlib import Path


WORKDIR = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
EPHEMERAL = Path("/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS")
RAW_COUNTS = EPHEMERAL / "counts" / "OCCAMS_RNAseq_GRCh37_raw_counts.txt"
RAW_SUMMARY = Path(f"{RAW_COUNTS}.summary")
MAPPING = EPHEMERAL / "subject_bam_mapping.csv"
MATCHED_SUMMARY = EPHEMERAL / "matched_samples_summary.csv"
MASTER_METADATA = EPHEMERAL / "OCCAMS_metadata.csv"

OUT_ROOT = WORKDIR / "ref_outs" / "OCCAMS"
COUNTS_OUT = OUT_ROOT / "counts" / "OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz"
GTF_PATH = OUT_ROOT / "source_data" / "GENCODE_v19_GRCh37" / "gencode.v19.annotation.gtf.gz"
CIBERSORTX_OUT = OUT_ROOT / "counts" / "OCCAMS_RNAseq_GRCh37_CIBERSORTx_Mixture.txt"
TABLE_DIR = OUT_ROOT / "tables"
SOURCE_DIR = OUT_ROOT / "source_data"
LOG_DIR = OUT_ROOT / "logs"
UPDATE_DIR = WORKDIR / "updates" / "new_updates" / "summaries"
QC_OUT = TABLE_DIR / "OCCAMS_RNAseq_GRCh37_featurecounts_qc.csv"
METADATA_ROWS_OUT = TABLE_DIR / "OCCAMS_count_matrix_metadata_rows.csv"
SUMMARY_OUT = LOG_DIR / "occams_count_matrix_build_summary.txt"
UPDATE_SUMMARY = UPDATE_DIR / "occams_count_matrix_build_summary.txt"


def load_mapping():
    by_bam = {}
    rows = []
    with MAPPING.open(newline="") as handle:
        for row in csv.DictReader(handle):
            bam_name = Path(row["BAM_Path"]).name
            if bam_name in by_bam:
                raise ValueError(f"Duplicate BAM basename in mapping: {bam_name}")
            by_bam[bam_name] = row["Subject_ID"].strip()
            rows.append(row)
    if len(rows) != 282 or len(set(by_bam.values())) != 282:
        raise ValueError(
            f"Expected 282 unique subject mappings; observed {len(rows)} rows and "
            f"{len(set(by_bam.values()))} subjects"
        )
    return by_bam, rows


def load_metadata_rows():
    with MASTER_METADATA.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if not fieldnames or fieldnames[0] != "SHA IDs":
        raise ValueError("OCCAMS metadata first column is not 'SHA IDs'")
    return rows, fieldnames


def read_featurecounts_header():
    with RAW_COUNTS.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                return line.rstrip("\n").split("\t")
    raise ValueError("No featureCounts header found")


def load_gene_map():
    gene_map = {}
    with gzip.open(GTF_PATH, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) > 8 and fields[2] == "gene":
                attrs = fields[8]
                m_id = re.search(r'gene_id "([^"]+)"', attrs)
                m_name = re.search(r'gene_name "([^"]+)"', attrs)
                if m_id and m_name:
                    gene_map[m_id.group(1)] = m_name.group(1)
    return gene_map


def write_matrix(by_bam, gene_map):
    header = read_featurecounts_header()
    if header[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
        raise ValueError(f"Unexpected featureCounts annotation columns: {header[:6]}")
    bam_names = [Path(value).name for value in header[6:]]
    missing = sorted(set(bam_names) - set(by_bam))
    if missing:
        raise ValueError(f"Count columns missing from subject mapping: {missing[:5]}")
    subjects = [by_bam[name] for name in bam_names]
    if len(subjects) != 282 or len(set(subjects)) != 282:
        raise ValueError(
            f"Expected 282 unique matrix subjects; observed {len(subjects)} columns "
            f"and {len(set(subjects))} unique subjects"
        )

    n_genes = 0
    cibersortx_counts = defaultdict(lambda: [0] * len(subjects))
    
    with RAW_COUNTS.open() as source, gzip.open(COUNTS_OUT, "wt", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["Geneid"] + subjects)
        data_started = False
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not data_started:
                data_started = True
                continue
            if len(fields) != 6 + len(subjects):
                raise ValueError(
                    f"Malformed featureCounts row {n_genes + 1}: "
                    f"expected {6 + len(subjects)} fields, observed {len(fields)}"
                )
            writer.writerow([fields[0]] + fields[6:])
            
            gene_id = fields[0]
            if gene_id in gene_map:
                symbol = gene_map[gene_id]
                if symbol:
                    counts = cibersortx_counts[symbol]
                    for i in range(len(subjects)):
                        counts[i] += int(fields[6 + i])
                        
            n_genes += 1

    with CIBERSORTX_OUT.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["GeneSymbol"] + subjects)
        for symbol, counts in cibersortx_counts.items():
            writer.writerow([symbol] + counts)

    return subjects, n_genes


def write_qc(by_bam):
    with RAW_SUMMARY.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "Status" not in rows[0]:
        raise ValueError("Malformed featureCounts summary")

    sample_columns = [column for column in rows[0] if column != "Status"]
    output_rows = []
    for sample_path in sample_columns:
        bam_name = Path(sample_path).name
        subject = by_bam.get(bam_name)
        if subject is None:
            raise ValueError(f"QC sample missing from subject mapping: {bam_name}")
        total = sum(int(row[sample_path]) for row in rows)
        for row in rows:
            count = int(row[sample_path])
            output_rows.append(
                {
                    "Subject_ID": subject,
                    "BAM_Name": bam_name,
                    "Status": row["Status"],
                    "Count": count,
                    "Percent": 100 * count / total if total else 0,
                }
            )

    with QC_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Subject_ID", "BAM_Name", "Status", "Count", "Percent"],
        )
        writer.writeheader()
        writer.writerows(output_rows)


def write_metadata_rows(subjects, metadata_rows, metadata_fields):
    subject_set = set(subjects)
    metadata_ids = {row["SHA IDs"].strip() for row in metadata_rows}
    missing = sorted(subject_set - metadata_ids)
    if missing:
        raise ValueError(f"Matrix subjects absent from OCCAMS metadata: {missing[:5]}")
    selected = [row for row in metadata_rows if row["SHA IDs"].strip() in subject_set]
    with METADATA_ROWS_OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fields)
        writer.writeheader()
        writer.writerows(selected)
    return len(selected), len({row["SHA IDs"].strip() for row in selected})


def main():
    force = os.environ.get("SCREF_FORCE_REBUILD", "FALSE").upper() == "TRUE"
    required = [RAW_COUNTS, RAW_SUMMARY, MAPPING, MATCHED_SUMMARY, MASTER_METADATA]
    missing_inputs = [str(path) for path in required if not path.is_file()]
    if missing_inputs:
        raise FileNotFoundError(f"Missing required inputs: {missing_inputs}")

    for directory in [COUNTS_OUT.parent, TABLE_DIR, SOURCE_DIR, LOG_DIR, UPDATE_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    reusable = [COUNTS_OUT, QC_OUT, METADATA_ROWS_OUT, SUMMARY_OUT, CIBERSORTX_OUT]
    if not force and all(path.is_file() and path.stat().st_size > 0 for path in reusable):
        print("Complete OCCAMS count-matrix outputs already exist; set SCREF_FORCE_REBUILD=TRUE to rebuild.")
        return

    by_bam, mapping_rows = load_mapping()
    metadata_rows, metadata_fields = load_metadata_rows()
    gene_map = load_gene_map()
    subjects, n_genes = write_matrix(by_bam, gene_map)
    write_qc(by_bam)
    n_metadata_rows, n_metadata_subjects = write_metadata_rows(
        subjects, metadata_rows, metadata_fields
    )

    shutil.copy2(MASTER_METADATA, SOURCE_DIR / "OCCAMS_metadata.csv")
    shutil.copy2(MAPPING, TABLE_DIR / "subject_bam_mapping.csv")
    shutil.copy2(MATCHED_SUMMARY, TABLE_DIR / "matched_samples_summary.csv")

    summary_lines = [
        "OCCAMS corrected GRCh37 count-matrix build",
        f"genes={n_genes}",
        f"matrix_subjects={len(subjects)}",
        f"unique_matrix_subjects={len(set(subjects))}",
        f"metadata_rows_for_matrix_subjects={n_metadata_rows}",
        f"unique_metadata_subjects_for_matrix={n_metadata_subjects}",
        f"count_matrix={COUNTS_OUT}",
        "result=PASS",
    ]
    summary_text = "\n".join(summary_lines) + "\n"
    SUMMARY_OUT.write_text(summary_text)
    UPDATE_SUMMARY.write_text(summary_text)
    print(summary_text, end="")


if __name__ == "__main__":
    main()
