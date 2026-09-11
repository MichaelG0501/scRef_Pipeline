#!/usr/bin/env python3
"""
Status: terminal
Script: analysis/OCCAMS/occams_validate_final_outputs.py
Description: Exhaustively validate the canonical OCCAMS matrix, subject mapping,
  featureCounts QC, metadata coverage, GTF checksum, and persistent deliverables.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs:
  - /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz
  - /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/OCCAMS/tables/{subject_bam_mapping,OCCAMS_RNAseq_GRCh37_featurecounts_qc,OCCAMS_count_matrix_metadata_rows}.csv
  - /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/OCCAMS/source_data/{OCCAMS_metadata.csv,GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz}
Outputs:
  - tables/: ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
  - logs/: ref_outs/OCCAMS/logs/occams_final_output_validation_summary.txt
  - reports/: updates/new_updates/summaries/occams_final_output_validation_summary.txt
Cache/replot: Deterministic exhaustive validation; always rebuilt. No figures.
Run: /opt/pbs/bin/qsub analysis/OCCAMS/occams_validate_final_outputs.sh
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
import gzip
import hashlib
import statistics
from collections import defaultdict
from pathlib import Path


WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
ROOT = WD / "ref_outs" / "OCCAMS"
MATRIX = ROOT / "counts" / "OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz"
MAPPING = ROOT / "tables" / "subject_bam_mapping.csv"
QC = ROOT / "tables" / "OCCAMS_RNAseq_GRCh37_featurecounts_qc.csv"
QC_SUMMARY = ROOT / "tables" / "OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv"
MATRIX_METADATA = ROOT / "tables" / "OCCAMS_count_matrix_metadata_rows.csv"
MASTER_METADATA = ROOT / "source_data" / "OCCAMS_metadata.csv"
GTF = ROOT / "source_data" / "GENCODE_v19_GRCh37" / "gencode.v19.annotation.gtf.gz"
SUMMARY = ROOT / "logs" / "occams_final_output_validation_summary.txt"
UPDATE_SUMMARY = WD / "updates" / "new_updates" / "summaries" / "occams_final_output_validation_summary.txt"
EXPECTED_GTF_MD5 = "bd83e28270e595d3bde6bfcb21c9748f"


def md5(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def main():
    required = [MATRIX, MAPPING, QC, MATRIX_METADATA, MASTER_METADATA, GTF]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"Missing or empty required outputs: {missing}")

    with MAPPING.open(newline="") as handle:
        mapping_rows = list(csv.DictReader(handle))
    mapping_subjects = [row["Subject_ID"].strip() for row in mapping_rows]
    if len(mapping_subjects) != 282 or len(set(mapping_subjects)) != 282:
        raise ValueError("Subject mapping is not 282 unique subjects")

    gene_ids = set()
    sample_totals = None
    n_genes = 0
    with gzip.open(str(MATRIX), "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        matrix_subjects = header[1:]
        if len(matrix_subjects) != 282 or len(set(matrix_subjects)) != 282:
            raise ValueError("Matrix header is not 282 unique subjects")
        if set(matrix_subjects) != set(mapping_subjects):
            raise ValueError("Matrix subjects do not equal subject mapping")
        sample_totals = [0] * len(matrix_subjects)
        for row in reader:
            if len(row) != 283:
                raise ValueError(f"Matrix row {n_genes + 2} has {len(row)} columns")
            gene = row[0]
            if gene in gene_ids:
                raise ValueError(f"Duplicate gene identifier: {gene}")
            gene_ids.add(gene)
            try:
                values = [int(value) for value in row[1:]]
            except ValueError:
                raise ValueError(f"Non-integer count at matrix row {n_genes + 2}")
            if any(value < 0 for value in values):
                raise ValueError(f"Negative count at matrix row {n_genes + 2}")
            sample_totals = [total + value for total, value in zip(sample_totals, values)]
            n_genes += 1
    if n_genes != 57820 or any(total <= 0 for total in sample_totals):
        raise ValueError(f"Unexpected genes or empty sample: genes={n_genes}")

    with MASTER_METADATA.open(encoding="utf-8", errors="replace", newline="") as handle:
        master_rows = list(csv.DictReader(handle))
    master_subjects = {row["SHA IDs"].strip() for row in master_rows}
    if not set(matrix_subjects).issubset(master_subjects):
        raise ValueError("A matrix subject is absent from master metadata")
    with MATRIX_METADATA.open(encoding="utf-8", errors="replace", newline="") as handle:
        selected_rows = list(csv.DictReader(handle))
    selected_subjects = {row["SHA IDs"].strip() for row in selected_rows}
    if selected_subjects != set(matrix_subjects):
        raise ValueError("Selected metadata rows do not cover exactly the matrix subjects")

    qc_by_subject = defaultdict(dict)
    pct_sums = defaultdict(float)
    with QC.open(newline="") as handle:
        qc_rows = list(csv.DictReader(handle))
    for row in qc_rows:
        subject = row["Subject_ID"].strip()
        status = row["Status"].strip()
        if status in qc_by_subject[subject]:
            raise ValueError(f"Duplicate QC status for subject {subject}: {status}")
        qc_by_subject[subject][status] = float(row["Percent"])
        pct_sums[subject] += float(row["Percent"])
    if set(qc_by_subject) != set(matrix_subjects):
        raise ValueError("QC subjects do not equal matrix subjects")
    if any(abs(total - 100) > 0.01 for total in pct_sums.values()):
        raise ValueError("FeatureCounts QC percentages do not sum to 100")
    required_statuses = {"Assigned", "Unassigned_NoFeatures", "Unassigned_MultiMapping"}
    if any(not required_statuses.issubset(statuses) for statuses in qc_by_subject.values()):
        raise ValueError("Required featureCounts QC statuses are missing")

    assigned = [statuses["Assigned"] for statuses in qc_by_subject.values()]
    no_features = [statuses["Unassigned_NoFeatures"] for statuses in qc_by_subject.values()]
    multimapping = [statuses["Unassigned_MultiMapping"] for statuses in qc_by_subject.values()]
    compatible = [a + m for a, m in zip(assigned, multimapping)]
    qc_summary_rows = []
    for subject in matrix_subjects:
        statuses = qc_by_subject[subject]
        annotation_compatible = (
            statuses["Assigned"] + statuses["Unassigned_MultiMapping"]
        )
        qc_summary_rows.append(
            {
                "Subject_ID": subject,
                "Assigned_Percent": statuses["Assigned"],
                "No_Features_Percent": statuses["Unassigned_NoFeatures"],
                "Multi_Mapping_Percent": statuses["Unassigned_MultiMapping"],
                "Annotation_Compatible_Percent": annotation_compatible,
                "Annotation_Compatible_Below_50_Flag": (
                    "yes" if annotation_compatible < 50 else "no"
                ),
            }
        )
    with QC_SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_summary_rows[0]))
        writer.writeheader()
        writer.writerows(qc_summary_rows)
    gtf_md5 = md5(GTF)
    if gtf_md5 != EXPECTED_GTF_MD5:
        raise ValueError(f"GTF MD5 mismatch: {gtf_md5}")

    lines = [
        "OCCAMS final output validation",
        f"matrix_bytes={MATRIX.stat().st_size}",
        f"genes={n_genes}",
        f"matrix_subjects={len(matrix_subjects)}",
        f"unique_matrix_subjects={len(set(matrix_subjects))}",
        f"subjects_with_positive_total_counts={sum(total > 0 for total in sample_totals)}",
        f"metadata_rows_for_matrix_subjects={len(selected_rows)}",
        f"unique_metadata_subjects_for_matrix={len(selected_subjects)}",
        f"qc_rows={len(qc_rows)}",
        f"assigned_pct_min={min(assigned):.6f}",
        f"assigned_pct_median={statistics.median(assigned):.6f}",
        f"assigned_pct_max={max(assigned):.6f}",
        f"no_features_pct_median={statistics.median(no_features):.6f}",
        f"annotation_compatible_pct_min={min(compatible):.6f}",
        f"annotation_compatible_pct_median={statistics.median(compatible):.6f}",
        f"subjects_annotation_compatible_below_50={sum(value < 50 for value in compatible)}",
        f"lowest_annotation_compatible_subject={min(qc_summary_rows, key=lambda row: row['Annotation_Compatible_Percent'])['Subject_ID']}",
        f"gtf_md5={gtf_md5}",
        "result=PASS",
    ]
    text = "\n".join(lines) + "\n"
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    UPDATE_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(text)
    UPDATE_SUMMARY.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
