#!/usr/bin/env python3
"""
Status: terminal
Script: analysis/OCCAMS/occams_metadata_coverage_audit.py
Description: Audit downloaded OCCAMS EGAF BAMs against EGA-to-subject mappings,
  the OCCAMS clinical metadata, and the 282-subject BAM mapping.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs:
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/EGAF*/*.bam
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/matched_samples_summary.csv
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/subject_bam_mapping.csv
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/OCCAMS_metadata.csv
Outputs:
  - tables/: ref_outs/OCCAMS/tables/occams_download_metadata_coverage.csv
  - logs/: ref_outs/OCCAMS/logs/occams_metadata_coverage_audit_summary.txt
  - reports/: updates/new_updates/summaries/occams_metadata_coverage_audit_summary.txt
Cache/replot: Deterministic lightweight audit; always rebuilt. Replot-only is not applicable.
Run: /opt/pbs/bin/qsub analysis/OCCAMS/occams_metadata_coverage_audit.sh
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path


WORKDIR = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
EPHEMERAL = Path("/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS")
MATCHED = EPHEMERAL / "matched_samples_summary.csv"
MAPPING = EPHEMERAL / "subject_bam_mapping.csv"
METADATA = EPHEMERAL / "OCCAMS_metadata.csv"
TABLE_DIR = WORKDIR / "ref_outs" / "OCCAMS" / "tables"
LOG_DIR = WORKDIR / "ref_outs" / "OCCAMS" / "logs"
UPDATE_DIR = WORKDIR / "updates" / "new_updates" / "summaries"
AUDIT_OUT = TABLE_DIR / "occams_download_metadata_coverage.csv"
SUMMARY_OUT = LOG_DIR / "occams_metadata_coverage_audit_summary.txt"
UPDATE_SUMMARY = UPDATE_DIR / "occams_metadata_coverage_audit_summary.txt"


def main():
    for path in [MATCHED, MAPPING, METADATA]:
        if not path.is_file():
            raise FileNotFoundError(path)
    for directory in [TABLE_DIR, LOG_DIR, UPDATE_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    with METADATA.open(encoding="utf-8", errors="replace", newline="") as handle:
        metadata_rows = list(csv.DictReader(handle))
    metadata_ids = {row["SHA IDs"].strip() for row in metadata_rows}
    metadata_counts = Counter(row["SHA IDs"].strip() for row in metadata_rows)

    with MATCHED.open(newline="") as handle:
        matched_rows = list(csv.DictReader(handle))
    target_rows = [
        row
        for row in matched_rows
        if row["Library_Strategy"].strip() == "RNA-Seq"
        and row["Phenotype"].strip() == "tumor"
    ]
    target_by_egaf = {row["File_Accession"].strip(): row for row in target_rows}
    if len(target_by_egaf) != len(target_rows):
        raise ValueError("RNA-seq tumour summary contains duplicate EGAF accessions")

    downloaded = {}
    for accession_dir in sorted(EPHEMERAL.glob("EGAF*")):
        if not accession_dir.is_dir():
            continue
        bams = sorted(accession_dir.glob("*.bam"))
        if bams:
            downloaded[accession_dir.name] = bams

    downloaded_subjects = {
        target_by_egaf[egaf]["Subject_ID"].strip()
        for egaf in downloaded
        if egaf in target_by_egaf
    }

    audit_rows = []
    all_accessions = sorted(set(target_by_egaf) | set(downloaded))
    for egaf in all_accessions:
        row = target_by_egaf.get(egaf)
        subject = row["Subject_ID"].strip() if row else ""
        bam_paths = downloaded.get(egaf, [])
        audit_rows.append(
            {
                "File_Accession": egaf,
                "Dataset_Accession": row["Dataset_Accession"].strip() if row else "",
                "Subject_ID": subject,
                "Downloaded": "yes" if bam_paths else "no",
                "Downloaded_BAM_Count": len(bam_paths),
                "BAM_Path": str(bam_paths[0]) if len(bam_paths) == 1 else "",
                "Matched_RNASeq_Tumor_Row": "yes" if row else "no",
                "Metadata_Entry_Available": "yes" if subject in metadata_ids else "no",
                "Metadata_Row_Count": metadata_counts.get(subject, 0),
                "Subject_Represented_By_Download": (
                    "yes" if subject and subject in downloaded_subjects else "no"
                ),
            }
        )

    fieldnames = list(audit_rows[0])
    with AUDIT_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(audit_rows)

    with MAPPING.open(newline="") as handle:
        mapping_rows = list(csv.DictReader(handle))
    mapped_subjects = {row["Subject_ID"].strip() for row in mapping_rows}

    downloaded_without_target = sorted(set(downloaded) - set(target_by_egaf))
    downloaded_subjects_without_metadata = sorted(downloaded_subjects - metadata_ids)
    target_subjects = {row["Subject_ID"].strip() for row in target_rows}
    target_subjects_without_download = sorted(target_subjects - downloaded_subjects)
    mapped_subjects_without_metadata = sorted(mapped_subjects - metadata_ids)
    unavailable_egafs = sorted(set(target_by_egaf) - set(downloaded))

    passed = (
        len(downloaded) == 300
        and all(len(paths) == 1 for paths in downloaded.values())
        and not downloaded_without_target
        and not downloaded_subjects_without_metadata
        and len(target_subjects) == 282
        and not target_subjects_without_download
        and len(mapping_rows) == 282
        and len(mapped_subjects) == 282
        and not mapped_subjects_without_metadata
    )

    lines = [
        "OCCAMS downloaded-sample metadata coverage audit",
        f"metadata_rows={len(metadata_rows)}",
        f"unique_metadata_subjects={len(metadata_ids)}",
        f"matched_rnaseq_tumor_egafs={len(target_rows)}",
        f"matched_rnaseq_tumor_subjects={len(target_subjects)}",
        f"downloaded_egaf_bams={len(downloaded)}",
        f"downloaded_subjects={len(downloaded_subjects)}",
        f"downloaded_egafs_without_rnaseq_tumor_mapping={len(downloaded_without_target)}",
        f"downloaded_subjects_without_metadata={len(downloaded_subjects_without_metadata)}",
        f"target_subjects_without_any_download={len(target_subjects_without_download)}",
        f"mapped_subjects={len(mapped_subjects)}",
        f"mapped_subjects_without_metadata={len(mapped_subjects_without_metadata)}",
        f"authorized_but_unavailable_egafs={len(unavailable_egafs)}",
        f"authorized_but_unavailable_egaf_ids={';'.join(unavailable_egafs)}",
        f"result={'PASS' if passed else 'FAIL'}",
    ]
    text = "\n".join(lines) + "\n"
    SUMMARY_OUT.write_text(text)
    UPDATE_SUMMARY.write_text(text)
    print(text, end="")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
