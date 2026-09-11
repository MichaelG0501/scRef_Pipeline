#!/usr/bin/env python3
"""
Status: active
Script: analysis/OCCAMS/occams_match_metadata.py
Description: Join archived EGA sample/file/run metadata to OCCAMS SHA subject IDs
  and create the deterministic all-assay file manifest used for BAM selection.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs:
  - ref_outs/OCCAMS/source_data/OCCAMS_metadata.csv
  - ref_outs/OCCAMS/source_data/ega_metadata_files/*_{samples,sample_file,study_experiment_run_sample}.csv
Outputs:
  - tables/: ref_outs/OCCAMS/tables/matched_samples_summary.csv
  - ephemeral compatibility copy: /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/matched_samples_summary.csv
  - logs/: ref_outs/OCCAMS/logs/occams_match_metadata_summary.txt
Cache/replot: Always deterministic and rebuilt; no figures and no replot mode.
Run: /opt/pbs/bin/qsub analysis/OCCAMS/occams_match_metadata.sh
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
import shutil
from collections import defaultdict
from pathlib import Path


WORKDIR = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
SOURCE = WORKDIR / "ref_outs" / "OCCAMS" / "source_data"
EGA_DIR = SOURCE / "ega_metadata_files"
MASTER = SOURCE / "OCCAMS_metadata.csv"
OUT = WORKDIR / "ref_outs" / "OCCAMS" / "tables" / "matched_samples_summary.csv"
EPHEMERAL_OUT = Path(
    "/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/matched_samples_summary.csv"
)
SUMMARY = WORKDIR / "ref_outs" / "OCCAMS" / "logs" / "occams_match_metadata_summary.txt"


def infer_phenotype(row):
    phenotype = row.get("phenotype", "unknown").lower().strip()
    combined = f"{row.get('title', '')} {row.get('description', '')}".lower()
    if "non-tumor" in phenotype or phenotype == "normal":
        return "normal"
    if any(value in phenotype for value in ["tumor", "tumour", "oac", "metastasis"]):
        return "tumor"
    if "barrett" in phenotype or phenotype == "bo":
        return "barretts"
    if "lymph" in phenotype:
        return "lymph"
    if any(value in combined for value in ["tumor", "tumour", "oac"]):
        return "tumor"
    if "barrett" in combined:
        return "barretts"
    if "normal" in combined:
        return "normal"
    if "lymph" in combined:
        return "lymph"
    return "unknown"


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    EPHEMERAL_OUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    with MASTER.open(encoding="utf-8", errors="replace", newline="") as handle:
        master_ids = {row["SHA IDs"].strip() for row in csv.DictReader(handle)}

    all_rows = []
    for samples_file in sorted(EGA_DIR.glob("*_samples.csv")):
        dataset = samples_file.name.split("_")[0]
        sample_file_map = EGA_DIR / f"{dataset}_sample_file.csv"
        run_map = EGA_DIR / f"{dataset}_study_experiment_run_sample.csv"
        for required in [sample_file_map, run_map]:
            if not required.is_file():
                raise FileNotFoundError(required)

        sample_info = {}
        with samples_file.open(encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle):
                matches = {value.strip() for value in row.values() if value and value.strip() in master_ids}
                if len(matches) > 1:
                    raise ValueError(f"Multiple OCCAMS subject IDs in one EGA sample row: {dataset}")
                if matches:
                    sample_info[row.get("accession_id", "")] = (
                        next(iter(matches)),
                        infer_phenotype(row),
                    )

        sample_files = defaultdict(list)
        with sample_file_map.open(encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle):
                sample = row.get("sample_accession_id", "")
                file_accession = row.get("file_accession_id", "")
                if sample and file_accession:
                    sample_files[sample].append(file_accession)

        with run_map.open(encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle):
                sample = row.get("sample_accession_id", "")
                if sample not in sample_info:
                    continue
                subject, phenotype = sample_info[sample]
                for file_accession in sample_files.get(sample, []):
                    all_rows.append(
                        {
                            "File_Accession": file_accession,
                            "Dataset_Accession": dataset,
                            "Subject_ID": subject,
                            "Library_Strategy": row.get("library_strategy", "UNKNOWN"),
                            "Phenotype": phenotype,
                        }
                    )

    grouped = defaultdict(list)
    for row in all_rows:
        grouped[row["File_Accession"]].append(row)
    deduplicated = []
    for file_accession, rows in grouped.items():
        biological = {
            (row["Subject_ID"], row["Library_Strategy"], row["Phenotype"])
            for row in rows
        }
        if len(biological) != 1:
            raise ValueError(f"Conflicting mappings for {file_accession}: {sorted(biological)}")
        deduplicated.append(max(rows, key=lambda row: row["Dataset_Accession"]))

    fieldnames = [
        "File_Accession",
        "Dataset_Accession",
        "Subject_ID",
        "Library_Strategy",
        "Phenotype",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(
            sorted(
                deduplicated,
                key=lambda row: (
                    row["Library_Strategy"],
                    row["Phenotype"],
                    row["Subject_ID"],
                    row["File_Accession"],
                ),
            )
        )
    shutil.copy2(OUT, EPHEMERAL_OUT)

    rnaseq_tumour = [
        row
        for row in deduplicated
        if row["Library_Strategy"] == "RNA-Seq" and row["Phenotype"] == "tumor"
    ]
    subjects = {row["Subject_ID"] for row in rnaseq_tumour}
    if len(rnaseq_tumour) != 302 or len(subjects) != 282:
        raise ValueError(
            f"Expected 302 RNA-seq tumour files/282 subjects; observed "
            f"{len(rnaseq_tumour)}/{len(subjects)}"
        )
    text = (
        "OCCAMS EGA metadata match\n"
        f"deduplicated_file_rows={len(deduplicated)}\n"
        f"rnaseq_tumour_file_rows={len(rnaseq_tumour)}\n"
        f"rnaseq_tumour_subjects={len(subjects)}\n"
        "result=PASS\n"
    )
    SUMMARY.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
