#!/usr/bin/env python3
"""
Status: active
Script: analysis/OCCAMS/occams_prepare_subject_bams.py
Description: Resolve downloaded RNA-seq tumour EGAF BAMs to 282 OCCAMS subjects,
  generate deterministic symlink/merge commands, and persist the subject mapping.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs:
  - ref_outs/OCCAMS/tables/matched_samples_summary.csv
  - /rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/EGAF*/*.bam
Outputs:
  - intermediate/: ephemeral `merge_commands.sh` and `merged_bams/*.bam`
  - tables/: ref_outs/OCCAMS/tables/subject_bam_mapping.csv and
    ref_outs/OCCAMS/tables/occams_unavailable_egafs.csv
  - ephemeral compatibility copy: subject_bam_mapping.csv
  - logs/: ref_outs/OCCAMS/logs/occams_prepare_subject_bams_summary.txt
Cache/replot: Commands/mapping are deterministically rebuilt; BAM creation is handled
  by the merge wrapper. Replot-only is not applicable.
Run: /opt/pbs/bin/qsub analysis/OCCAMS/occams_merge_subject_bams.sh
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
import os
import shlex
import shutil
from collections import defaultdict
from pathlib import Path


WORKDIR = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
EPHEMERAL = Path("/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS")
MATCHED = WORKDIR / "ref_outs" / "OCCAMS" / "tables" / "matched_samples_summary.csv"
MERGED = EPHEMERAL / "merged_bams"
COMMANDS = EPHEMERAL / "merge_commands.sh"
LIVE_MAPPING = WORKDIR / "ref_outs" / "OCCAMS" / "tables" / "subject_bam_mapping.csv"
EPHEMERAL_MAPPING = EPHEMERAL / "subject_bam_mapping.csv"
UNAVAILABLE = WORKDIR / "ref_outs" / "OCCAMS" / "tables" / "occams_unavailable_egafs.csv"
SUMMARY = WORKDIR / "ref_outs" / "OCCAMS" / "logs" / "occams_prepare_subject_bams_summary.txt"


def main():
    MERGED.mkdir(parents=True, exist_ok=True)
    LIVE_MAPPING.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with MATCHED.open(newline="") as handle:
        targets = [
            row
            for row in csv.DictReader(handle)
            if row["Library_Strategy"].strip() == "RNA-Seq"
            and row["Phenotype"].strip() == "tumor"
        ]

    subject_bams = defaultdict(list)
    unavailable_rows = []
    for row in targets:
        accession = row["File_Accession"].strip()
        subject = row["Subject_ID"].strip()
        bams = sorted((EPHEMERAL / accession).glob("*.bam"))
        if len(bams) == 1:
            subject_bams[subject].append(bams[0])
        elif len(bams) == 0:
            unavailable_rows.append(
                {"File_Accession": accession, "Subject_ID": subject, "Reason": "not_downloaded"}
            )
        else:
            raise ValueError(f"Expected at most one BAM in {accession}; observed {len(bams)}")

    target_subjects = {row["Subject_ID"].strip() for row in targets}
    missing_subjects = sorted(target_subjects - set(subject_bams))
    if missing_subjects:
        raise ValueError(f"Subjects with no downloaded BAM: {missing_subjects[:5]}")
    if len(subject_bams) != 282:
        raise ValueError(f"Expected 282 subjects; observed {len(subject_bams)}")

    short_ids = {subject: subject[:16] for subject in subject_bams}
    if len(set(short_ids.values())) != len(short_ids):
        raise ValueError("A 16-character subject-ID collision was detected")

    mapping_rows = []
    command_lines = ["#!/bin/bash", "set -euo pipefail", ""]
    for subject in sorted(subject_bams):
        bams = sorted(set(subject_bams[subject]))
        destination = MERGED / f"{short_ids[subject]}.bam"
        if len(bams) == 1:
            command_lines.append(
                f"ln -sfn {shlex.quote(str(bams[0]))} {shlex.quote(str(destination))}"
            )
        else:
            inputs = " ".join(shlex.quote(str(path)) for path in bams)
            command_lines.append(
                f"samtools merge -@ 4 -f {shlex.quote(str(destination))} {inputs}"
            )
        mapping_rows.append(
            {
                "Subject_ID": subject,
                "Short_ID": short_ids[subject],
                "BAM_Path": str(destination),
                "N_Input_Files": len(bams),
            }
        )

    COMMANDS.write_text("\n".join(command_lines) + "\n")
    os.chmod(COMMANDS, 0o750)
    fields = ["Subject_ID", "Short_ID", "BAM_Path", "N_Input_Files"]
    with LIVE_MAPPING.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(mapping_rows)
    shutil.copy2(LIVE_MAPPING, EPHEMERAL_MAPPING)
    with UNAVAILABLE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["File_Accession", "Subject_ID", "Reason"])
        writer.writeheader()
        writer.writerows(unavailable_rows)

    multi = sum(row["N_Input_Files"] > 1 for row in mapping_rows)
    text = (
        "OCCAMS subject BAM preparation\n"
        f"subjects=282\n"
        f"single_file_subjects={282 - multi}\n"
        f"multi_file_subjects={multi}\n"
        f"unavailable_egafs={len(unavailable_rows)}\n"
        "short_id_collisions=0\n"
        "result=PASS\n"
    )
    SUMMARY.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
