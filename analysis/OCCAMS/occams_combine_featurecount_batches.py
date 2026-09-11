#!/usr/bin/env python3
"""
Status: active
Script: analysis/OCCAMS/occams_combine_featurecount_batches.py
Description: Combine eight validated featureCounts batches into the canonical
  282-column ephemeral table and summary in subject-mapping order.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs: ref_outs/OCCAMS/tables/{subject_bam_mapping,occams_featurecounts_batch_manifest}.csv
  and ephemeral counts/grch37_batches/OCCAMS_RNAseq_GRCh37_batch_{01..08}.txt[.summary].
Outputs: ephemeral counts/OCCAMS_RNAseq_GRCh37_raw_counts.txt[.summary].
Cache/replot: Rebuilt when SCREF_FORCE_REBUILD=TRUE; no figures or replot mode.
Run: /opt/pbs/bin/qsub analysis/OCCAMS/occams_combine_featurecount_batches.sh
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
import os
from pathlib import Path


WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
ROOT = Path("/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/counts")
BATCH_DIR = ROOT / "grch37_batches"
MAPPING = WD / "ref_outs" / "OCCAMS" / "tables" / "subject_bam_mapping.csv"
MANIFEST = WD / "ref_outs" / "OCCAMS" / "tables" / "occams_featurecounts_batch_manifest.csv"
OUT = ROOT / "OCCAMS_RNAseq_GRCh37_raw_counts.txt"


def open_counts(path):
    handle = path.open()
    for line in handle:
        if not line.startswith("#"):
            return handle, line.rstrip("\n").split("\t")
    raise ValueError(f"No header in {path}")


def main():
    with MAPPING.open(newline="") as handle:
        mapping = list(csv.DictReader(handle))
    desired_bams = [Path(row["BAM_Path"]).name for row in mapping]
    if len(desired_bams) != 282 or len(set(desired_bams)) != 282:
        raise ValueError("Expected 282 unique desired BAMs")

    paths = [BATCH_DIR / f"OCCAMS_RNAseq_GRCh37_batch_{index:02d}.txt" for index in range(1, 9)]
    opened = [open_counts(path) for path in paths]
    handles = [value[0] for value in opened]
    headers = [value[1] for value in opened]
    locations = {}
    for batch_index, header in enumerate(headers):
        for column_index, sample_path in enumerate(header[6:], start=6):
            bam = Path(sample_path).name
            if bam in locations:
                raise ValueError(f"BAM occurs in multiple batches: {bam}")
            locations[bam] = (batch_index, column_index, sample_path)
    if set(locations) != set(desired_bams):
        raise ValueError("Combined batch BAM set does not equal the 282-subject mapping")

    with OUT.open("w") as target:
        target.write("# Combined from eight size-balanced GRCh37 featureCounts batches\n")
        target.write("\t".join(headers[0][:6] + [locations[bam][2] for bam in desired_bams]) + "\n")
        n_genes = 0
        while True:
            lines = [handle.readline() for handle in handles]
            if not any(lines):
                break
            if not all(lines):
                raise ValueError("Batch count tables have unequal gene-row counts")
            fields = [line.rstrip("\n").split("\t") for line in lines]
            annotations = [value[:6] for value in fields]
            if any(value != annotations[0] for value in annotations[1:]):
                raise ValueError(f"Batch annotation mismatch at gene row {n_genes + 1}")
            counts = [fields[locations[bam][0]][locations[bam][1]] for bam in desired_bams]
            target.write("\t".join(annotations[0] + counts) + "\n")
            n_genes += 1
    for handle in handles:
        handle.close()

    summaries = []
    for path in paths:
        with Path(f"{path}.summary").open(newline="") as handle:
            summaries.append(list(csv.DictReader(handle, delimiter="\t")))
    statuses = [row["Status"] for row in summaries[0]]
    if any([row["Status"] for row in summary] != statuses for summary in summaries[1:]):
        raise ValueError("Batch summary status rows differ")
    summary_locations = {}
    for batch_index, summary in enumerate(summaries):
        for sample_path in summary[0]:
            if sample_path != "Status":
                summary_locations[Path(sample_path).name] = (batch_index, sample_path)
    with Path(f"{OUT}.summary").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Status"] + [locations[bam][2] for bam in desired_bams])
        for row_index, status in enumerate(statuses):
            writer.writerow(
                [status]
                + [summaries[summary_locations[bam][0]][row_index][summary_locations[bam][1]] for bam in desired_bams]
            )
    if n_genes <= 0:
        raise ValueError("No gene rows combined")
    print(f"genes={n_genes} samples=282 result=PASS")


if __name__ == "__main__":
    main()
