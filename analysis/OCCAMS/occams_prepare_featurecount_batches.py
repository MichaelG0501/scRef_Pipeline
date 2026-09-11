#!/usr/bin/env python3
"""
Status: active
Script: analysis/OCCAMS/occams_prepare_featurecount_batches.py
Description: Greedily balance the 282 mapped OCCAMS subject BAMs by byte size
  across eight reproducible featureCounts batches.
Methodology: analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md
Inputs: ref_outs/OCCAMS/tables/subject_bam_mapping.csv and its mapped ephemeral BAMs.
Outputs: tables/ ref_outs/OCCAMS/tables/occams_featurecounts_batch_manifest.csv;
  logs/ ref_outs/OCCAMS/logs/occams_featurecounts_batch_plan_summary.txt.
Cache/replot: Deterministically rebuilt; no figures or replot mode.
Run: Called by analysis/OCCAMS/occams_submit_featurecounts_batches.sh.
Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
"""

import csv
from pathlib import Path


WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
MAPPING = WD / "ref_outs" / "OCCAMS" / "tables" / "subject_bam_mapping.csv"
MANIFEST = WD / "ref_outs" / "OCCAMS" / "tables" / "occams_featurecounts_batch_manifest.csv"
SUMMARY = WD / "ref_outs" / "OCCAMS" / "logs" / "occams_featurecounts_batch_plan_summary.txt"
N_BATCHES = 8


def main():
    with MAPPING.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 282 or len({row["Subject_ID"] for row in rows}) != 282:
        raise ValueError("Expected 282 unique subject-BAM mapping rows")
    for row in rows:
        path = Path(row["BAM_Path"])
        if not path.exists():
            raise FileNotFoundError(path)
        row["BAM_Bytes"] = path.stat().st_size

    totals = [0] * N_BATCHES
    batches = [[] for _ in range(N_BATCHES)]
    for row in sorted(rows, key=lambda value: (-value["BAM_Bytes"], value["Subject_ID"])):
        batch = min(range(N_BATCHES), key=lambda index: (totals[index], index))
        batches[batch].append(row)
        totals[batch] += row["BAM_Bytes"]

    output_rows = []
    for batch, batch_rows in enumerate(batches, start=1):
        for row in sorted(batch_rows, key=lambda value: value["Subject_ID"]):
            output_rows.append(
                {
                    "Batch": batch,
                    "Subject_ID": row["Subject_ID"],
                    "Short_ID": row["Short_ID"],
                    "BAM_Path": row["BAM_Path"],
                    "BAM_Bytes": row["BAM_Bytes"],
                }
            )
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    lines = ["OCCAMS featureCounts batch plan", "batches=8", "subjects=282"]
    for index, (batch_rows, total) in enumerate(zip(batches, totals), start=1):
        lines.append(f"batch_{index:02d}_subjects={len(batch_rows)}")
        lines.append(f"batch_{index:02d}_bytes={total}")
    lines.append("result=PASS")
    text = "\n".join(lines) + "\n"
    SUMMARY.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
