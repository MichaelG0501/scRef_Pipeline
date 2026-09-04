#!/usr/bin/env python3
####################
# scatlas_velocity_prepare_refs.py
# Status: active upstream
# Script: analysis/trajectory/scatlas_velocity_prepare_refs.py
# Methodology: analysis/methodology/trajectory/scatlas_velocity_methodology.md
# Inputs: external GRCh38-2024-A genes.gtf.gz; optional SCATLAS_REPEATMASKER_GTF,
#   otherwise UCSC hg38 rmsk.txt.gz download.
# Outputs: ref_outs/Auto_velocity_scATLAS/ref/{genes.GRCh38-2024-A.gtf,repeatmasker.hg38.gtf,rmsk.hg38.txt.gz}.
# Cache/replot: reuses non-empty references; downloads/rebuilds only missing files.
# Run: qsub analysis/trajectory/scatlas_velocity_submit.sh (reference preparation stage)
# Conda env: velocity
#
# Prepare reference GTF and RepeatMasker GTF for velocyto on scATLAS Cell
# Ranger GRCh38 BAMs.
####################

import csv
import gzip
import os
import shutil
import urllib.request
from pathlib import Path
from typing import List


WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
OUT = WD / "ref_outs" / "Auto_velocity_scATLAS"
REF_IN = Path("/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def write_gene_gtf() -> Path:
    out = OUT / "ref" / "genes.GRCh38-2024-A.gtf"
    if out.exists() and out.stat().st_size > 0:
        return out

    out.parent.mkdir(parents=True, exist_ok=True)
    with open_text(REF_IN) as inp, open(out, "w") as handle:
        shutil.copyfileobj(inp, handle)
    return out


def local_repeatmasker_candidates() -> List[Path]:
    env = os.environ.get("SCATLAS_REPEATMASKER_GTF", "")
    candidates = [Path(env)] if env else []
    candidates.extend([
        OUT / "ref" / "repeatmasker.hg38.gtf",
        Path("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_velocity_PDO/ref/repeatmasker.hg38.gtf"),
    ])
    return candidates


def download_repeatmasker() -> Path:
    raw = OUT / "ref" / "rmsk.hg38.txt.gz"
    if raw.exists() and raw.stat().st_size > 0:
        return raw

    raw.parent.mkdir(parents=True, exist_ok=True)
    tmp = raw.with_suffix(".tmp")
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
    urllib.request.urlretrieve(url, tmp)
    shutil.move(tmp, raw)
    return raw


def write_repeatmasker_gtf() -> Path:
    out = OUT / "ref" / "repeatmasker.hg38.gtf"
    for candidate in local_repeatmasker_candidates():
        if candidate.exists() and candidate.stat().st_size > 0:
            if candidate.resolve() != out.resolve():
                out.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(candidate, out)
            return out

    raw = download_repeatmasker()
    with gzip.open(raw, "rt") as inp, open(out, "w") as handle:
        handle.write("# UCSC hg38 rmsk converted to GTF for CellRanger GRCh38 BAMs\n")
        for row in csv.reader(inp, delimiter="\t"):
            if len(row) < 17:
                continue
            chrom = row[5]
            start = int(row[6]) + 1
            end = int(row[7])
            strand = row[9] if row[9] in {"+", "-"} else "."
            rep_name = row[10].replace('"', "'")
            rep_class = row[11].replace('"', "'")
            rep_family = row[12].replace('"', "'")
            attrs = (
                f'gene_id "{rep_name}"; transcript_id "{rep_name}"; '
                f'rep_class "{rep_class}"; rep_family "{rep_family}";'
            )
            handle.write(f"{chrom}\tUCSC_rmsk\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
    return out


def main() -> None:
    for subdir in ["ref", "coord", "looms", "h5ad", "figures", "tables", "logs", "tmp_sort", "tmp_pages"]:
        (OUT / subdir).mkdir(parents=True, exist_ok=True)
    gene_gtf = write_gene_gtf()
    repeat_gtf = write_repeatmasker_gtf()
    print(f"Gene GTF: {gene_gtf}")
    print(f"RepeatMasker GTF: {repeat_gtf}")


if __name__ == "__main__":
    main()
