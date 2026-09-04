#!/usr/bin/env python3
####################
# scatlas_velocyto_run.py
# Status: active upstream
# Script: analysis/trajectory/scatlas_velocyto_run.py
# Methodology: analysis/methodology/trajectory/scatlas_velocity_methodology.md
# Inputs: CLI-supplied Cell Ranger BAM, filtered-barcode file, genes GTF, and RepeatMasker mask GTF.
# Outputs: sample loom under ref_outs/Auto_velocity_scATLAS/looms/<sample>/ when invoked by the PBS wrapper.
# Cache/replot: velocyto rebuild; the wrapper stages large work in ephemeral and retains the final loom live.
# Run: analysis/trajectory/scatlas_velocity_run_velocyto.sh <sample arguments>
# Conda env: velocity
#
# Run velocyto.py on Cell Ranger BAMs, accepting standard CB/UB tags.
####################

import argparse
import logging
from pathlib import Path
from typing import Optional

import pysam
import velocyto as vcy
from velocyto.commands._run import _run


def _peek_cellranger(self, bamfile: str, lines: int = 1000) -> None:
    logging.debug(f"Peeking into {bamfile}")
    fin = pysam.AlignmentFile(bamfile)
    cellranger = 0
    parse_count = 0
    dropseq = 0
    failed = 0
    umi_tag: Optional[str] = None
    for i, read in enumerate(fin):
        if read.is_unmapped:
            continue
        if read.has_tag("CB") and read.has_tag("UB"):
            cellranger += 1
            umi_tag = "UB"
        elif read.has_tag("CB") and read.has_tag("pN"):
            parse_count += 1
            umi_tag = "pN"
        elif read.has_tag("XC") and read.has_tag("XM"):
            dropseq += 1
            umi_tag = "XM"
        else:
            logging.warning(f"Not found cell and umi barcode in entry {i} of the bam file")
            failed += 1
        if cellranger > lines:
            self.cellbarcode_str = "CB"
            self.umibarcode_str = "UB"
            break
        if parse_count > lines:
            self.cellbarcode_str = "CB"
            self.umibarcode_str = "pN"
            break
        if dropseq > lines:
            self.cellbarcode_str = "XC"
            self.umibarcode_str = "XM"
            break
        if failed > 5 * lines:
            raise IOError("The BAM does not contain usable CB+UB, CB+pN, or XC+XM barcode tags.")
    fin.close()
    if self.umibarcode_str == "NULL_UB" and umi_tag is None:
        raise IOError("No usable UMI tag detected while peeking into the BAM.")


def _peek_umi_only(self, bamfile: str, lines: int = 30) -> None:
    logging.debug(f"Peeking into {bamfile}")
    fin = pysam.AlignmentFile(bamfile)
    cellranger = 0
    parse_count = 0
    dropseq = 0
    failed = 0
    for i, read in enumerate(fin):
        if read.is_unmapped:
            continue
        if read.has_tag("UB"):
            cellranger += 1
        elif read.has_tag("pN"):
            parse_count += 1
        elif read.has_tag("XM"):
            dropseq += 1
        else:
            logging.warning(f"Not found umi barcode in entry {i} of the bam file")
            failed += 1
        if cellranger > lines:
            self.umibarcode_str = "UB"
            break
        if parse_count > lines:
            self.umibarcode_str = "pN"
            break
        if dropseq > lines:
            self.umibarcode_str = "XM"
            break
        if failed > 5 * lines:
            raise IOError("The BAM does not contain usable UB, pN, or XM UMI tags.")
    fin.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run velocyto.py on scATLAS Cell Ranger BAMs.")
    parser.add_argument("--bcfile", "-b", required=True)
    parser.add_argument("--outputfolder", "-o", required=True)
    parser.add_argument("--sampleid", "-e", required=True)
    parser.add_argument("--mask", "-m", required=True)
    parser.add_argument("--logic", "-l", default="Default")
    parser.add_argument("--samtools-threads", "-@", type=int, default=8)
    parser.add_argument("--samtools-memory", type=int, default=7000)
    parser.add_argument("--dtype", "-t", default="uint32")
    parser.add_argument("--dump", "-d", default="0")
    parser.add_argument("--verbose", "-v", action="count", default=0)
    parser.add_argument("bamfile")
    parser.add_argument("gtffile")
    args = parser.parse_args()

    Path(args.outputfolder).mkdir(parents=True, exist_ok=True)
    vcy.ExInCounter.peek = _peek_cellranger
    vcy.ExInCounter.peek_umi_only = _peek_umi_only

    _run(
        bamfile=(args.bamfile,),
        gtffile=args.gtffile,
        bcfile=args.bcfile,
        outputfolder=args.outputfolder,
        sampleid=args.sampleid,
        metadatatable=None,
        repmask=args.mask,
        onefilepercell=False,
        logic=args.logic,
        without_umi=False,
        umi_extension="no",
        multimap=False,
        test=False,
        samtools_threads=args.samtools_threads,
        samtools_memory=args.samtools_memory,
        loom_numeric_dtype=args.dtype,
        dump=args.dump,
        verbose=args.verbose,
        additional_ca={},
    )


if __name__ == "__main__":
    main()
