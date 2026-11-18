# Status: active
# Script: analysis/GSE_Krish/geo_oac_bulk_download.py
# Description: Download original processed GEO OAC RNA-seq matrices and complete SOFT metadata.
# Methodology: analysis/methodology/GSE_Krish/gse_krish_methodology.md; pinned downloads with gzip/SHA256 validation.
# Inputs: HTTPS GEO series GSE207526 only; Europe PMC PMC9496882 fullTextXML;
#   analysis/GSE_Krish/GSE207526_GSE112926_meta.txt (copied unchanged).
# Outputs: ref_outs/GSE_Krish/reference/GSE207526/*,
#   reference/GSE207526_GSE112926_meta.txt, reference/PMC9496882.xml,
#   reports/source_manifest.csv, reports/download_complete.json.
# Cache: existing complete files are validated and reused; empty/corrupt files cause immediate failure.
#   Existing pinned downloads are validated and reused. No replot behavior (no figures).
# Run: python analysis/GSE_Krish/geo_oac_bulk_download.py
# PBS: job_scripts/GSE_Krish/geo_oac_bulk_download.sh; 1 CPU, 4 GB, 30 min; dmtcp.
# Inventory: analysis/ANALYSIS_MAP.md and analysis/GSE_Krish/geo_oac_bulk_inventory.md.
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path('/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline')
OUT = ROOT / 'ref_outs/GSE_Krish'


def curl(*args):
    return subprocess.check_output(['curl', '--fail', '--silent', '--show-error', '--location',
                                    '--connect-timeout', '20', '--max-time', '600', *args])


def validate(path):
    if path.stat().st_size == 0:
        raise ValueError(f'Empty existing file; preservation restriction prevents replacement: {path}')
    if path.name.endswith('.gz'):
        with gzip.open(path, 'rb') as src:
            while src.read(1024 * 1024):
                pass
    elif path.suffix == '.xml':
        ET.parse(path)


PINNED = {
    'GSE207526_family.soft.gz': (6489, 'e0232b8e692c48d3931bda1cd1e0e9b979819aaa2091f95739dfee7ea5ad0c1a'),
    'GSE207526_110.EAC.and.10.Normal.for.GSEA.txt.gz': (3319962, '024cf307f5301c5382ae653acc3e0e345ec77fb4b9d7942e6c7bba60b2c34fb0'),
}


def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    expected_size, expected_sha = PINNED[path.name]
    if not path.exists():
        # Validate the completed response before creating its persistent destination.
        data = curl(url)
        assert len(data) == expected_size, f'Download size mismatch: {url}'
        assert hashlib.sha256(data).hexdigest() == expected_sha, f'Download SHA256 mismatch: {url}'
        with path.open('xb') as dst:
            dst.write(data)
    validate(path)
    assert path.stat().st_size == expected_size, f'Existing source size mismatch: {path}'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_sha, f'Existing source SHA256 mismatch: {path}'
    print('Validated', path.name, flush=True)
    return dict(url=url, file=str(path.relative_to(OUT)), bytes=expected_size, sha256=expected_sha)


def main():
    for tier in ['reference', 'tables', 'reports']:
        (OUT / tier).mkdir(parents=True, exist_ok=True)
    manifest = []
    for accession, matrix in [
        ('GSE207526', 'GSE207526_110.EAC.and.10.Normal.for.GSEA.txt.gz'),
    ]:
        base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/'
        for relative in [f'soft/{accession}_family.soft.gz', f'suppl/{matrix}']:
            manifest.append(download(base + relative + '?download=1', OUT / 'reference' / accession / Path(relative).name))
    paper_url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9496882/fullTextXML'
    paper = OUT / 'reference/PMC9496882.xml'
    if not paper.exists():
        with paper.open('xb') as dst:
            dst.write(curl(paper_url))
    validate(paper)
    manifest.append(dict(url=paper_url, file=str(paper.relative_to(OUT)), bytes=paper.stat().st_size,
                         sha256=hashlib.sha256(paper.read_bytes()).hexdigest()))
    original = ROOT / 'analysis/GSE_Krish/GSE207526_GSE112926_meta.txt'
    copied = OUT / 'reference/GSE207526_GSE112926_meta.txt'
    if not copied.exists():
        with copied.open('xb') as dst:
            dst.write(original.read_bytes())
    assert copied.read_bytes() == original.read_bytes()
    manifest.append(dict(url=str(original), file=str(copied.relative_to(OUT)), bytes=copied.stat().st_size,
                         sha256=hashlib.sha256(copied.read_bytes()).hexdigest()))
    manifest_path = OUT / 'reports/source_manifest.csv'
    if not manifest_path.exists():
        with manifest_path.open('x', newline='') as dst:
            writer = csv.DictWriter(dst, fieldnames=['url', 'file', 'bytes', 'sha256'])
            writer.writeheader()
            writer.writerows(manifest)
    else:
        with manifest_path.open(newline='') as src:
            previous = list(csv.DictReader(src))
        assert [{**x, 'bytes': str(x['bytes'])} for x in manifest] == previous
    marker = OUT / 'reports/download_complete.json'
    if not marker.exists():
        with marker.open('x') as dst:
            json.dump(dict(status='verified', files=len(manifest)), dst)
    print('ALL_DOWNLOADS_VALIDATED', flush=True)


if __name__ == '__main__':
    main()
