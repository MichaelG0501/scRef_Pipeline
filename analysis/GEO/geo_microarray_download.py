#!/usr/bin/env python3
# Status: active
# Script: analysis/GEO/geo_microarray_download.py
# Description: Download the exact 76 GSE19417 two-channel supplementary files and platform/clinical sources.
# Methodology: analysis/methodology/GEO/geo_microarray_methodology.md
# Inputs: ref_outs/geo_survival/Auto_GSE19417_meta.csv (GEO sample/channel inventory).
# Downloads: exact sample supplementary URLs in that inventory; pinned GPL4372 annotation and clinical supplement.
# Outputs: ref_outs/GEO/reference/GSE19417/*; tables/download_manifest.csv.
# Cache: existing gzip files must pass full decompression and recorded SHA256; no alternative URL or retry.
# Run/PBS: python analysis/GEO/geo_microarray_download.py; job_scripts/GEO/geo_microarray_download.sh;
#   dmtcp, 1 CPU, 4 GB, 2 hours. No plotting.
from pathlib import Path
import csv,gzip,hashlib,subprocess,concurrent.futures
ROOT=Path('/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline')
OUT=ROOT/'ref_outs/GEO';REF=OUT/'reference/GSE19417'
REF.mkdir(parents=True,exist_ok=True);(OUT/'tables').mkdir(exist_ok=True)
rows=list(csv.DictReader((ROOT/'ref_outs/geo_survival/Auto_GSE19417_meta.csv').open()))
assert len(rows)==76 and len({r['sample_geo_accession'] for r in rows})==76
sources=[(r['sample_geo_accession'],r['supplementary_file'].replace('ftp://','https://')) for r in rows]
sources += [('platform','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/annot/GPL4372.annot.gz'),
 ('clinical','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19417/suppl/GSE19417_De-identified_Data_76-patients_13Oct2010.txt.gz')]
manifest_path=OUT/'tables/download_manifest.csv'
old={r['sample']:r for r in csv.DictReader(manifest_path.open())} if manifest_path.exists() else {}
def download(item):
 sample,url=item;p=REF/url.split('/')[-1]
 if not p.exists():subprocess.run(['curl','--fail','--location','--connect-timeout','20','--max-time','600','--output',str(p),url],check=True)
 with gzip.open(p,'rb') as f:
  while f.read(1024*1024):pass
 h=hashlib.sha256(p.read_bytes()).hexdigest()
 if sample in old:assert h==old[sample]['sha256']
 return dict(sample=sample,url=url,path=str(p.relative_to(ROOT)),sha256=h,bytes=p.stat().st_size)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:manifest=list(pool.map(download,sources))
with manifest_path.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(manifest[0]));w.writeheader();w.writerows(manifest)
print('GEO_DOWNLOAD_COMPLETE:',len(manifest),'validated gzip sources',flush=True)
