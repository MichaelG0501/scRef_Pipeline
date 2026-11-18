# Status: active
# Script: analysis/GSE_Krish/geo_oac_bulk_metadata_prep.py
# Description: Clean GSE207526 clinical metadata and map exactly 110 malignant-biopsy samples.
# Methodology: analysis/methodology/GSE_Krish/gse_krish_methodology.md.
# Inputs: analysis/GSE_Krish/GSE207526_GSE112926_meta.txt;
#   ref_outs/GSE_Krish/reference/GSE207526/GSE207526_family.soft.gz.
# Outputs: ref_outs/GSE_Krish/tables/{clinical_metadata_clean,geo_sample_mapping,cohort_summary,
#   metadata_source_audit,excluded_normal_samples}.csv; reports/metadata_validation.json.
# Cache: deterministic tables are refreshed. Original data are never changed. No expression imputation.
# Run: python analysis/GSE_Krish/geo_oac_bulk_metadata_prep.py
# PBS: job_scripts/GSE_Krish/geo_oac_bulk_metadata_prep.sh; dmtcp, 1 CPU, 4 GB, 30 minutes.
import csv
import gzip
import json
from pathlib import Path

ROOT = Path('/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline')
OUT = ROOT / 'ref_outs/GSE_Krish'
SOURCE = ROOT / 'analysis/GSE_Krish/GSE207526_GSE112926_meta.txt'


def save_csv(name, rows):
    with (OUT / 'tables' / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    for tier in ['tables', 'reports']:
        (OUT / tier).mkdir(parents=True, exist_ok=True)
    with SOURCE.open(newline='') as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t') if r.get('FileName.GenomeScan')]
    clinical = []
    for r in rows:
        grade = r['res.spec.Mandard.score...']
        assert grade in ['NA', '1', '2', '3', '4', '5']
        status = r['OS.status']
        assert status in ['NA', 'Deceased', 'Alive or censored']
        clinical.append(dict(patient_id=r['FileName.GenomeScan'], sex=r['sex'],
            age_band=r['Age.Diagnosis.5yr'], mandard_score='' if grade == 'NA' else int(grade),
            mandard_response='' if grade == 'NA' else ('TRG1–3' if int(grade) <= 3 else 'TRG4–5'),
            treatment_group='Surgery only' if grade == 'NA' else 'Neoadjuvant',
            treatment_group_basis='User-defined: missing Mandard = surgery only; observed Mandard = neoadjuvant',
            os_days=int(r['OS.days']), os_event='' if status == 'NA' else int(status == 'Deceased'),
            os_status='' if status == 'NA' else status, baseline_stage_detailed=r['AJCC_baseline_detailed']))
    meta = {r['patient_id']: r for r in clinical}
    assert len(meta) == len(clinical) == 110
    samples, current = [], None
    with gzip.open(OUT / 'reference/GSE207526/GSE207526_family.soft.gz', 'rt') as f:
        for line in f:
            if line.startswith('^SAMPLE = '):
                current = dict(gsm=line.strip().split(' = ', 1)[1]); samples.append(current)
            elif line.startswith('^'):
                current = None
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                k, v = line.strip().split(' = ', 1); k = k.removeprefix('!Sample_')
                current[k] = current[k] + ' | ' + v if k in current else v
    mapped, normals = [], []
    for s in samples:
        sid = s['title'].removeprefix('sample ')
        row = dict(sample_id=sid, gsm=s['gsm'], title=s['title'], tissue=s['characteristics_ch1'])
        if sid in meta:
            assert s['characteristics_ch1'].startswith('tissue: esophageal adenocarcinoma |')
            mapped.append(row)
        else:
            assert 'healthy squamous cell tissue' in s['characteristics_ch1']; normals.append(row)
    assert len(mapped) == 110 and len(normals) == 10
    assert set(meta) == {r['sample_id'] for r in mapped}
    summaries = []
    for name in ['All cohort', 'Surgery only', 'Neoadjuvant']:
        group = clinical if name == 'All cohort' else [r for r in clinical if r['treatment_group'] == name]
        summaries.append(dict(cohort=name, selected=len(group), os_eligible=sum(r['os_event'] != '' and r['os_days'] > 0 for r in group),
                              deaths=sum(r['os_event'] == 1 for r in group)))
    audit = [dict(variable=v, source_column=c, source_file=str(SOURCE.relative_to(ROOT)), definition=d) for v, c, d in [
        ('sex', 'sex', 'Recorded sex'), ('age_band', 'Age.Diagnosis.5yr', 'Recorded five-year age band; exact age not inferred'),
        ('mandard_score', 'res.spec.Mandard.score...', 'Recorded grade 1–5; neoadjuvant clinical analyses only'),
        ('baseline_stage_detailed', 'AJCC_baseline_detailed', 'Recorded baseline stage; AJCC edition not supplied'),
        ('os_days', 'OS.days', 'Recorded OS/follow-up duration; time origin not supplied'),
        ('os_event', 'OS.status', 'Deceased=1; alive/censored=0; missing remains missing'),
        ('treatment_group', 'res.spec.Mandard.score...', 'User-defined proxy groups, not independently documented treatments')]]
    save_csv('clinical_metadata_clean.csv', clinical); save_csv('geo_sample_mapping.csv', mapped)
    save_csv('excluded_normal_samples.csv', normals); save_csv('cohort_summary.csv', summaries)
    save_csv('metadata_source_audit.csv', audit)
    with (OUT / 'reports/metadata_validation.json').open('w') as f:
        json.dump(dict(malignant_samples=110, excluded_normal_samples=10, cohorts=summaries,
                       missing_os_status_ids=[r['patient_id'] for r in clinical if r['os_event'] == '']), f, indent=2)
    print(json.dumps(summaries), flush=True)


if __name__ == '__main__':
    main()
