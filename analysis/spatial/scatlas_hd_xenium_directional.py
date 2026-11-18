#!/usr/bin/env python
####################
# Analysis registry
# Status: active
# Script: analysis/spatial/scatlas_hd_xenium_directional.py
# Description: scATLAS FDR significance followed by descriptive mean-effect
#   direction checks in both Visium HD and Xenium; no spatial significance gate.
# Methodology: analysis/methodology/spatial/scatlas_hd_xenium_directional_methodology.md
# Inputs under ref_outs/: non_malignant_mp_correlations/{01_cancer_mps_cross_only,
#   03_cancer_states_cross_only,05_cancer_mps_vs_whole_celltypes,
#   06_cancer_states_vs_whole_celltypes}/Auto_celltype_correlations_all.csv;
#   visium_hd_outs/spatial_interactions/tables/Auto_spatial_interactions_pooled.csv.gz;
#   spatial_cross_platform/tables/Auto_{within_platform_meta.csv,
#   harmonised_interactions_by_donor.csv.gz,harmonised_interactions_by_sample.csv.gz};
#   spatial_cross_platform/logs/Auto_input_provenance.json;
#   xenium_scatlas_states/spatial_neighbourhoods/tables/Auto_signature_coverage.csv;
#   Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv;
#   Yates_spatial/source_data/Zenodo_upload/processed/Xenium.h5ad (header only).
#   analysis/spatial/spatial_interaction_sample_design.tsv supplies donor identity.
# Outputs under ref_outs/spatial_scatlas_hd_xenium_directional/:
#   tables/{all_pairs,hd_confirmed_pairs,xenium_directionally_consistent_pairs,
#   donor_evidence,selection_summary,signature_coverage,figure_source,sample_figure_source}.csv;
#   figures/Auto_final/01_scatlas_hd_xenium_evidence.pdf;
#   reports/{CURRENT_RESULTS.md,validation.csv}; logs/input_manifest.json;
#   updates/new_updates/summaries/scatlas_hd_xenium_directional.csv.
# Cache/replot: SCREF_FORCE_REBUILD=TRUE recomputes this table-level screen only;
#   SCREF_REPLOT_ONLY=TRUE requires matching persisted input sizes/mtimes and method.
#   No rescoring, spatial permutations, downloads, or silent fallback.
# Run: /opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE
#   job_scripts/spatial/scatlas_hd_xenium_directional.sh
# Environment: bidcell_temp Python; dmtcp initialized before activation.
# Resources: 1 CPU, 8 GB, 00:30:00; scheduler logs in project root.
####################
import json
import os
from pathlib import Path
from textwrap import fill

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

####################
ROOT = Path('/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline')
OUT = ROOT / 'ref_outs/spatial_scatlas_hd_xenium_directional'
HD = ROOT / 'ref_outs/visium_hd_outs/spatial_interactions'
META = ROOT / 'ref_outs/spatial_cross_platform'
XENIUM = ROOT / 'ref_outs/xenium_scatlas_states/spatial_neighbourhoods'
GROUPING = ROOT / 'ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv'
H5 = ROOT / 'ref_outs/Yates_spatial/source_data/Zenodo_upload/processed/Xenium.h5ad'
DESIGN = ROOT / 'analysis/spatial/spatial_interaction_sample_design.tsv'
METHOD = '2026-11-17_scatlas_fdr_binned_hd_both_spatial_mean_directions'
MODES = {
    '01_cancer_mps_vs_tme_mps': ('01_cancer_mps_cross_only', 'Malignant MPs and TME programmes'),
    '02_cancer_states_vs_tme_mps': ('03_cancer_states_cross_only', 'Malignant states and TME programmes'),
    '03_cancer_mps_vs_whole_celltypes': ('05_cancer_mps_vs_whole_celltypes', 'Malignant MPs and TME cell types'),
    '04_cancer_states_vs_whole_celltypes': ('06_cancer_states_vs_whole_celltypes', 'Malignant states and TME cell types'),
}
STATES = ['Classic proliferation', 'Squamous-to-intestinal', 'Glandular-to-intestinal',
          'Stress-adaptive', 'Cancer-cell immune mimicry']
WHOLE = {'fibroblast': 'fibroblast', 'endothelial': 'endothelial', 'macrophage': 'macrophage',
         'nk': 'nk.cell', 'plasma': 'plasma', 'cd4': 't.cell', 'cd8': 't.cell'}
KEYS = ['mode_id', 'source_feature', 'target_feature']
RAW = {mode: ROOT / 'ref_outs/non_malignant_mp_correlations' / subdir / 'Auto_celltype_correlations_all.csv'
       for mode, (subdir, _) in MODES.items()}
INPUTS = list(RAW.values()) + [GROUPING, H5, DESIGN,
    HD / 'tables/Auto_spatial_interactions_pooled.csv.gz',
    HD / 'tables/Auto_spatial_interaction_parameters.csv',
    META / 'tables/Auto_within_platform_meta.csv',
    META / 'tables/Auto_harmonised_interactions_by_sample.csv.gz',
    META / 'tables/Auto_harmonised_interactions_by_donor.csv.gz',
    META / 'logs/Auto_input_provenance.json', XENIUM / 'tables/Auto_signature_coverage.csv']


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    p = values[ok]
    if np.any((p < 0) | (p > 1)):
        raise ValueError('Invalid probability')
    order = np.argsort(p)
    if len(p):
        adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
        inverse = np.empty(len(p))
        inverse[order] = np.minimum(adjusted, 1)
        result[ok] = inverse
    return result


def manifest():
    return {str(p.relative_to(ROOT)): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
            for p in INPUTS}


def check_upstream():
    hd_params = pd.read_csv(HD / 'tables/Auto_spatial_interaction_parameters.csv').set_index('parameter')['value']
    if float(hd_params['bin_size_um']) != 16:
        raise ValueError('The current directional screen requires 16-um binned HD interactions')
    recorded = json.loads((META / 'logs/Auto_input_provenance.json').read_text())
    if recorded['version'] != '2026-10-27_donor_aware_separate_tile_origins':
        raise ValueError('Requires current spatial inference version')
    for name, old in recorded['inputs'].items():
        path = Path(name)
        if path.suffix == '.py':
            continue  # Plot-only edits to the producer do not change its inference version.
        stat = path.stat()
        if stat.st_size != old['size_bytes'] or stat.st_mtime_ns != old['mtime_ns']:
            raise ValueError(f'Upstream meta-analysis is stale: {path}')
    design = pd.read_csv(DESIGN, sep='\t')
    if design.groupby('platform')['donor_id'].nunique().to_dict() != {'visium_hd': 2, 'xenium': 6}:
        raise ValueError('Expected two HD and six Xenium donors')


def read_atlas(grouping):
    aliases = dict(zip(grouping['description'], grouping['mp']))
    aliases.update({mp: mp for mp in grouping['mp']})
    aliases.update({'replication-stress-associated cell cycling': 'MP13+',
                    'Stress-reactive columnar epithelium': 'MP6+',
                    'Epithelial antiviral interferon response': 'MP11+',
                    'Immune-interactive glandular progenitor': 'MP17',
                    'Hypoxic inflammatory adaptive plasticity': 'MP12',
                    'T/NK-like cancer-cell immune mimicry': 'MP15'})
    state_aliases = dict(zip(STATES, STATES))
    state_aliases.update({'Basal to intestinal metaplasia': STATES[1],
                          'SMG to intestinal metaplasia': STATES[2], 'Stress adaptive': STATES[3]})
    rows = []
    for mode, path in RAW.items():
        raw = pd.read_csv(path)
        raw = raw.loc[(raw['compartment1'] == 'cancer') ^ (raw['compartment2'] == 'cancer')].copy()
        if raw.empty or raw['edge_id'].duplicated().any():
            raise ValueError(f'Empty or duplicated cancer-TME hypotheses: {mode}')
        cancer_first = raw['compartment1'] == 'cancer'
        source = raw['mp1_name'].where(cancer_first, raw['mp2_name'])
        target = raw['mp2_name'].where(cancer_first, raw['mp1_name'])
        parent = raw['celltype2_display'].where(cancer_first, raw['celltype1_display'])
        is_mp = mode.startswith(('01_', '03_'))
        whole = mode.startswith(('03_', '04_'))
        source = source.map(aliases if is_mp else state_aliases)
        spatial_parent = parent.map(WHOLE) if whole else parent
        if source.isna().any() or spatial_parent.isna().any():
            raise ValueError(f'Unknown source or target label in {mode}')
        rows.append(pd.DataFrame({
            'mode_id': mode, 'scrna_edge_id': raw['edge_id'], 'source_feature': source,
            'source_kind': 'cancer_mp' if is_mp else 'cancer_state',
            'scrna_target_compartment': parent, 'target_compartment': spatial_parent,
            'target_kind': 'whole_celltype' if whole else 'tme_mp',
            'target_feature': spatial_parent if whole else spatial_parent + '__' + target,
            'scatlas_rho': raw['spearman_r'], 'scatlas_p': raw['spearman_p'],
            'scatlas_fdr': bh(raw['spearman_p']), 'scatlas_shared_samples': raw['shared_sample_n'],
            'scatlas_family_n': int(np.isfinite(raw['spearman_p']).sum()),
        }))
    return pd.concat(rows, ignore_index=True)


def build(grouping):
    atlas = read_atlas(grouping)
    hd = pd.read_csv(HD / 'tables/Auto_spatial_interactions_pooled.csv.gz')
    hd = hd.loc[hd['ring'] == 1].copy()
    for _, page in hd.groupby('mode_id'):
        if not np.allclose(bh(page['pooled_p']), page['pooled_fdr'], equal_nan=True):
            raise ValueError('HD full-family FDR differs from stored pooled P values')
    hd_columns = {'pooled_log2_enrichment': 'hd_mean_log2_enrichment', 'pooled_z': 'hd_z',
                  'pooled_p': 'hd_p', 'pooled_fdr': 'hd_fdr', 'n_donors_tested': 'hd_n_donors',
                  'direction_concordant': 'hd_all_donors_direction_concordant',
                  'n_discordant_paired_donors': 'hd_n_discordant_paired_donors'}
    atlas = atlas.merge(hd[KEYS + list(hd_columns)].rename(columns=hd_columns), on=KEYS,
                        how='left', validate='many_to_one')
    within = pd.read_csv(META / 'tables/Auto_within_platform_meta.csv')
    x = within.loc[(within['platform'] == 'xenium') & (within['ring'] == 1)].copy()
    x_columns = {'mean_log2_enrichment': 'xenium_mean_log2_enrichment', 'pooled_z': 'xenium_z',
                 'pooled_p': 'xenium_p', 'pooled_fdr': 'xenium_fdr', 'n_donors_tested': 'xenium_n_donors',
                 'n_enriched_donors': 'xenium_n_enriched', 'n_depleted_donors': 'xenium_n_depleted'}
    atlas = atlas.merge(x[KEYS + list(x_columns)].rename(columns=x_columns), on=KEYS,
                        how='left', validate='many_to_one')
    atlas['scatlas_significant'] = atlas['scatlas_fdr'] < .05
    atlas['hd_significant'] = (atlas['hd_fdr'] < .05) & (atlas['hd_n_donors'] >= 2)
    atlas['hd_direction_available'] = (atlas['hd_n_donors'] >= 2) & np.isfinite(atlas['hd_mean_log2_enrichment'])
    atlas['scatlas_hd_direction_match'] = atlas['hd_direction_available'] & (atlas['scatlas_rho'] * atlas['hd_mean_log2_enrichment'] > 0)
    atlas['hd_z_matches_mean'] = atlas['hd_z'] * atlas['hd_mean_log2_enrichment'] > 0
    atlas['hd_directionally_consistent'] = atlas['scatlas_significant'] & atlas['scatlas_hd_direction_match']
    # Retain this stricter diagnostic export; it does not determine final hits.
    atlas['hd_confirmed'] = atlas['hd_directionally_consistent'] & atlas['hd_significant'] & atlas['hd_z_matches_mean']
    atlas['xenium_direction_available'] = (atlas['xenium_n_donors'] >= 2) & np.isfinite(atlas['xenium_mean_log2_enrichment'])
    atlas['xenium_same_direction'] = atlas['xenium_direction_available'] & (atlas['scatlas_rho'] * atlas['xenium_mean_log2_enrichment'] > 0)
    atlas['passes_directional_rule'] = atlas['hd_directionally_consistent'] & atlas['xenium_same_direction']
    atlas['xenium_significant_same_direction'] = atlas['xenium_same_direction'] & (atlas['xenium_fdr'] < .05) & (atlas['scatlas_rho'] * atlas['xenium_z'] > 0)
    atlas['xenium_status'] = np.select([
        ~atlas['xenium_direction_available'], atlas['xenium_same_direction'],
        atlas['xenium_mean_log2_enrichment'] == 0], ['Unavailable', 'Same direction', 'Zero effect'], default='Opposite direction')
    atlas['xenium_unavailable_reason'] = np.where(atlas['xenium_direction_available'], '',
        np.where(atlas['target_compartment'].isin(['nk', 'nk.cell', 'plasma']),
                 'No matched Xenium parent annotation', 'Fewer than two eligible Xenium donors'))
    atlas['xenium_donors_matching_scatlas'] = np.where(atlas['scatlas_rho'] > 0, atlas['xenium_n_enriched'], atlas['xenium_n_depleted'])
    atlas['xenium_z_matches_mean'] = atlas['xenium_z'] * atlas['xenium_mean_log2_enrichment'] > 0
    labels = grouping.set_index('mp')['plot_label'].to_dict()
    atlas['source_label'] = [labels[s] if kind == 'cancer_mp' else s
                             for s, kind in zip(atlas['source_feature'], atlas['source_kind'])]
    names = {'cd4': 'CD4 T', 'cd8': 'CD8 T', 'nk': 'NK', 'nk.cell': 'NK cells', 't.cell': 'T cells',
             'fibroblast': 'Fibroblast', 'endothelial': 'Endothelial', 'macrophage': 'Macrophage', 'plasma': 'Plasma'}
    atlas['target_label'] = [names[parent] + (' · ' + target.split('__')[1] if kind == 'tme_mp' else '')
        for parent, target, kind in zip(atlas['target_compartment'], atlas['target_feature'], atlas['target_kind'])]
    atlas['subtype_resolution_note'] = np.where((atlas['target_kind'] == 'whole_celltype') &
        atlas['scrna_target_compartment'].isin(['cd4', 'cd8']), 'Whole T-cell spatial evidence; CD4/CD8 specificity unresolved', '')
    atlas['pair_label'] = atlas['source_label'] + ' / ' + atlas['target_label']
    coverage = pd.read_csv(XENIUM / 'tables/Auto_signature_coverage.csv').set_index('feature_id')
    atlas['xenium_source_signature_coverage'] = [coverage.loc['cancer__' + s, 'coverage_fraction'] if kind == 'cancer_mp' else np.nan
        for s, kind in zip(atlas['source_feature'], atlas['source_kind'])]
    atlas['xenium_target_signature_coverage'] = [coverage.loc[t, 'coverage_fraction'] if kind == 'tme_mp' and t in coverage.index else np.nan
        for t, kind in zip(atlas['target_feature'], atlas['target_kind'])]
    order = {name: i for i, name in enumerate(grouping['mp'].tolist() + STATES)}
    atlas['source_order'] = atlas['source_feature'].map(order)
    atlas = atlas.sort_values(['mode_id', 'source_order', 'target_feature', 'scrna_edge_id'])
    donors = pd.read_csv(META / 'tables/Auto_harmonised_interactions_by_donor.csv.gz')
    donors = donors.loc[donors['ring'] == 1]
    if donors.duplicated(KEYS + ['platform', 'donor_id']).any():
        raise ValueError('Duplicate spatial hypothesis within donor')
    # CD4/CD8 whole-type atlas hypotheses deliberately share a spatial T-cell test.
    selected_donors = atlas.loc[atlas['scatlas_significant'], ['scrna_edge_id'] + KEYS].merge(
        donors, on=KEYS, how='left', validate='many_to_many')
    # Reconstruct the Xenium mean used for the descriptive direction check.
    means = donors.loc[donors['platform'] == 'xenium'].groupby(KEYS)['log2_enrichment'].mean().rename('reconstructed_mean')
    check = atlas.merge(means, on=KEYS, how='left', validate='many_to_one')
    if not np.allclose(check['xenium_mean_log2_enrichment'], check['reconstructed_mean'], equal_nan=True):
        raise ValueError('Xenium mean does not match independent-donor source table')
    atlas.to_csv(OUT / 'tables/all_pairs.csv', index=False)
    atlas.loc[atlas['hd_confirmed']].to_csv(OUT / 'tables/hd_confirmed_pairs.csv', index=False)
    atlas.loc[atlas['passes_directional_rule']].to_csv(OUT / 'tables/xenium_directionally_consistent_pairs.csv', index=False)
    selected_donors.to_csv(OUT / 'tables/donor_evidence.csv', index=False)
    coverage.reset_index().to_csv(OUT / 'tables/signature_coverage.csv', index=False)
    return atlas


def plot(atlas):
    """Two landscape pages per populated mode: platform means and section effects."""
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.titlesize': 13,
                         'axes.labelsize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
                         'pdf.fonttype': 42})
    colours = {'Enrichment': '#B2182B', 'Depletion': '#2166AC'}
    selected = atlas.loc[atlas['passes_directional_rule']].copy()
    selected['direction'] = np.where(selected['scatlas_rho'] > 0, 'Enrichment', 'Depletion')
    selected['direction_order'] = np.where(selected['scatlas_rho'] > 0, 0, 1)
    selected = selected.sort_values(['mode_id', 'direction_order', 'source_order', 'target_feature', 'scrna_edge_id'])
    selected['figure_order'] = selected.groupby('mode_id').cumcount()
    source_abbreviations = {
        'MP3+ Basal-columnar invasive epithelium': 'MP3+ Basal-columnar invasive epi.',
        'MP11+ Epithelial type I interferon response': 'MP11+ Epithelial IFN-I response',
    }
    target_abbreviations = {'Fibroblast': 'Fib.', 'Macrophage': 'Mac.', 'Endothelial': 'Endo.'}
    selected['pair_axis_label'] = [source_abbreviations.get(row.source_label, row.source_label) + ' / ' +
        target_abbreviations.get(row.target_label.split(' · ')[0], row.target_label.split(' · ')[0]) +
        (' ' + row.target_label.split(' · ')[1] if ' · ' in row.target_label else '') +
        (' [scATLAS ' + row.scrna_target_compartment.upper() + ']' if row.subtype_resolution_note else '')
        for row in selected.itertuples()]
    selected.to_csv(OUT / 'tables/figure_source.csv', index=False)
    sections = pd.read_csv(META / 'tables/Auto_harmonised_interactions_by_sample.csv.gz')
    sections = sections.loc[sections['ring'] == 1].copy()
    if sections.duplicated(KEYS + ['sample_id']).any():
        raise ValueError('Duplicate interaction within spatial section')
    eligible = sections['eligible_pair'].astype(str).str.lower()
    if not eligible.isin(['true', 'false']).all():
        raise ValueError('Non-boolean section eligibility')
    sections['display_log2_enrichment'] = sections['log2_enrichment'].where(eligible == 'true')
    section_source = selected[['scrna_edge_id', 'figure_order', 'direction'] + KEYS].merge(
        sections, on=KEYS, how='left', validate='many_to_many')
    section_source.to_csv(OUT / 'tables/sample_figure_source.csv', index=False)
    with PdfPages(OUT / 'figures/Auto_final/01_scatlas_hd_xenium_evidence.pdf') as pdf:
        if selected.empty:
            fig, ax = plt.subplots(figsize=(300 / 25.4, 100 / 25.4))
            ax.axis('off')
            ax.text(.5, .6, 'No scATLAS-significant pairs have matching trends in both spatial datasets', ha='center')
            pdf.savefig(fig)
            plt.close(fig)
        for mode, page in selected.groupby('mode_id', sort=True):
            height = 42 + 10 * len(page)
            labels = page['pair_axis_label'].tolist()
            row_colours = [colours[value] for value in page['direction']]
            split = int((page['direction'] == 'Enrichment').sum()) - .5
            fig, axes = plt.subplots(1, 3, figsize=(300 / 25.4, height / 25.4), sharey=True)
            fig.subplots_adjust(left=.45, right=.97, top=1 - 20 / height, bottom=22 / height, wspace=.35)
            spatial_limit = max(.15, np.nanmax(np.abs(page[['hd_mean_log2_enrichment', 'xenium_mean_log2_enrichment']].to_numpy())) * 1.2)
            for col, (field, title, xlabel) in enumerate([
                ('scatlas_rho', 'scATLAS', 'Spearman ρ'), ('hd_mean_log2_enrichment', 'Visium HD', 'Mean log₂ enrichment'),
                ('xenium_mean_log2_enrichment', 'Xenium', 'Mean log₂ enrichment')]):
                ax = axes[col]
                ax.axvline(0, color='#999999', lw=.8)
                ax.scatter(page[field], range(len(page)), s=38, c=row_colours, zorder=3)
                if -.5 < split < len(page) - .5:
                    ax.axhline(split, color='#999999', lw=.7, linestyle='--')
                ax.set_title(title, fontweight='bold')
                ax.set_xlabel(xlabel)
                ax.set_xlim((-1, 1) if col == 0 else (-spatial_limit, spatial_limit))
                ax.set_ylim(len(page) - .5, -.5)
                ax.set_yticks(range(len(page)), labels)
                ax.tick_params(axis='y', length=0, pad=8)
                ax.spines[['top', 'right', 'left']].set_visible(False)
                if col:
                    ax.tick_params(labelleft=False)
            fig.suptitle(MODES[mode][1], x=.025, ha='left', y=1 - 3 / height, fontsize=15, fontweight='bold')
            handles = [Line2D([], [], marker='o', ls='none', color=colour, markersize=6, label=direction)
                       for direction, colour in colours.items()]
            fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.65, .01), frameon=False, ncol=2, fontsize=12)
            pdf.savefig(fig)
            plt.close(fig)

            sample_page = section_source.loc[section_source['mode_id'] == mode]
            sample_limit = max(.5, np.ceil(np.nanmax(np.abs(sample_page['display_log2_enrichment'])) * 2) / 2)
            fig, axes = plt.subplots(1, 3, figsize=(300 / 25.4, height / 25.4), sharey=True)
            fig.subplots_adjust(left=.45, right=.97, top=1 - 20 / height, bottom=25 / height, wspace=.35)
            atlas_ax = axes[0]
            atlas_ax.axvline(0, color='#999999', lw=.8)
            atlas_ax.scatter(page['scatlas_rho'], range(len(page)), s=38, c=row_colours, zorder=3)
            if -.5 < split < len(page) - .5:
                atlas_ax.axhline(split, color='#999999', lw=.7, linestyle='--')
            atlas_ax.set_title('scATLAS', fontweight='bold')
            atlas_ax.set_xlabel('Spearman ρ')
            atlas_ax.set_xlim(-1, 1)
            atlas_ax.set_ylim(len(page) - .5, -.5)
            atlas_ax.set_yticks(range(len(page)), labels)
            atlas_ax.tick_params(axis='y', length=0, pad=8)
            atlas_ax.spines[['top', 'right', 'left']].set_visible(False)
            for ax, platform, title, mean_field in zip(axes[1:], ['visium_hd', 'xenium'],
                    ['Visium HD', 'Xenium'], ['hd_mean_log2_enrichment', 'xenium_mean_log2_enrichment']):
                ax.axvline(0, color='#999999', lw=.8)
                for y, row in enumerate(page.itertuples()):
                    points = sample_page.loc[(sample_page['scrna_edge_id'] == row.scrna_edge_id) &
                        (sample_page['platform'] == platform), 'display_log2_enrichment'].dropna()
                    ax.scatter(points, np.full(len(points), y), s=15, color='#777777', alpha=.8, zorder=3)
                    ax.scatter(getattr(row, mean_field), y, s=65, marker='D', color=colours[row.direction],
                               edgecolors='white', linewidths=.7, zorder=4)
                if -.5 < split < len(page) - .5:
                    ax.axhline(split, color='#999999', lw=.7, linestyle='--')
                ax.set_title(title, fontweight='bold')
                ax.set_xlim(-sample_limit * 1.08, sample_limit * 1.08)
                ax.set_ylim(len(page) - .5, -.5)
                ax.set_yticks(range(len(page)), labels)
                ax.tick_params(axis='y', length=0, pad=8)
                ax.set_xlabel('log₂ enrichment')
                ax.spines[['top', 'right', 'left']].set_visible(False)
            for ax in axes[1:]:
                ax.tick_params(labelleft=False)
            fig.suptitle('scATLAS and individual spatial samples', x=.025,
                         y=1 - 3 / height, ha='left', fontsize=15, fontweight='bold')
            handles = [Line2D([], [], marker='o', ls='none', color='#777777', markersize=4,
                              label='Individual sample')]
            handles += [Line2D([], [], marker='D', ls='none', color=colour, markersize=7,
                               label='Mean: ' + direction.lower()) for direction, colour in colours.items()]
            fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.64, .025),
                       frameon=False, ncol=3, fontsize=12)
            fig.text(.45, .005, 'HD mean gives FFPEA1/FFPED1 one donor weight.', fontsize=12)
            pdf.savefig(fig)
            plt.close(fig)


def report(atlas, panel_genes):
    rows = []
    for mode, page in atlas.groupby('mode_id', sort=True):
        significant = page.loc[page['scatlas_significant']]
        hd_match = significant.loc[significant['scatlas_hd_direction_match']]
        passed = page.loc[page['passes_directional_rule']]
        rows.append({'mode_id': mode, 'scatlas_pairs': len(page), 'scatlas_fdr_0_05': len(significant),
            'hd_same_direction': len(hd_match),
            'hd_unavailable': int((~significant['hd_direction_available']).sum()),
            'both_spatial_same_direction': len(passed),
            'same_direction_unique_spatial_tests': len(passed[KEYS].drop_duplicates()),
            'xenium_opposite_after_hd_match': int((hd_match['xenium_status'] == 'Opposite direction').sum()),
            'xenium_unavailable_after_hd_match': int((hd_match['xenium_status'] == 'Unavailable').sum()),
            'final_hits_with_hd_fdr_0_05': int(passed['hd_significant'].sum()),
            'final_hits_with_xenium_fdr_0_05': int(passed['xenium_significant_same_direction'].sum())})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / 'tables/selection_summary.csv', index=False)
    summary_path = ROOT / 'updates/new_updates/summaries/scatlas_hd_xenium_directional.csv'
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    if summary_path.stat().st_size >= 100000:
        raise ValueError('Compact summary exceeds 100 KB')
    lines = ['# scATLAS-significant interactions with concordant spatial trends', '',
        'The final hit rule requires scATLAS Spearman BH FDR < 0.05 and matching nonzero mean donor '
        'log2-enrichment directions in both Visium HD and Xenium. Each spatial platform must have at least '
        'two eligible independent donors. Neither spatial P/FDR nor spatial Z direction is a selection gate. '
        'Spatial P/FDR, inferential directions and donor heterogeneity remain available as diagnostics. '
        'These are scATLAS-significant, spatially concordant hits; direction agreement does not establish '
        'statistical significance in either spatial dataset.', '',
        'scATLAS BH correction covers every finite cancer–TME test within each of four comparison modes, '
        'before spatial filtering. It is not recomputed on successful hits. Only the primary ±16-µm '
        'axis-distance window is used. No favourable radius, section or donor is selected. '
        'The three HD sections enter through two donor weights, preserving the FFPEA1/FFPED1 pairing. '
        'Means are calculated across eligible donors; agreement does not require all individual donors '
        'to agree. Missing evidence and exactly zero effects cannot pass.', '',
        'This secondary exploratory rule was revised after examining earlier results. Its final list '
        'does not carry a formal spatial replication FDR guarantee. Near-zero means can pass a sign-only '
        'check; effect magnitudes, donor directions and signature coverage should guide interpretation. '
        'scATLAS sample-level correlations and local spatial enrichment measure different quantities.', '',
        f'The Xenium matrix contains {panel_genes:,} genes. Its six donor fields include 16 tissue regions. '
        'NK and plasma lack matched Xenium parent annotations. Whole CD4/CD8 hypotheses map to a shared '
        'spatial T-cell test; both original hypotheses and unique spatial-test counts are retained. '
        'All malignant MP labels use the current full descriptions. TME programmes retain compartment-specific IDs.', '',
        f'Final hits: {int(atlas.passes_directional_rule.sum())}. '
        f'Of these, {int((atlas.passes_directional_rule & atlas.hd_significant).sum())} also have HD FDR < 0.05 '
        f'and {int((atlas.passes_directional_rule & atlas.xenium_significant_same_direction).sum())} '
        'have Xenium FDR < 0.05 with a matching inferential direction.', '',
        '## Current figures and source data', '',
        '- [Final hit figure](../figures/Auto_final/01_scatlas_hd_xenium_evidence.pdf): all final hits, grouped by enrichment/depletion; platform means followed by scATLAS correlations alongside small individual-sample points and the same spatial donor-weighted means for each populated mode. Labels are black and one line per pair; abbreviations retain full names in the source table.',
        '- [Final hit table](../tables/xenium_directionally_consistent_pairs.csv)',
        '- [Complete hypotheses, failures and unavailable evidence](../tables/all_pairs.csv)',
        '- [Stricter HD-significance diagnostic subset](../tables/hd_confirmed_pairs.csv): not the final selection criterion.',
        '- [Selection summary](../tables/selection_summary.csv)',
        '- [Donor evidence for every scATLAS-significant pair](../tables/donor_evidence.csv)',
        '- [Signature coverage](../tables/signature_coverage.csv)',
        '- [Figure source and row order](../tables/figure_source.csv)',
        '- [Individual sample figure source](../tables/sample_figure_source.csv)',
        '- [Production validation](validation.csv)', '',
        '## Counts', '', '```', summary.to_string(index=False), '```', '']
    (OUT / 'reports/CURRENT_RESULTS.md').write_text('\n'.join(lines))
    print(summary.to_string(index=False), flush=True)
    print(atlas.loc[atlas['passes_directional_rule'], ['pair_label', 'scatlas_fdr', 'hd_fdr', 'xenium_mean_log2_enrichment', 'xenium_status']].to_string(index=False), flush=True)


def main():
    if not os.getenv('PBS_JOBID'):
        raise ValueError('Run within a PBS allocation')
    flags = [os.getenv(name, 'FALSE') for name in ['SCREF_FORCE_REBUILD', 'SCREF_REPLOT_ONLY']]
    if any(flag not in ['TRUE', 'FALSE'] for flag in flags) or flags.count('TRUE') != 1:
        raise ValueError('Set exactly one rebuild or replot flag to TRUE')
    for tier in ['tables', 'reports', 'logs', 'figures/Auto_final']:
        (OUT / tier).mkdir(parents=True, exist_ok=True)
    check_upstream()
    current = manifest()
    grouping = pd.read_csv(GROUPING)
    if len(grouping) != 17 or grouping['mp'].duplicated().any() or grouping['plot_label'].isna().any():
        raise ValueError('Expected the current 17-MP label collection')
    with h5py.File(H5) as handle:
        panel_genes = int(handle['layers/counts'].attrs['shape'][1])
    if flags[1] == 'TRUE':
        saved = json.loads((OUT / 'logs/input_manifest.json').read_text())
        if saved['inputs'] != current or saved['method'] != METHOD:
            raise ValueError('Inputs changed: explicit screen rebuild required')
        atlas = pd.read_csv(OUT / 'tables/all_pairs.csv').fillna({'subtype_resolution_note': ''})
    else:
        atlas = build(grouping)
        (OUT / 'logs/input_manifest.json').write_text(json.dumps({'method': METHOD, 'inputs': current, 'panel_genes': panel_genes}, indent=2) + '\n')
    plot(atlas)
    report(atlas, panel_genes)
    if manifest() != current:
        raise ValueError('An input changed during this run')
    checks = pd.DataFrame({'check': ['Current upstream provenance', 'Two HD and six Xenium donors',
        'Full-family BH reconstruction', 'Independent-donor Xenium means', 'Inputs unchanged during analysis',
        'Unique scATLAS edge identities'], 'passed': True})
    if atlas.duplicated(['mode_id', 'scrna_edge_id']).any():
        raise ValueError('Duplicated atlas edge')
    checks.to_csv(OUT / 'reports/validation.csv', index=False)
    print('Completed:', OUT, flush=True)


if __name__ == '__main__':
    main()
####################
