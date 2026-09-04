# Developmental Reference Methodology

`developmental.R` builds the persistent developmental TERM2GENE/TERM2NAME
references consumed by centred MP enrichment. Source workbooks and the adult
stomach Seurat object are external, immutable inputs under
`/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/`;
the exact files are listed in the script header.

Early-embryogenesis, normal-development long-format, adult-oesophagus, adult-
stomach, and Barrett's marker tables retain genes with adjusted P, q value, or
FDR `<0.05`. Normal-development short marker lists and the published
organogenesis lists have no per-gene significance field and are retained after
non-empty/duplicate filtering. Adult stomach markers are recalculated with
Seurat `FindAllMarkers(only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)` and
then restricted to adjusted P `<0.05`. Organogenesis is exported at both the
specified major-system and endoderm/subcluster resolutions. Term ordering is
explicit and preserved as factors before all references are combined.

The combined reference is written to
`ref_outs/developmental_reference/enrich_dev.rds`; separately replotable stage
references are written to `ref_outs/developmental_reference/per_stage/`.
Centred refinement step 04 reads only these live per-stage files. Updating any
source therefore requires rerunning `developmental.R`, followed by centred step
04 and its downstream enrichment exports. Limitations are inherited from the
heterogeneous published marker-selection methods and manual term harmonization;
the script does not treat gene-list length as equivalent evidence strength.
