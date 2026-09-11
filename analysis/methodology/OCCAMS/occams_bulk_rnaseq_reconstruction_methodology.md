# OCCAMS bulk RNA-seq reconstruction methodology

## Scope and run order

This workflow reconstructs a subject-level raw gene-count matrix from controlled-access OCCAMS BAMs. The operational order is:

1. download EGA dataset metadata and match EGA file accessions to the OCCAMS SHA subject identifiers;
2. select `Library_Strategy == "RNA-Seq"` and `Phenotype == "tumor"`;
3. download BAM payloads to ephemeral storage;
4. audit metadata coverage and BAM reference headers;
5. symlink single-file subjects or merge multiple BAMs belonging to the same subject;
6. quantify exon-overlapping fragments with featureCounts and build the canonical live matrix.

Raw BAMs and the reconstructable featureCounts text file remain under `/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/`. The GTF, final compressed matrix, clinical metadata snapshot, EGA mappings, subject mapping, QC table, and run summaries are canonical under `ref_outs/OCCAMS/`.

## Credentials and controlled data

The EGA credential JSON is external to the Git repository at `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/ega.json`, with owner-only mode `600`. Download scripts reference that exact path and never embed, print, or copy credential values into `scRef_Pipeline`. Raw controlled-access BAMs are not copied to live or Git.

## Metadata matching and filtering

The master clinical file is `OCCAMS_metadata.csv`; its `SHA IDs` field is the subject identifier. EGA sample records are scanned for exact full-field matches to this identifier. Sample-to-file and study/experiment/run/sample tables then provide EGAF and library-strategy links. Phenotype is normalized from the EGA phenotype field, falling back to title/description only when necessary. Exact RNA-seq/tumour labels define the count cohort.

An EGAF duplicated across EGA datasets is accepted only when subject, library strategy, and phenotype agree; the lexicographically greatest dataset accession is retained to reproduce the original deterministic resolution. Metadata rows are not collapsed: the clinical file has repeated SHA identifiers, so all rows for count-matrix subjects are preserved in the live metadata-row table. Downstream analyses must define a clinical-row aggregation rule appropriate to their question rather than silently selecting one row.

The audit unit is first the downloaded EGAF payload and then the unique subject. Production expectations are 300 downloaded BAM payloads, 302 authorized RNA-seq/tumour EGAF rows, and 282 subjects. EGAF00004939962 and EGAF00004939963 returned EGA 403 responses in the August 2026 download, but both belong to subject `7f03ba2fd9cadcebb17de87c65d1e2b005a844d58549fc1b161c510e9f811b05`, which is represented by other downloaded BAMs. A passing audit requires every downloaded BAM to map to exactly one target row and metadata subject, and every one of the 282 target subjects to have at least one downloaded BAM.

## BAM preparation and genome-build validation

Single-file subjects are symlinked into the ephemeral `merged_bams/` directory. Same-subject multi-file BAMs are merged with SAMtools and coordinate sorted/indexed where required. The 16-character output stem is retained for compatibility with existing files only after checking that it is unique across all 282 full SHA identifiers; the full identifier is preserved in `subject_bam_mapping.csv` and in the final matrix header.

Genome build is inferred from the primary chromosome 1 sequence length in each raw BAM header: 249,250,621 is GRCh37 and 248,956,422 is GRCh38. The production audit found all 300 BAMs to be GRCh37. Of these, 156 use UCSC-style `chr` names and 144 omit `chr`; none have an unknown build. Three subjects combine GRCh37 BAMs from both naming families. This is not a coordinate-build conflict: SAMtools preserves both reference names, and featureCounts resolves them with the explicit `grch37_contig_aliases.csv` mapping. FASTQ reconstruction or realignment is therefore unnecessary for the current files.

## Quantification

The pinned annotation is the GENCODE v19 comprehensive chromosome annotation for GRCh37.p13, downloaded from the official GENCODE/EMBL-EBI release directory. Its required MD5 is `bd83e28270e595d3bde6bfcb21c9748f`. The 282 BAMs are greedily assigned by byte size to eight batches to reduce the walltime critical path; columns are then recombined in the canonical subject-mapping order with exact annotation-row equality checks. featureCounts 2.0.6 uses:

- eight threads;
- paired-end fragment counting (`-p --countReadPairs`);
- unstranded assignment (`-s 0`), which is conservative across the multiple contributing OCCAMS sequencing series;
- exon features (`-t exon`) summarized to `gene_id` (`-g gene_id`);
- the tracked `chr1,1` through `chr22,22`, `chrX,X`, `chrY,Y`, and `chrM,MT` aliases;
- default exclusion of multimapping and ambiguous fragments.

The assignment/statistical unit is a read pair (fragment), and the matrix unit is integer fragments per GENCODE gene per unique OCCAMS subject. No expression-level filtering, normalization, batch correction, or clinical-row aggregation is applied at reconstruction.

Before the full run, representative `chr`, non-`chr`, and mixed-name BAMs are quantified. Validation reports assigned, no-feature, multimapping, and annotation-compatible (`Assigned + Unassigned_MultiMapping`) percentages. The preflight passes when every test BAM has at least 50% annotation-compatible alignments and at most 40% `Unassigned_NoFeatures`. Multimapping is reported separately because it is an input-alignment property, not evidence of assembly mismatch.

## Cache semantics and outputs

The reference download and batch quantification reuse complete existing outputs. `SCREF_FORCE_REBUILD=TRUE` forces regeneration. The per-batch and combined raw featureCounts tables are ephemeral, fully reconstructable same-workflow caches. The canonical downstream matrix is the gzipped live TSV with full SHA subject IDs. Persistent QC retains every featureCounts status/count/percentage per subject. Replot-only mode is not applicable because this workflow makes no figures.

## Limitations and validation

The supplied BAMs were produced by multiple studies and aligners. Multimapping fractions vary substantially and are intentionally not rescued into gene counts. Three same-subject merges contain both GRCh37 contig naming conventions; the explicit alias test guards this edge case. Header-based assembly inference verifies reference coordinate length but cannot recover every aligner parameter or library-preparation detail. One retained subject (`90687ad9d09f7677ab6fc83f75b78bd3fed8a9934e071f071306007c79277dab`, EGAF00001809022) has only 30.11% annotation-compatible alignments and 69.34% `Unassigned_NoFeatures` despite a verified GRCh37 header. It is flagged in the sample QC summary and should be tested in sensitivity analyses rather than silently removed from the complete cohort. Clinical metadata contain repeated subject identifiers and require analysis-specific aggregation downstream.

A production run is complete only when PBS exit status is zero; the raw table and summary each contain 282 sample columns; all columns map uniquely to 282 full SHA identifiers; all identifiers occur in the metadata; the persistent compressed matrix and QC/source tables are non-empty; and the compact run summaries report `PASS`.
