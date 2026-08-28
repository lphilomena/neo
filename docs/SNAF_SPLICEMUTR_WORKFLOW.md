# SNAF and SpliceMutr workflow

This workflow is used by `open-neo run` when paired tumor RNA FASTQ files are
provided. STAR alignment and junction extraction run first. SNAF and
SpliceMutr are independent optional discovery branches; missing optional
evidence is reported as `UNASSESSED`, not negative.

## SNAF

Required assets:

- SNAF Python environment (`SNAF_PYTHON`)
- official SNAF reference root (`NEOAG_SNAF_DB`)
- AltAnalyze container (`NEOAG_ALTANALYZE_IMAGE`)
- patient HLA class-I alleles

The SNAF reference root must contain:

```text
Alt91_db/Hs_Ensembl_exon_add_col.txt
Alt91_db/mRNA-ExonIDs.txt
Alt91_db/Hs_gene-seq-2000_flank.fa
controls/GTEx_junction_counts.h5ad
```

Execution path:

```text
RNA FASTQ -> STAR BAM -> AltAnalyze junction matrix -> SNAF/GTEx filtering
          -> snaf_candidates.tsv -> splice consensus and safety review
```

The AltAnalyze step creates a technical duplicate only to satisfy its legacy
two-group matrix builder. The duplicate column is removed before SNAF and does
not count as biological replication.

Verify the installation:

```bash
export NEOAG_SNAF_DB="$OPEN_NEO_REFERENCE_ROOT/data/snaf/reference/data"
export SNAF_PYTHON="$OPEN_NEO_WORK_ROOT/envs/neoag-snaf/bin/python"
export NEOAG_ALTANALYZE_IMAGE="neoag-altanalyze:snaf"
bash scripts/verify_snaf_splicemutr_assets.sh
```

## SpliceMutr

The release includes a production cohort workflow at
`workflows/splicemutr/SpliceMutr.smk`. It runs STAR-junction conversion,
LeafCutter differential intron usage, official SpliceMutr transcript/ORF
reconstruction and exports `splicemutr_candidates.tsv` for NeoAg.

Copy and edit both templates before execution:

```text
configs/workflows/splicemutr.cohort.example.yaml
configs/workflows/splicemutr.samples.example.tsv
```

The sample sheet requires `sample_id`, `role` (`normal` or `tumor`) and the
matching STAR `SJ.out.tab`. Normal samples are always written first as the
LeafCutter reference group. Production defaults require at least two normal
samples and one tumor sample; `allow_low_power: true` is available for method
development but is marked `LOW_POWER` and must not be called cross-validated.

Run directly:

```bash
splicemutr-neoag workflow workflows/splicemutr/SpliceMutr.smk \
  --configfile configs/workflows/splicemutr.cohort.yaml --cores 16
```

Or through Skill2:

```bash
open-neo run ... \
  --splicemutr-config configs/workflows/splicemutr.cohort.yaml \
  --splicemutr-samples configs/workflows/splicemutr.samples.tsv
```

Do not label a single tumor sample as cross-validated SpliceMutr evidence.
Without a compatible cohort or normal-control junction matrix, the branch is
`UNASSESSED`. RegTools junction support and GTEx normal-junction review still
run and remain available to the downstream evidence model.

Recommended controls are assay- and genome-build-matched normal RNA junctions.
Tissue-specific GTEx files, such as liver, are additional safety backgrounds;
they do not replace the complete SNAF GTEx database or a SpliceMutr cohort.

## Skill2 example

```bash
open-neo run \
  --sample-manifest configs/sample.yaml \
  --snaf-db "$OPEN_NEO_REFERENCE_ROOT/data/snaf/reference/data" \
  --snaf-python "$OPEN_NEO_WORK_ROOT/envs/neoag-snaf/bin/python" \
  --altanalyze-image neoag-altanalyze:snaf \
  --normal-junctions "$OPEN_NEO_REFERENCE_ROOT/data/normal/junctions/normal_junctions.recount3_gtex_v8_grch38.tsv" \
  --mode execute --approved --outdir results/CASE001
```

The junction catalog above is a recount3 re-quantification of GTEx v8 RNA-seq,
not the separately versioned GTEx v11 expression reference. Its filename,
metadata and run manifest must retain `recount3_GTEx_v8_GRCh38`. Catalog
non-membership is reported as `NOT_LISTED_IN_NORMAL_CATALOG`; it is never
treated as an adequate-coverage negative or as proof of tumor specificity.

When no compatible normal RNA cohort is available, structural novelty and
tumor specificity remain separate. An altered junction-spanning sequence may
be retained for discovery, but its tumor-specificity state is
`UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT`, normal safety is at most N1 and
the final evidence grade is capped at R3.

## Formal Splice gate

Formal candidates must pass the ordered biological funnel before HLA
presentation is interpreted as candidate support:

```text
exact junction -> unique reads -> total junction coverage -> PSI ->
coverage-qualified matched normal -> coverage-aware normal cohort ->
annotated normal isoform exclusion -> ORF/frame -> NMD ->
junction-spanning peptide -> normal proteome/transcriptome -> HLA presentation
```

`FORMAL_SPLICE_CANDIDATE` is reserved for events that pass every prerequisite.
Rows with missing fields remain in `EXPLORATION_EVIDENCE_INCOMPLETE`; optional
HLA predictions may be retained for technical review, but are capped at R3 and
are not counted as formal presentation passes in patient reports.

To override the bundled workflow with a reviewed site-specific implementation, add:

```bash
  --splicemutr-workflow configs/workflows/splicemutr.site.smk
```
