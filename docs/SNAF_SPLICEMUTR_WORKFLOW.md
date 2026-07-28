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

SpliceMutr is cohort-oriented. A site-reviewed Snakemake workflow must be
provided through `--splicemutr-workflow`; it must consume the BAM and junction
paths exported by the wrapper and write `splicemutr_candidates.tsv` to
`NEOAG_SPLICEMUTR_OUTDIR`.

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
  --normal-junctions "$OPEN_NEO_REFERENCE_ROOT/data/normal/junctions/normal_junctions.gtex_v11_combined.tsv" \
  --mode execute --approved --outdir results/CASE001
```

For a licensed or locally reviewed SpliceMutr deployment, add:

```bash
  --splicemutr-workflow configs/workflows/splicemutr.site.smk
```
