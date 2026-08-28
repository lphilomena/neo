# WES SNV Phase 1

Tumor-normal **WES** somatic SNV calling with GATK Mutect2 + FilterMutectCalls, then pVACseq and score.

## Pipeline

```mermaid
flowchart LR
  A[tumor/normal BAM] --> B[Mutect2 paired]
  B --> C[FilterMutectCalls]
  C --> D[somatic.vcf.gz]
  D --> E[pVACseq + binding]
  E --> F[score]
```

Phase 1 scope:

- WES capture intervals (`-L intervals.bed`)
- Optional gnomAD / PoN for filtering
- Fixture mode: skip BAM, use `data/fixtures_snv/mini_somatic.vcf` + upstream stub

## CLI

```bash
source conf/tools.env.sh

# Fixture E2E (no BAM, no GATK)
neoag snv-run-full-wes \
  --sample-id SNVMINI \
  --profile default \
  --hla data/fixtures_snv/hla.txt \
  --somatic-vcf data/fixtures_snv/mini_somatic.vcf \
  --tumor-sample-name SNVMINI_TUMOR \
  --normal-sample-name SNVMINI_NORMAL \
  --outdir results/SNVMINI_snv_wes_e2e

# Real calling (requires GATK + indexed reference + BAMs)
neoag snv-call-wes \
  --sample-id P001 \
  --tumor-bam P001.tumor.bam \
  --normal-bam P001.normal.bam \
  --reference-fasta GRCh38.fa \
  --intervals-bed wes_capture.bed \
  --tumor-sample-name P001_TUMOR \
  --normal-sample-name P001_NORMAL \
  --outdir results/P001_calling

neoag snv-run-full-wes \
  --sample-id P001 \
  --hla hla.txt \
  --tumor-bam P001.tumor.bam \
  --normal-bam P001.normal.bam \
  --reference-fasta GRCh38.fa \
  --intervals-bed wes_capture.bed \
  --tumor-sample-name P001_TUMOR \
  --normal-sample-name P001_NORMAL \
  --no-upstream-stub \
  --outdir results/P001_snv_wes
```

## Nextflow

Fixture smoke (stub upstream):

```bash
nextflow run workflows/snv_phase1_wes_fixture.nf -c conf/snv_wes_demo.config
```

Full WES (requires BAM + GATK on PATH):

```bash
nextflow run workflows/snv_phase1_wes.nf \
  --sample_id P001 \
  --tumor_bam tumor.bam --normal_bam normal.bam \
  --reference_fasta GRCh38.fa \
  --reference_dict GRCh38.dict \
  --intervals_bed capture.bed \
  --tumor_sample_name TUMOR --normal_sample_name NORMAL \
  -c conf/snv_wes_demo.config
```

## Fixtures

| File | Purpose |
|------|---------|
| `data/fixtures_snv/mini_ref.fa` | Mini reference (+ `.fai`, `.dict`) |
| `data/fixtures_snv/wes_capture.bed` | Capture intervals |
| `data/fixtures_snv/mini_somatic.vcf` | Stand-in Mutect2 filtered VCF |
| `conf/run.snv_wes_fixture.toml` | Static upstream stub config |

## External resources (production)

| Resource | Required for real runs |
|----------|------------------------|
| GRCh38.fa + .fai + .dict | Yes |
| WES capture BED | Yes |
| gnomAD VCF | Recommended for FilterMutectCalls |
| Panel of Normals | Recommended for exome |

Install GATK4 and ensure `gatk` is on PATH after `source conf/tools.env.sh`.

## Limitations

- Phase 1 does not build PoN or run BQSR/indels realignment.
- pVACseq expects annotated/combined VCF for production; stub mode bypasses real pVAC.
- Not clinical-grade variant calling; validate calls independently.

See also: [SV_PHASE1_5_WES.md](SV_PHASE1_5_WES.md) for SV on WES.
