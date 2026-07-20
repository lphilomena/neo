# SV Phase 1 WGS Neoantigen Adapter

This extension adds first-phase support for structural-variant-derived neoantigens to `neoag`.

It is designed for **paired tumor-normal WGS**. For WES samples, see [SV_PHASE1_5_WES.md](SV_PHASE1_5_WES.md).

## What this phase does

Input:

- one or more somatic SV VCFs from Manta, SvABA, GRIDSS2, DELLY, or an already merged SV VCF
- GRCh38/GENCODE-compatible FASTA + GTF
- patient HLA class-I alleles
- optional expression, RNA junction, normal-expression, and normal-HLA-ligand evidence

Output:

- `parsed/raw_events.tsv`
- `parsed/raw_peptides.tsv`
- `sv/sv_events.full.tsv`
- `sv/sv_protein_reconstruction.tsv`
- `sv/sv_mutant_proteins.fa`
- `sv/sv_event_to_peptide.tsv`
- `provenance.sv_phase1.json`

The first two files are the existing v0.3 input schema. The `sv/` sidecar files preserve SV-specific provenance.

## Important limitations

This code is a first-phase computational triage adapter. It does not claim clinical-grade SV reconstruction.

Key limitations:

1. BND/fusion reconstruction is heuristic CDS prefix/suffix reconstruction.
2. Complex rearrangement graph reconstruction is not implemented.
3. RNA junction evidence is optional, but DNA-only candidates should be interpreted conservatively.
4. MHC binding is not performed by `sv-build-raw`; run NetMHCpan/MHCflurry afterward using the existing project B tools.
5. Stub or random predictor output should not be used for real ranking.

## Minimal command

```bash
PYTHONPATH=src python -m neoag.cli sv-build-raw \
  --sample-id P001 \
  --profile sv_wgs_phase1 \
  --sv-vcf manta.somaticSV.vcf.gz svaba.somatic.sv.vcf gridss.vcf.gz \
  --callers Manta SvABA GRIDSS2 \
  --reference-fasta /ref/GRCh38.fa \
  --gencode-gtf /ref/gencode.v44.annotation.gtf \
  --hla HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*40:01,HLA-C*03:04,HLA-C*07:02 \
  --expression gene_expression.tsv \
  --rna-junctions sv_rna_junction_reads.tsv \
  --normal-expression normal_expression.tsv \
  --normal-hla-ligands normal_hla_ligands.tsv \
  --outdir results/P001_sv_phase1
```

## Smoke test command

```bash
PYTHONPATH=src python -m neoag.cli sv-build-raw \
  --sample-id SVMINI \
  --profile sv_wgs_phase1 \
  --sv-vcf data/fixtures_sv/mini_sv.vcf \
  --callers GRIDSS2 \
  --reference-fasta data/fixtures_sv/mini_ref.fa \
  --gencode-gtf data/fixtures_sv/mini.gtf \
  --hla data/fixtures_sv/hla.txt \
  --expression data/fixtures_sv/expression.tsv \
  --rna-junctions data/fixtures_sv/rna_junctions.tsv \
  --normal-expression data/fixtures_sv/normal_expression.tsv \
  --normal-hla-ligands data/fixtures_sv/normal_hla_ligands.tsv \
  --outdir results/SVMINI_sv_phase1
```

Then inspect:

```bash
column -t -s $'\t' results/SVMINI_sv_phase1/sv/sv_events.full.tsv | less -S
column -t -s $'\t' results/SVMINI_sv_phase1/parsed/raw_peptides.tsv | less -S
```

## Downstream scoring

After `sv-build-raw`, run binding + score in one step:

```bash
source conf/tools.env.sh

# End-to-end (adapter + scoring) on fixture SV VCF
neoag sv-run-full \
  --sample-id SVMINI \
  --profile sv_wgs_phase1 \
  --sv-vcf data/fixtures_sv/mini_sv.vcf \
  --callers GRIDSS2 \
  --reference-fasta data/fixtures_sv/mini_ref.fa \
  --gencode-gtf data/fixtures_sv/mini.gtf \
  --hla data/fixtures_sv/hla.txt \
  --expression data/fixtures_sv/expression.tsv \
  --rna-junctions data/fixtures_sv/rna_junctions.tsv \
  --normal-expression data/fixtures_sv/normal_expression.tsv \
  --normal-hla-ligands data/fixtures_sv/normal_hla_ligands.tsv \
  --binding-stub \
  --outdir results/SVMINI_sv_e2e

# Or score an existing sv-build-raw output
neoag sv-score \
  --sample-id SVMINI \
  --profile sv_wgs_phase1 \
  --sv-outdir results/SVMINI_sv_phase1/adapter \
  --binding-stub \
  --outdir results/SVMINI_sv_scored
```

Nextflow (fixture smoke test):

```bash
nextflow run workflows/sv_phase1_fixture.nf -c conf/sv_demo.config
```

WGS workflow with scoring enabled by default:

```bash
nextflow run workflows/sv_phase1_wgs.nf \
  --sample_id P001 \
  --tumor_bam tumor.bam --normal_bam normal.bam \
  --reference_fasta GRCh38.fa --gencode_gtf gencode.gtf \
  --hla 'HLA-A*02:01,HLA-B*07:02' \
  --outdir results/P001_sv_phase1 \
  --binding_stub false
```

Use `--run_scoring false` on `sv_phase1_wgs.nf` for adapter-only mode.

## Added files

- `src/neoag/sv/`
- `profiles/sv_wgs_phase1.toml`
- `conf/run.sv_wgs_phase1.example.toml`
- `workflows/sv_phase1_wgs.nf`
- `modules/sv_manta/`
- `modules/sv_svaba/`
- `modules/sv_gridss/`
- `modules/sv_normalize_merge/`
- `modules/sv_build_raw/`
- `data/fixtures_sv/`
- `tests/test_sv_phase1.py`
- `tests/test_sv_appm_safety.py`

## Quality checks included

The package test suite passes with:

```bash
PYTHONPATH=src pytest -q
```

At packaging time this returned: `32 passed`.
