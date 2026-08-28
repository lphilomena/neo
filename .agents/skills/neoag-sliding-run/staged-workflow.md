# Staged NeoAg Sliding-window Workflow

Use this reference when the user wants to split the SNV/InDel sliding-window workflow into VEP annotation, short peptide generation, and downstream scoring/ranking/reporting.

Replace placeholders before running:

- `<sample_id>`
- `<profile>`
- `<variants_vcf>`
- `<annotated_vcf>`
- `<tumor_sample_name>`
- `<hla_csv>`
- `<outdir>`
- `<normal_proteome_fasta>`

## 0. Environment Readiness

Before staged execution, verify the project package and tool environment:

```bash
source conf/tools.env.sh
python -m pip install -e .
neoag check-tools
```

If `neoag` or `neoag.tools` cannot be imported, rerun the editable install from the repository root. If VEP or NetMHCpan is missing and the run needs real local tools, follow the bootstrap instructions in `SKILL.md` before continuing.

## 1. VEP Annotation

Skip this stage only if the input VCF already contains VEP `CSQ` annotations.

```bash
source conf/tools.env.sh

neoag vep-annotate \
  --input-vcf <variants_vcf> \
  --output-vcf <outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz \
  --sample-id <sample_id> \
  --fasta "$NEOAG_REFERENCE_FASTA" \
  --cache-dir "$NEOAG_VEP_CACHE" \
  --plugins-dir "$NEOAG_VEP_PLUGINS" \
  --fork 4
```

Expected output:

- `<outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz`

Quick checks:

```bash
test -s <outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz
test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm"
test -f "$NEOAG_VEP_PLUGINS/Frameshift.pm"
```

## 2. Short Peptide Generation

Use the VEP-annotated VCF from stage 1, or the original VCF if it already has `CSQ`.

```bash
source conf/tools.env.sh

neoag extract-variant-peptides \
  --input-vcf <annotated_vcf> \
  --output <outdir>/upstream/tools/variant_peptides.tsv \
  --sample-id <sample_id> \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --hla-alleles <hla_csv> \
  --tumor-sample-name <tumor_sample_name> \
  --normal-proteome-fasta <normal_proteome_fasta> \
  --filter-normal-proteome
```

If no normal proteome FASTA is available, omit `--normal-proteome-fasta` and `--filter-normal-proteome`, and report that normal-proteome safety filtering was not applied.

Expected outputs from the initial peptide extraction step:

- `<outdir>/upstream/tools/variant_peptides.tsv`
- `<outdir>/upstream/tools/variant_peptides.annotated.tsv` (initial catalog sidecar; refresh after scoring in stage 3)
- `<outdir>/upstream/parsed/raw_events.tsv`
- `<outdir>/upstream/parsed/raw_peptides.tsv`

Quick checks:

```bash
test -s <outdir>/upstream/tools/variant_peptides.tsv
test -s <outdir>/upstream/parsed/raw_events.tsv
test -s <outdir>/upstream/parsed/raw_peptides.tsv
```

## 3. Tool Scoring, Ranking, Validation, Reports

Use this route when peptide-HLA binding still needs to be generated:

```bash
source conf/tools.env.sh

neoag peptide-predict \
  --raw-peptides <outdir>/upstream/parsed/raw_peptides.tsv \
  --outdir <outdir>/presentation \
  --sample-id <sample_id> \
  --tools netmhcpan,mhcflurry
```

Then run ranking and reports from intermediates:

```bash
neoag run \
  --outdir <outdir> \
  --sample-id <sample_id> \
  --profile <profile> \
  --raw-events <outdir>/upstream/parsed/raw_events.tsv \
  --raw-peptides <outdir>/upstream/parsed/raw_peptides.tsv \
  --netmhcpan <outdir>/presentation/netmhcpan.xls \
  --mhcflurry <outdir>/presentation/mhcflurry.csv \
  --normal-expression resources/normal_expression.example.tsv \
  --normal-hla-ligands resources/normal_hla_ligands.example.tsv
```

If binding predictor outputs already exist, skip `peptide-predict` and provide the existing `--netmhcpan` and `--mhcflurry` paths to `run`.

After scoring, refresh the annotated peptide catalog with NetMHCpan/MHCflurry and optional immunogenicity evidence. This makes `variant_peptides.annotated.tsv` the final review table rather than only the extraction sidecar.

```bash
PYTHONPATH=src python .agents/skills/neoag-sliding-run/scripts/refresh_variant_peptides_annotated.py \
  --variant-peptides <outdir>/upstream/tools/variant_peptides.tsv \
  --output <outdir>/upstream/tools/variant_peptides.annotated.tsv \
  --hla-alleles <hla_csv> \
  --netmhcpan-xls <outdir>/presentation/netmhcpan.xls \
  --mhcflurry-csv <outdir>/presentation/mhcflurry.csv \
  --netmhcstabpan-tsv <outdir>/presentation/netmhcstabpan.tsv \
  --prime-tsv <outdir>/tools/prime.tsv \
  --bigmhc-im-tsv <outdir>/tools/bigmhc_im.tsv \
  --iedb-immunogenicity-tsv <outdir>/presentation/iedb_immunogenicity.tsv
```

The refresh script safely ignores optional evidence files that do not exist. Required inputs are `variant_peptides.tsv`, HLA alleles, and at least one predictor output if the user expects tool score columns to be populated.

Expected outputs:

- `<outdir>/scoring/ranked_events.tsv`
- `<outdir>/scoring/ranked_peptides.tsv`
- `<outdir>/scoring/validation_plan.tsv`
- `<outdir>/upstream/tools/variant_peptides.annotated.tsv` (refreshed with tool scores)
- `<outdir>/reports/evidence_report.html`
- `<outdir>/reports/evidence_report.patient.html`
- `<outdir>/reports/evidence_report.technical.html`

Quick checks:

```bash
test -s <outdir>/scoring/ranked_events.tsv
test -s <outdir>/scoring/ranked_peptides.tsv
test -s <outdir>/scoring/validation_plan.tsv
test -s <outdir>/upstream/tools/variant_peptides.annotated.tsv
test -s <outdir>/reports/evidence_report.technical.html
```

## Failure Handling

- VEP/cache/plugin errors: rerun only stage 1 after fixing `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_PLUGINS`, or `NEOAG_REFERENCE_FASTA`.
- Peptide extraction errors: rerun only stage 2 after checking CSQ annotations, tumor sample name, HLA alleles, and normal proteome path.
- NetMHCpan/MHCflurry errors: rerun stage 3 only; do not repeat VEP annotation or peptide extraction unless their inputs changed.
- For smoke tests without licensed tools, use the one-command `run-full` path with `[tools].stub = true`, or use existing fixture predictor outputs with `run`.
