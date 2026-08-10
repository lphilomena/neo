# NeoAg Event Pipeline v0.5.0 Formal Splice Provenance Layer

NeoAg Event Pipeline is a research-oriented neoantigen prioritization pipeline. It converts SNV/InDel, fusion, splice, structural-variant, and peptide-only candidates into standardized event and peptide-HLA tables, then layers presentation, APPM, CCF, safety, immune-escape, validation-plan, and report evidence.

This package is a lightweight online release. It includes source code, CLI entry points, Nextflow workflows, tests, fixtures, profiles, setup scripts, and documentation. It does not bundle large references, licensed tools, conda environments, cached work directories, real patient data, or production results.

Important boundary: the pipeline produces computational triage and validation-planning outputs. It does not make clinical diagnoses, clinical resistance calls, or validated treatment recommendations.

## v0.5.0 Formal Splice Provenance Layer

v0.5.0 upgrades the v0.4.4 exact-junction repair into a formal, referentially intact splice provenance model:

```text
junction → splice event → transcript hypothesis → ORF → peptide origin → peptide-HLA presentation
```

The release adds canonical entity registries, event-to-junction and peptide-origin link tables, conservative normal-background state handling, independent evidence-group consensus, strict pVACbind FASTA-index mapping, and compatibility projections back to `raw_events.tsv`, `raw_peptides.tsv`, and `rna_junction_evidence.tsv`. The production entry point is `scripts/run_splice_provenance_v050.sh`; the Python CLI is `neoag-splice-layer`. See `docs/V050_SPLICE_PROVENANCE_LAYER.md` and `CHANGELOG_V050_SPLICE_PROVENANCE_LAYER.md`.

## v0.4.4 Exact Junction and Provenance Repair

v0.4.4 introduces canonical splice-junction identities, removes gene/nearest-locus read transfer, separates caller-provided counts from verified exact-junction support, preserves every source row during production merges, and requires exact canonical agreement for splice cross-domain confirmation. See `CHANGELOG_V044_SPLICE_PROVENANCE.md` and `docs/V044_SPLICE_PROVENANCE.md`.

## What It Does

The pipeline can:

- Parse pVACtools-like SNV/fusion/splice outputs into `raw_events.tsv` and `raw_peptides.tsv`.
- Generate sliding-window variant peptides from VEP-annotated VCFs, with optional automatic VEP annotation when CSQ annotations are missing.
- Score MHC presentation evidence from NetMHCpan, MHCflurry, and optional stability/immunogenicity tools.
- Build APPM 2.0 evidence, including input completeness, conflicts, peptide modifiers, and immune-context annotations.
- Estimate CCF/clonality from purity, CNV, and VAF context.
- Build peptide safety evidence from normal expression, normal ligandome, normal junction, matched-normal, and reference-proteome context.
- Build immune-escape evidence from HLA LOH, APPM, CCF, B2M/JAK/APM context, and related evidence tables.
- Generate long-peptide and minigene validation designs for frameshift, splice, exon-junction, fusion, and SV candidates.
- Produce both patient-facing and technical HTML reports.
- Run fixture workflows through the CLI or the included Nextflow wrappers.

The `.tsv` suffix in ranked outputs is a schema-compatibility label. It is not the software version. The current release is v0.5.0 and writes schema-compatible tables so older downstream scripts can keep reading the same filenames.

## Agent Skills And Coordinator

This release includes a repo-scoped agent skills pack under `.agents/skills/` and a lightweight coordinator CLI:

```bash
neoag-agent --message "compare recommendation and NetMHCpan42 rankings" --result-dir results/sample --outdir work/agent_plan
```

Default mode is dry-run planning. Add `--execute` for supported low-impact skills. See `docs/AGENT_SKILLS_P0_P1.md` for the skill list, expected inputs, outputs, and interpretation boundaries.

For new-machine deployment and production execution, use the three Open-Neo macro Skills below. They replace the older script-by-script migration notes.

## New Machine Install And Run: Three Macro Skills

Use the three public Open-Neo macro Skills for new-machine deployment and case
execution. The machine-readable manifests remain the source of truth:

- `.agents/skills/open-neo-install-check/SKILL.md`: Skill1, machine setup,
  reference/tool discovery, approved install/repair, Doctor, smoke tests and
  production-readiness checks.
- `.agents/skills/open-neo-run/SKILL.md`: Skill2, input QC, route selection,
  Gateway-controlled execution/resume, multi-tool evidence generation and
  weighted plus evidence-consensus ranking.
- `.agents/skills/open-neo-review/SKILL.md`: Skill3, read-only result review,
  experiment-priority tables and patient/technical reports.

### Skill1: install and verify the machine

Clone the release branch on the target machine:

```bash
mkdir -p /home/na/project
git clone --branch na0707_upload_release \
  https://github.com/lphilomena/neo.git \
  /home/na/project/neo
cd /home/na/project/neo
```

If the console script is not yet on `PATH`, use the module entrypoint from the
project environment:

```bash
export PYTHONPATH="$PWD/src"
alias open-neo='/home/na/miniforge3/envs/neoag-tools/bin/python -m neoag.open_neo.cli'
```

Run the install-check macro. For production use, the default installer profile
is `all-open`; it installs the open production tool set where licenses permit,
synchronizes full production reference assets, and writes configured manifests
under `configs/local/`.

```bash
open-neo install-check \
  --project-root "$PWD" \
  --deployment-tier full \
  --mode install \
  --installer-profile all-open \
  --asset-source-host na@10.200.50.134 \
  --asset-source-root /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git \
  --tools-root /home/na/project/open-neo-deploy/env_tool \
  --reference-root /home/na/project/open-neo-deploy/refs \
  --conda-base /home/na/miniforge3 \
  --allow-download \
  --approved \
  --outdir work/install-check-full
```

Use `--installer-profile minimal` for review/core use and
`--installer-profile standard` for the lighter production main path. Re-running
Skill1 is safe: installed tools, synchronized assets and PASS checkpoints are
reused when signatures still match. If installation was interrupted, resume it
instead of starting from scratch:

```bash
open-neo install-check \
  --project-root "$PWD" \
  --deployment-tier full \
  --mode resume \
  --installer-profile all-open \
  --approved \
  --outdir work/install-check-full
```

Full installs default to no wall-clock timeout. If an operator supplies
`--install-timeout SECONDS`, interruption or timeout terminates the whole
installer process group before writing a controlled checkpoint.

### Skill2: run a case through Gateway

Start a local NeoAg Gateway on the target machine. Keep it bound to
`127.0.0.1` unless there is a reviewed reason to expose it.

```bash
cd /home/na/project/neo
source /home/na/miniforge3/bin/activate neoag-tools

mkdir -p work/neoag_gateway
nohup env PYTHONPATH="$PWD/src" \
  python -m neoag.controlled_execution.gateway \
  --host 127.0.0.1 \
  --port 8000 \
  --project-root "$PWD" \
  --outdir "$PWD/work/neoag_gateway" \
  --allowed-root "$PWD/work" \
  --allowed-root /mnt/zzbnew/Public/neoag_results \
  > work/neoag_gateway/gateway.log 2>&1 &

curl -s http://127.0.0.1:8000/health
```

Submit a full run from a sample manifest:

```bash
open-neo run \
  --sample-manifest configs/local/sample.yaml \
  --tools-manifest configs/local/tools_manifest.configured.yaml \
  --reference-manifest configs/local/reference_manifest.configured.yaml \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode execute \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

Direct BAM/VCF/RNA FASTQ inputs are also supported when a manifest is not yet
available:

```bash
open-neo run \
  --sample-id CASE001 \
  --tumor-dna-bam /path/to/tumor.bam \
  --normal-dna-bam /path/to/normal.bam \
  --somatic-vcf /path/to/somatic.pass.vcf.gz \
  --tumor-rna-fastq /path/to/R1.fq.gz /path/to/R2.fq.gz \
  --rna-threads 12 \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode execute \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

Skill2 runs domain tools concurrently where safe. For paired GRCh38 DNA, the
purity/CNV evidence stage can run FACETS, Sequenza and PURPLE and then build a
cross-tool purity/ploidy recommendation. HLA LOH waits for a non-single-tool
purity consensus before launching LOHHLA and SpecHLA in parallel. If both HLA
LOH tools produce usable output the report is labelled `dual_tool_consensus`;
if one fails or has no output, the run continues with `single_tool_result`
evidence explicitly recorded in `hla_loh_tool_status.tsv`,
`hla_loh_summary.json`, `hla_loh_review.md` and `recommended_hla_loh.tsv`.

Resume an interrupted case with:

```bash
open-neo run \
  --result-dir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode resume \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

### Skill3: review and report results

After Skill2 finishes ranking, run Skill3 read-only on the result directory:

```bash
open-neo review \
  --result-dir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --reports patient,technical,onepage \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001/review
```

Skill3 does not rerun heavy tools. It checks result integrity, compares weighted
and evidence-consensus rankings, emits event-level experiment-priority tables,
and writes bounded patient/technical reports. Missing or single-tool evidence
is reported as partial evidence, never as a negative biological result.

## Quick Start

Run these commands from the project root:

```bash
python -m pip install -e .
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
```

Important demo outputs include:

- `work/demo_v043/scoring/ranked_peptides.tsv`
- `work/demo_v043/scoring/ranked_peptides.evidence_consensus.tsv` (parallel review ranking; does not replace the primary ranking)
- `work/demo_v043/scoring/ranked_events.evidence_consensus.tsv`
- `work/demo_v043/scoring/evidence_states.tsv`
- `work/demo_v043/scoring/evidence_conflicts.tsv`
- `work/demo_v043/scoring/weighted_vs_consensus_comparison.tsv`
- `work/demo_v043/scoring/ranked_events.tsv`
- `work/demo_v043/scoring/validation_plan.tsv`
- `work/demo_v043/reports/evidence_report.html`
- `work/demo_v043/reports/evidence_report.patient.html`
- `work/demo_v043/reports/evidence_report.technical.html`
- `work/demo_v043/appm/appm_summary.tsv`
- `work/demo_v043/appm/appm_peptide_modifiers.tsv`
- `work/demo_v043/clonality/ccf_lite.tsv`
- `work/demo_v043/safety/peptide_safety.tsv`
- `work/demo_v043/immune_escape/peptide_escape_flags.tsv`

### Parallel Evidence-consensus Ranking

Every pipeline run keeps the existing weighted `ranked_peptides.tsv`. After
building `comprehensive_peptide_evidence.tsv`, it independently writes peptide
and event consensus rankings, normalized evidence states, and a row-level
weighted-versus-consensus comparison. Missing evidence remains explicit and is
not described as a negative biological result. Candidates receive an R1-R4
evidence grade and a Pareto front calculated within the same event track.

To generate it for an existing ranked table:

```bash
neoag evidence-rank \
  --comprehensive-evidence results/sample/scoring/comprehensive_peptide_evidence.tsv \
  --weighted-baseline results/sample/scoring/ranked_peptides.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v1.toml \
  --provenance results/sample/provenance.json \
  --outdir results/sample/scoring/evidence_consensus \
  --mode parallel --track all \
  --emit-event-ranking --compare-weighted --deterministic
```

Agents should use the thin public Skill2 wrapper instead of reimplementing the
ranking logic:

```bash
neoag-skill run open-neo-run \
  --outdir results/sample/scoring/evidence_consensus \
  --arg comprehensive_evidence=results/sample/scoring/comprehensive_peptide_evidence.tsv \
  --arg weighted_baseline=results/sample/scoring/ranked_peptides.tsv
```

`open-neo-run` and the compatibility skill `neoag-ranking` both invoke the
production `neoag evidence-rank` CLI. The Skill layer contains no independent
R1-R4, Pareto, hard-fail, priority-cap, or event-deduplication implementation.

Compare any two rankings with the generalized audit command:

```bash
neoag-ranking-compare \
  --left results/sample/scoring/ranked_peptides.weighted_baseline.tsv \
  --left-name weighted_baseline \
  --right results/sample/scoring/ranked_peptides.evidence_consensus.tsv \
  --right-name evidence_consensus \
  --outdir results/sample/scoring/ranking_comparison
```

The comparison includes Top10/20/50/100 overlap, Spearman correlation,
promotions/demotions, high-ranking hard-fail candidates, event-type and HLA
composition, evidence conflicts, missing evidence, and manual-review cases.

The first phase does not expose `--replace-primary-ranking`. See
`docs/EVIDENCE_CONSENSUS.md` for field definitions and interpretation.

For tests:

```bash
python -m pip install -e '.[test]'
pytest -q
```

The default test command intentionally skips integration, benchmark, and external-tool tests.

## Common Run Commands

### Prepare The Environment

For fixture-only development:

```bash
python -m pip install -e '.[test]'
pytest -q
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
```

For runs that need external tools:

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag check-tools
```

For a smaller development/test environment:

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
```

### Complete Production Workflow

The production flow runs HLA typing, purity/CNV/CCF, HLA LOH/APPM, and RNA expression in parallel, then merges SNV/InDel, fusion, and splice peptides before presentation prediction and ranking. Fusion and splice events are not considered complete until their junction-derived peptides pass NetMHCpan/MHCflurry and the shared immunogenicity and safety layers.

See `docs/PRODUCTION_WORKFLOW.md`. Use `configs/workflows/production_workflow.example.toml` with `neoag-production-run` for full stage execution, or `conf/run.production_multisource.example.toml` for merge-only runs. Missing fusion/splice sources are reported as `LOW_CONFIDENCE` without stopping the run; missing required presentation predictors still fail clearly.

The following script runs the complete manifest-driven workflow. On first use it creates a private manifest and stops so that sample paths, HLA typing, BAM/FASTQ/VCF inputs, references, and licensed tools can be configured before execution.

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/root/neo/src/na0707_upload_release}"
SAMPLE_ID="${SAMPLE_ID:-SAMPLE001}"
MANIFEST="${MANIFEST:-${PROJECT_ROOT}/configs/workflows/production_workflow.private.toml}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/results/${SAMPLE_ID}_production}"

cd "$PROJECT_ROOT"
source conf/tools.env.sh
python -m pip install --no-deps -e .

if [[ ! -f "$MANIFEST" ]]; then
  cp configs/workflows/production_workflow.example.toml "$MANIFEST"
  echo "Created $MANIFEST"
  echo "Edit every /path/... value and stage command, then rerun this script."
  exit 2
fi

# Environment/tool preflight, followed by a command-only production plan.
neoag check-tools
neoag-production-run \
  --project-root "$PROJECT_ROOT" \
  --manifest "$MANIFEST" \
  --outdir "$OUTDIR"

# Execute missing stages and reuse stages whose declared outputs already exist.
neoag-production-run \
  --project-root "$PROJECT_ROOT" \
  --manifest "$MANIFEST" \
  --outdir "$OUTDIR" \
  --execute

# Minimum acceptance checks for the unified production result.
test -s "$OUTDIR/production_run_summary.json"
test -s "$OUTDIR/production_stage_status.tsv"
test -s "$OUTDIR/final/scoring/ranked_peptides.tsv"
test -s "$OUTDIR/final/scoring/comprehensive_peptide_evidence.tsv"
test -s "$OUTDIR/final/reports/evidence_report.html"

python -m json.tool "$OUTDIR/production_run_summary.json"
```

Do not add `--force` for routine resumes. It is only for intentionally recomputing stages whose declared outputs already exist. A final `LOW_CONFIDENCE` status means ranking completed but one or more optional sources, commonly fusion or splice, were unavailable or failed; inspect `production_stage_status.tsv` before interpreting the ranking.

### Run From Existing pVAC-like Tables

Use this when you already have pVACseq/pVACfuse/pVACsplice-like aggregated tables:

```bash
neoag run \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --pvac data/fixtures/pvacseq_aggregated.tsv \
  --immunogenicity-stub
```

### Run From Pre-built Raw Intermediates

Use this when `parsed/raw_events.tsv` and `parsed/raw_peptides.tsv` already exist:

```bash
neoag run \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --netmhcpan results/sample/presentation/netmhcpan.xls \
  --mhcflurry results/sample/presentation/mhcflurry.csv \
  --expression results/sample/parsed/expression.tsv \
  --hla-loh results/sample/tools/hla_loh.tsv \
  --purity results/sample/tools/purity.tsv \
  --cnv results/sample/tools/cnv_segments.tsv
```

### Sliding-window Variant Peptides To Ranking

Use this path when you have a somatic SNV/InDel VCF and want to generate mutant peptides by sliding window, predict peptide-HLA presentation, and rank event/peptide candidates.

If the VCF already contains VEP `CSQ` annotations, the pipeline uses it directly. If `CSQ` is missing, `run-full` will run VEP annotation first when VEP, cache, reference FASTA, and plugins are configured.

```bash
cat > conf/run.sliding.private.toml <<'TOML'
[sample]
id = "SAMPLE001"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "/path/to/sample.somatic.pass.vcf.gz"
tumor_sample_name = "TUMOR"
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
extract_appm_from_vcf = false
normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
TOML

neoag run-full \
  --config conf/run.sliding.private.toml \
  --outdir results/SAMPLE001_sliding
```

Key outputs:

- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv`
- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.annotated.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_events.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_peptides.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_events.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_peptides.tsv`
- `results/SAMPLE001_sliding/scoring/comprehensive_peptide_evidence.tsv` (all annotation, predictor, RNA, APPM, CCF, safety, escape, ranking, and validation evidence)
- `results/SAMPLE001_sliding/scoring/validation_plan.tsv`
- `results/SAMPLE001_sliding/reports/evidence_report.html`
- `results/SAMPLE001_sliding/reports/evidence_report.patient.html`
- `results/SAMPLE001_sliding/reports/evidence_report.technical.html`

When `inputs.purity_recommendation` points to `purity_recommendation.json`, the full
pipeline uses the multi-tool median for concordant calls. All tool values and the
consensus range/status are retained in CCF outputs; strong discordance forces low confidence.

Manual debug path for variant peptide extraction:

```bash
neoag extract-variant-peptides \
  --input-vcf /path/to/sample.vep.annotated.vcf.gz \
  --output results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv \
  --sample-id SAMPLE001 \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --hla-alleles HLA-A*02:01,HLA-B*07:02,HLA-C*07:02 \
  --tumor-sample-name TUMOR \
  --normal-proteome-fasta /path/to/Homo_sapiens.GRCh38.pep.all.fa \
  --filter-normal-proteome
```

For smoke tests without licensed predictors, set `stub = true` in the TOML or add `--immunogenicity-stub` to direct `run` calls. For production ranking, use real NetMHCpan/MHCflurry outputs and real normal-expression/normal-ligand evidence instead of fixture resources.

### Build Standard Evidence Sidecars

```bash
neoag build-evidence-layer \
  --outdir results/sample \
  --profile default \
  --sample-id SAMPLE001 \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --expression results/sample/parsed/gene_expression.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv \
  --fusion-evidence results/sample/parsed/fusion_evidence.tsv \
  --normal-expression resources/normal_expression.example.tsv \
  --normal-hla-ligands resources/normal_hla_ligands.example.tsv
```

### HLA LOH Conversion And Cross-check

```bash
neoag convert-lohhla \
  -i results/sample/tools/LOHHLA.HLAlossPrediction_CI.xls \
  -o results/sample/tools/lohhla.hla_loh.tsv

neoag convert-spechla \
  -i results/sample/tools/merge.hla.copy.txt \
  -o results/sample/tools/spechla.hla_loh.tsv

neoag crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
```

### Generate Reports

Generate the default combined report plus patient and technical audience-specific reports:

```bash
neoag report \
  --profile default \
  --ranked-events results/sample/scoring/ranked_events.tsv \
  --ranked-peptides results/sample/scoring/ranked_peptides.tsv \
  --appm-summary results/sample/appm/appm_summary.tsv \
  --validation-plan results/sample/scoring/validation_plan.tsv \
  --outdir results/sample \
  --audience both \
  --out results/sample/reports/evidence_report.html
```

### Nextflow Fixture Run

Use the project wrapper rather than calling `nextflow` directly. The wrapper prioritizes the current checkout's `bin/neoag`, sets project paths, and avoids writing Nextflow metadata into a root-owned location.

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

For command-specific options:

```bash
neoag <command> --help
```

## Workflow Dependency Matrix

| Workflow / command | Minimal inputs | Tools | Reference/data |
| --- | --- | --- | --- |
| Fixture demo: `neoag run-demo --outdir work/demo_v043 --sample-id DEMO001` | Bundled fixtures | None beyond Python package | Bundled fixtures/resources |
| Parsed pVAC results: `neoag run --outdir results/sample --sample-id SAMPLE001 --pvac data/fixtures/pvacseq_aggregated.tsv --immunogenicity-stub` | pVAC-like TSVs | None if inputs already exist | Optional normal expression/ligand tables |
| Raw intermediates: `neoag run --outdir results/sample --raw-events raw_events.tsv --raw-peptides raw_peptides.tsv` | `raw_events.tsv`, `raw_peptides.tsv` | NetMHCpan/MHCflurry outputs if provided; optional evidence tools | Optional expression, LOH, purity, CNV, normal evidence |
| Full upstream run: `neoag run-full --config conf/run.sample.private.toml --outdir results/sample` | Run config | Depends on enabled tools | Depends on enabled tools |
| Binding prediction only: `peptide-predict` | Peptide/HLA table | NetMHCpan, MHCflurry, PRIME/BigMHC/DeepImmuno as selected | HLA alleles; predictor model data |
| VEP annotation: `vep-annotate` | VCF | VEP | VEP cache, reference FASTA, plugins |
| Variant peptide extraction: `extract-variant-peptides` | VEP-annotated VCF | Python; optional VEP pre-step | Reference FASTA, optional normal proteome |
| WES SNV calling: `snv-call-wes` | Tumor/normal BAM | GATK4 | GRCh38 FASTA, gnomAD AF VCF, PoN, intervals as needed |
| WES SNV full: `snv-run-full-wes` | Somatic VCF or BAMs | GATK if BAM mode; pVAC/binding tools if enabled | GRCh38 FASTA, HLA, optional normal evidence |
| SV WGS raw build: `sv-build-raw` | SV VCF, FASTA, GTF, HLA | Python | Reference FASTA, GTF, HLA file |
| SV WES raw build: `sv-build-raw-wes` | SV VCF, FASTA, GTF, HLA, capture BED | Python | Reference FASTA, GTF, capture BED, HLA file |
| SV score: `sv-score` | Raw events/peptides | NetMHCpan/MHCflurry unless `--binding-stub` | HLA alleles, optional evidence tables |
| Long-read SV wrapper | FASTQ/BAM or Sniffles2 VCF | minimap2/samtools/Sniffles2 as selected | Reference FASTA, GTF, HLA |
| Fusion discovery | FASTQ/BAM or caller outputs | STAR-Fusion, FusionCatcher, Arriba, EasyFuse as selected | CTAT/EasyFuse/fusion caller references |
| Splice discovery | RNA BAM/junctions + annotated VCF | pVACsplice, RegTools, SNAF, SpliceMutr; optional ASNEO/NeoSplice/splice2neo | GRCh38 FASTA/GTF, genome-specific BSgenome/TxDb and STAR index for SpliceMutr, HLA alleles, normal-junction background when available |
| Immune escape evidence: `immune-escape` | Raw peptides, APPM/CCF/LOH evidence | Optional LOHHLA/FACETS upstream | HLA LOH, CNV, VEP/APM/JAK/B2M evidence |
| Nextflow fixture | Bundled pVAC fixture | Java/Nextflow runtime | Bundled fixtures; writable `NXF_HOME` |

## Tests

Default pytest runs fast unit tests only:

```bash
pytest -q
```

Run broader groups explicitly:

```bash
pytest -q --run-integration
pytest -q --run-benchmark
pytest -q --run-external
pytest -q --run-all
```

Marker form is also supported:

```bash
pytest -q -m unit
pytest -q -m integration --run-integration
pytest -q -m benchmark --run-benchmark
pytest -q -m external --run-external
```

This split prevents lightweight release users from accidentally running long Nextflow, benchmark, or external-tool tests with plain `pytest`.

## Common Errors And Fixes

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `open-neo` is not found | Console script is not on `PATH` yet | Use the module entrypoint shown in the Skill1 section, then rerun `open-neo install-check`. |
| Install, tool, reference, or smoke check fails | Machine is not production-ready yet | Rerun `open-neo install-check --mode resume --approved` and read `deployment_report.md` plus `production_run_readiness.tsv`. |
| Case execution stops or disconnects | Gateway job or production stage was interrupted | Rerun `open-neo run --mode resume --approved --gateway-url http://127.0.0.1:8000 --gateway-wait`. |
| Evidence is marked `UNASSESSED`, `PARTIAL`, or `single_tool_result` | One input/tool/reference was missing or only one tool produced usable evidence | Continue review, but keep the partial-evidence label in reports and experiment decisions. |

## Release Boundary

Do not commit or package:

- `.git`, `.venv`, `.nextflow`, `.pytest_cache`
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`
- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- `conf/*.private.toml`
- real patient data or sample identifiers
- licensed tool binaries
- large references such as `data/ref` and `data/vep`

Use `scripts/check_release_boundary.sh` before preparing an online release.

## Additional Documentation

- `docs/V043_CCF21.md`: CCF 2.1 notes.
- `docs/V042_P1_APPM_EXPLAINABILITY.md`: APPM explainability notes.
- `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`: safety and immune-escape evidence notes.
- `RELEASE.md`: release boundary and test summary.

## Interpretation Boundary

This pipeline is for research triage and validation planning. Ranked candidates should be reviewed with assay validation, disease context, HLA typing, tumor purity, expression/protein support, safety evidence, immune-escape context, and appropriate clinical or research governance.

### NetMHCpan 4.2c container runtime

For servers where the official NetMHCpan 4.2c binary cannot run because of `tcsh` or glibc compatibility, use the Docker/Apptainer runtime documented in [docs/NETMHCPAN_CONTAINER.md](docs/NETMHCPAN_CONTAINER.md). The image contains only OS dependencies; the licensed official `tools/netMHCpan` directory is mounted at runtime.

### Priority tool containers

Docker/Apptainer runtimes for NetMHCpan, NetMHCstabpan, HLA-LA, SpecHLA, PURPLE/AMBER/COBALT, and EasyFuse are documented in [docs/PRIORITY_TOOL_CONTAINERS.md](docs/PRIORITY_TOOL_CONTAINERS.md). These images contain only runtime dependencies; licensed tools and large reference data are mounted from host paths.

## LLM-assisted Coordinator P1

This release adds an optional LLM-assisted Coordinator layer on top of the P0 Skills Pack. The default mode is dependency-free and rule-based; installing the optional `agent-llm` extra enables LiteLLM/LangGraph integration.

Plan only:

```bash
neoag-llm-agent --message "compare recommendation and NetMHCpan42 rankings" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_plan --mode plan
```

Execute safe Skills:

```bash
neoag-llm-agent --message "compare recommendation and NetMHCpan42 rankings" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_execute --mode execute-safe
```

Local Qwen/vLLM through LiteLLM/OpenAI-compatible API:

```bash
neoag-llm-agent --message "update patient report" \
  --file evidence_report.v04x_latest.html \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_report --mode execute-safe \
  --llm-provider litellm --model openai/qwen3-32b \
  --api-base http://localhost:8000/v1 --api-key-env LOCAL_VLLM_API_KEY
```

The Coordinator does not replace Project B CLI/Nextflow. It plans and calls registered Skills; high-impact operations such as HPC submission, installation, deletion, and overwrite require explicit approval.

See `docs/LLM_COORDINATOR_P1.md` and `docs/MODEL_API_AND_AGENT_FRAMEWORK_SELECTION.md`.

- [Tool inventory](docs/TOOL_INVENTORY.md): external tools, Docker images, environment variables, references, and licensing boundaries.

## Skills Taxonomy A/B/C/D

This release includes an upgraded NeoAg Skills taxonomy organized into four categories:

- **A Entry adapter skills**: `neoag-vcf`, `neoag-fusion`, `neoag-splice`, `neoag-sv-wgs`, `neoag-sv-wes`, `neoag-peptide-csv`.
- **B Public evidence analysis skills**: `neoag-hla-typing-loh`, `neoag-presentation`, `neoag-expression`, `neoag-rna-evidence`, `neoag-ccf`, `neoag-appm-escape`, `neoag-safety`, `neoag-ranking`.
- **C Review/report/design skills**: `neoag-ranking-compare`, `neoag-experiment-design`, `neoag-patient-report`, `neoag-technical-report`, `neoag-concept-explainer`.
- **D Governance/execution-control skills**: `neoag-input-qc`, `neoag-doctor`, `neoag-tool-reference-qc`, `neoag-run-demo-and-smoke`, `neoag-pipeline-full`, `neoag-release-qc`, `neoag-gateway-submit`, `neoag-hpc-runner`.

Use:

```bash
neoag-skill list
neoag-skill describe neoag-vcf
neoag-skill validate --root . --outdir work/skill_validate
neoag-skill run neoag-peptide-csv --outdir work/peptides --arg peptide_csv=peptides.tsv
```

Skills are SOP wrappers. They do not make clinical decisions, do not include patient BAM/FASTQ/VCF or large references, and high-risk execution paths remain dry-run or human-approval gated.

## Open-Neo Macro Skills

Three public macro entrypoints compose the existing A/B/C/D Skills without duplicating their algorithms:

- `open-neo install-check`: portable deployment planning, approved installation/repair, Doctor, references, and smoke checks.
- `open-neo run`: manifest/CLI/directory input QC, Gateway-controlled production execution or resume, domain-level tool consensus, comprehensive evidence, weighted baseline, and evidence-consensus ranking.
- `open-neo review`: event-level review, an R1/R2 first-batch set, an R3 evidence-completion queue, and patient/technical reports.

For exactly one paired tumor RNA FASTQ input, `open-neo run` now generates the
standard `rna_fusion_splice_v1` production profile automatically. After Doctor
checks and Gateway approval, it can run FASTQ QC, STAR, Salmon, EasyFuse,
STAR-Fusion, Arriba, RegTools, optional reviewed SNAF/SpliceMutr workflows,
fusion/splice cross-validation, peptide generation and the two rankings.
Missing optional tools or normal-background references remain `UNASSESSED` or
`SAFETY_PARTIAL`; they are never reported as negative evidence.

`plan`, `verify`, and review operations are non-destructive. Installation, repair, execution, and resume require explicit approval. See [docs/OPEN_NEO_MACRO_SKILLS.md](docs/OPEN_NEO_MACRO_SKILLS.md).
## v0.5.1 three evidence chains

The formal splice layer now supports an RNA-driven ImmunoPepper + moPepGen branch, a DNA-causal splice2neo + EasyQuant + pVACsplice branch, and a separately audited normal-background + k4neo branch. See `README_v0.5.1.md` and `docs/V051_THREE_EVIDENCE_CHAINS.md`.

## Candidate source-chain confidence C1–C4 (v0.5.2)

NeoAg now emits an event-specific candidate source-chain confidence tier for SNV, InDel, Fusion and Splice candidates. C1–C4 are separate from R1–R4:

```text
C1: complete chain + independent/cross-modal confirmation
C2: complete strong computational/read chain, no orthogonal confirmation
C3: plausible but incomplete/low-power chain
C4: refuted or invalid event/ORF/peptide chain
```

Use compatibility mode first to preserve the existing EC rank while auditing source-chain completeness:

```bash
neoag source-chain --input all_tool_results.tsv \
  --output source_chain_confidence.tsv \
  --requirements-out source_chain_requirements.long.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v2_1_source_chain.toml
```

See `docs/SOURCE_CHAIN_CONFIDENCE_C1_C4.md` for event-specific rules and the integrated research profile.
