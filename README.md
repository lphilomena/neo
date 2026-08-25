# NeoAg Event Pipeline v0.5.0

NeoAg Event Pipeline is a research-oriented neoantigen prioritization pipeline. It converts SNV/InDel, fusion, splice, structural-variant, and peptide-only candidates into standardized event and peptide-HLA tables, then layers presentation, APPM, CCF, safety, immune-escape, validation-plan, and report evidence.

This package is a lightweight online release. It includes source code, CLI entry points, Nextflow workflows, tests, fixtures, profiles, setup scripts, and documentation. It does not bundle large references, licensed tools, conda environments, cached work directories, real patient data, or production results.

Important boundary: the pipeline produces computational triage and validation-planning outputs. It does not make clinical diagnoses, clinical resistance calls, or validated treatment recommendations.

## Project Function Overview

NeoAg Event Pipeline is an end-to-end research workflow for converting tumor DNA/RNA evidence into traceable neoantigen candidates, cross-tool consensus results, prioritized peptide-HLA records, and experiment-oriented reports. Its main functional layers are:

- **Input harmonization:** accepts paired tumor/normal DNA, tumor RNA, somatic VCF, and supported upstream tool outputs, then records input QC and routing decisions.
- **Cross-tool biological evidence:** combines HLA typing, purity/CNV, HLA LOH, expression, RNA allele support, fusion, splice, and structural-variant evidence while retaining disagreements and missing results.
- **Candidate reconstruction:** normalizes SNV/InDel, fusion, splice, and SV events into event, peptide, transcript/ORF, junction, and peptide-HLA provenance tables.
- **Presentation and safety assessment:** integrates MHC presentation, processing, immunogenicity, normal-proteome, normal-expression, ligandome, and normal-junction background evidence.
- **Prioritization and reporting:** produces weighted and evidence-consensus rankings, validation designs, audit artifacts, and patient-facing and technical reports.

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

## System Requirements

Requirements depend on assay size and enabled tools. The following baseline is intended for full production DNA/RNA analysis; review-only use can run with substantially fewer resources.

- **Operating system:** 64-bit Linux; Ubuntu 22.04 LTS is recommended. Ubuntu 24.04 is supported when legacy or licensed binaries are run through the configured compatible environment or container.
- **CPU:** 16 cores minimum; 32-64 cores recommended for WGS, RNA alignment, fusion calling, and parallel prediction.
- **Memory:** 64 GB minimum; 128 GB recommended. Allocate 256 GB or more for large WGS cases or memory-intensive fusion/HLA workflows.
- **Storage:** reserve about 100 GB for code and environments, 0.5-1.5 TB for references and container assets, and approximately 0.5-2 TB of working/output space per WGS-scale case. Fast local SSD scratch is recommended; shared storage is suitable for fixed assets and archived results.
- **Base software:** Bash, Git, SSH, rsync, curl or wget, tar, gzip, and standard build tools. Skill1 provisions the supported Python/Conda environments and checks optional Docker/Apptainer requirements.
- **Network and permissions:** outbound access to the code repository and public asset repository is required when downloading; licensed tools require operator-supplied installers/licenses. The execution account must be able to write to the project, scratch, Gateway, and approved output roots.
- **Input conventions:** use one consistent genome build, normally GRCh38. BAM/CRAM inputs must be indexed; VCF inputs should be normalized, compressed, and indexed where applicable; paired FASTQ files must have consistent sample and read-pair naming.
- **Gateway:** production execution requires NeoAg Gateway. Bind it to `127.0.0.1` by default and authorize only the intended result directories with `allowed-root` settings.

## New Machine Install And Run: Three Macro Skills

Use the three public Open-Neo macro Skills for new-machine deployment and case execution. The machine-readable manifests remain the source of truth. Production-required tools, cross-validation gates, and accepted upstream-result inputs are summarized in `docs/PRODUCTION_REQUIRED_TOOLS.md`. Production ranking requires both NetMHCstabpan stability and NetChop 3.1d cleavage evidence; missing results block production release rather than being silently skipped:

- `.agents/skills/open-neo-install-check/SKILL.md`: Skill1, machine setup, reference/tool discovery, approved install/repair, Doctor, smoke tests and production-readiness checks.
- `.agents/skills/open-neo-run/SKILL.md`: Skill2, input QC, route selection, Gateway-controlled execution/resume, multi-tool evidence generation and weighted plus evidence-consensus ranking.
- `.agents/skills/open-neo-review/SKILL.md`: Skill3, read-only result review, experiment-priority tables and patient/technical reports.

### Skill1: Install And Verify The Machine

Clone the release branch on the target machine:

```bash
mkdir -p /home/na/project
git clone --branch na0707_upload_release \
  https://github.com/lphilomena/neo.git \
  /home/na/project/neo
cd /home/na/project/neo
```

If the console script is not yet on `PATH`, use the module entrypoint from the project environment:

```bash
export PYTHONPATH="$PWD/src"
alias open-neo='/home/na/miniforge3/envs/neoag-tools/bin/python -m neoag.open_neo.cli'
```

Run the install-check macro. For production use, the default installer profile is `all-open`; it installs the open production tool set where licenses permit, synchronizes full production reference assets, and writes configured manifests under `configs/local/`.

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

Use `--installer-profile minimal` for review/core use and `--installer-profile standard` for the lighter production main path. Re-running Skill1 is safe: installed tools, synchronized assets and PASS checkpoints are reused when signatures still match. If installation was interrupted, resume it instead of starting from scratch:

```bash
open-neo install-check \
  --project-root "$PWD" \
  --deployment-tier full \
  --mode resume \
  --installer-profile all-open \
  --approved \
  --outdir work/install-check-full
```

Full installs default to no wall-clock timeout. If an operator supplies `--install-timeout SECONDS`, interruption or timeout terminates the whole installer process group before writing a controlled checkpoint.

### Skill2: Run A Case Through Gateway

Start a local NeoAg Gateway on the target machine. Keep it bound to `127.0.0.1` unless there is a reviewed reason to expose it.

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

Direct BAM/VCF/RNA FASTQ inputs are also supported when a manifest is not yet available:

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

Skill2 runs domain tools concurrently where safe. For paired GRCh38 DNA, the purity/CNV evidence stage can run FACETS, Sequenza and PURPLE and then build a cross-tool purity/ploidy recommendation. HLA LOH waits for a non-single-tool purity consensus before launching LOHHLA and SpecHLA in parallel. If both HLA LOH tools produce usable output the report is labelled `dual_tool_consensus`; if one fails or has no output, the run continues with `single_tool_result` evidence explicitly recorded in `hla_loh_tool_status.tsv`, `hla_loh_summary.json`, `hla_loh_review.md` and `recommended_hla_loh.tsv`.

Resume an interrupted case with:

```bash
open-neo run \
  --result-dir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode resume \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

### Skill3: Review And Report Results

After Skill2 finishes ranking, run Skill3 read-only on the result directory:

```bash
open-neo review \
  --result-dir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --reports patient,technical,onepage \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001/review
```

Skill3 does not rerun heavy tools. It checks result integrity, compares weighted and evidence-consensus rankings, emits event-level experiment-priority tables, and writes bounded patient/technical reports. Missing or single-tool evidence is reported as partial evidence, never as a negative biological result.

## Output Files By Function

The exact directory set depends on available inputs and enabled tools. The generated run/output manifests and provenance files are the authoritative inventory. Common production outputs are grouped below.

### Input QC, Planning, And Audit

- `input_qc/`: input inventory, missing/ambiguous input checks, and the selected route plan.
- `manifests/`: effective sample, tool, reference, capability, and production manifests used for the run.
- `production_stage_status.tsv` and `provenance.json`: stage status, reused/generated evidence, commands, source files, versions, and audit links.

### HLA, Purity/CNV, And Immune Escape

- `hla/consensus/`: normalized HLA typing results and cross-tool HLA consensus.
- `purity/consensus/`: per-tool purity/ploidy summaries, recommended purity, CNV segments, and consensus confidence.
- `hla_loh/`: LOHHLA and SpecHLA allele-level outputs plus HLA-LOH consensus; failed or single-tool evidence remains explicit.
- `appm/` and `immune_escape/`: antigen-processing/presentation machinery evidence and peptide-level immune-escape flags.

### Events, RNA Evidence, And Candidate Peptides

- `branches/snv/`, `branches/fusion/`, `branches/splice/`, and `branches/sv/`: source-specific raw events, raw peptides, normalized caller evidence, and branch consensus outputs.
- `rna/expression/`: gene- and transcript-level expression tables from supported quantification tools.
- `rna/rna_alt_vaf.tsv` and junction evidence: RNA depth, ALT reads, RNA VAF, exact-junction support, and source-to-event linkage where the required RNA evidence is available.
- `merged/raw_events.tsv` and `merged/raw_peptides.tsv`: the unified event and peptide inputs used by downstream evidence and ranking stages.

### Presentation, Safety, Ranking, And Reports

- `presentation/` or predictor-specific evidence tables: NetMHCpan, MHCflurry, NetMHCstabpan, NetChop, and configured supplementary predictor results, followed by presentation consensus.
- `safety/`: normal-proteome, normal-expression, ligandome, matched-normal, and normal-junction background assessments.
- `final/scoring/ranked_peptides.evidence_consensus.tsv`: the preferred final peptide review table with normalized evidence states, conflicts, confidence caps, and R-level prioritization.
- `final/scoring/ranked_events.evidence_consensus.tsv`: event-level deduplicated prioritization linking each event to its representative and alternative peptide-HLA candidates.
- `final/scoring/validation_plan.tsv`: short-peptide, long-peptide, minigene, and orthogonal-validation planning records.
- `reports/`: patient report, technical report, one-page summary, and supporting review tables generated by Skill3.
