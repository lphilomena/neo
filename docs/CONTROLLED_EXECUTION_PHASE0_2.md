# controlled-execution Phase 0-2 enhancement for Project B

This enhancement adds an manifest-driven controlled execution layer on top of the existing
Project B CLI, Skills Pack, and LLM-assisted Coordinator.

## Phase 0: release cleanup + NeoAg Doctor

Phase 0 adds read-only diagnostics and release-boundary auditing.

Commands:

```bash
neoag-doctor \
  --project-root . \
  --tools-manifest configs/controlled_execution/tools_manifest.example.yaml \
  --reference-manifest configs/controlled_execution/reference_manifest.example.yaml \
  --sample-manifest configs/controlled_execution/sample_manifest.example.yaml \
  --outdir work/doctor \
  --dry-run \
  --mini-smoke

neoag doctor \
  --project-root . \
  --tools-manifest configs/controlled_execution/tools_manifest.example.yaml \
  --reference-manifest configs/controlled_execution/reference_manifest.example.yaml \
  --sample-manifest configs/controlled_execution/sample_manifest.example.yaml \
  --outdir work/doctor \
  --dry-run \
  --mini-smoke

neoag-release-audit --root . --outdir work/release_audit
neoag release-audit --root . --outdir work/release_audit
```

Outputs:

- `doctor_status.json`
- `doctor_checks.tsv`
- `tool_status.tsv`
- `reference_status.tsv`
- `smoke_tests.tsv`
- `blocking_issues.tsv`
- `recommended_fixes.md`
- `doctor_summary.md`
- `release_audit/*`

Doctor is intentionally read-only. It checks package importability, core CLI
entrypoints, tool entrypoints, reference/sample manifest paths, release boundary,
and optional smoke tests. It does not install tools, download references, or
modify pipeline configuration.

Additional Phase 0 checks include:

- console-script exposure for `neoag-doctor`, `neoag-release-audit`,
  `neoag-pipeline-full`, and `neoag-gateway` with a clear `pip install -e .`
  hint when scripts are not on `PATH`;
- workflow readiness checks for Java/Nextflow/GATK Mutect2 help, `bin/`
  permissions, writable `work/` and `results/`, and local/slurm profile hints;
- tool-specific manifest/license notes for VEP, NetMHCpan, MHCflurry, LOHHLA,
  FACETS, PURPLE, ASCAT, Sequenza, Arriba, EasyFuse, SpecHLA, HLA-LA, and
  PyClone-VI;
- reference-specific checks for FASTA `.fai`/`.dict`, GENCODE GTF, VEP cache
  layout, FACETS SNP VCF `.tbi`, LOHHLA, ASCAT loci/alleles, Arriba and HLA
  references when declared;
- `--mini-smoke`, which performs read-only imports, CLI `--help` checks,
  NetMHCpan one-peptide/one-HLA smoke when available, MHCflurry CLI/model
  readiness, VEP tiny offline-cache smoke when references are declared, and
  config/reference readiness checks for FACETS, LOHHLA, Arriba and STAR-Fusion;
- release-boundary scanning that records but skips `.git`, `work`, `results`,
  `.nextflow`, `__pycache__`, `.pytest_cache` and similar generated directories
  by default. Use `--scan-generated-dirs` only when a deep scan of runtime
  outputs is required.

## Phase 1: NeoAg Gateway

Phase 1 adds a stdlib HTTP JSON gateway. It is intentionally dependency-light
and safe-by-default. Requests are validated against a route registry, assigned a
`job_id`, written to `jobs/{job_id}.json`, logged to `audit_log.jsonl`, and run
in a background thread. Callers should use `GET /job-status/{job_id}` to poll
completion.

Start:

```bash
neoag-gateway \
  --host 127.0.0.1 \
  --port 8000 \
  --project-root . \
  --outdir work/gateway

# Optional: allow additional write roots outside project/outdir.
neoag-gateway --allowed-root /path/to/allowed/output ...
```

Example requests:

```bash
curl -s http://127.0.0.1:8000/health

curl -s -X POST http://127.0.0.1:8000/doctor \
  -H 'Content-Type: application/json' \
  -d '{"outdir":"work/gateway_doctor","project_root":".","execute":false,"mini_smoke":true}'

curl -s -X POST http://127.0.0.1:8000/ranking-compare \
  -H 'Content-Type: application/json' \
  -d '{"recommendation":"ranked_peptides.recommendation.tsv","netmhcpan42":"ranked_peptides.netmhcpan42.tsv","outdir":"work/gateway_rank"}'

curl -s http://127.0.0.1:8000/job-status/job-xxxxxxxxxx
```

Routes:

- `GET /health`
- `POST /doctor`
- `POST /input-qc`
- `POST /ranking-compare`
- `POST /appm-review`
- `POST /ccf-review`
- `POST /patient-report`
- `POST /pipeline-full`
- `GET /job-status/{job_id}`

Risk policy:

- LOW: read-only/result-review tasks including `doctor` dry-run, `input-qc`,
  `ranking-compare`, `appm-review`, `ccf-review`, `patient-report` draft
  generation, and dry-run `pipeline-full`;
- MEDIUM: reserved for future document/PPT generation, result-directory writes
  beyond the gateway output area, and demo execution;
- HIGH: `pipeline-full execute=true`, `hpc_submit=true`, `overwrite=true`,
  `allow_overwrite=true`, install/delete intent flags. HIGH requests require
  `approved=true`; install/delete flags are rejected because this gateway does
  not implement those actions.

MVP audit fields recorded in job JSON and audit events include:

- `request_id`, `case_id`, `sample_id`, `user`, `timestamp`;
- `input_manifest`, `skill_or_pipeline`, `command_preview`;
- `risk_level`, `approval_status`, `job_id`, `output_dir`;
- `status`, `failure_reason`.

Safety controls:

- missing required fields and unknown fields return `BAD_REQUEST`;
- HIGH requests without approval return `APPROVAL_REQUIRED` and still receive a
  `job_id` for auditability;
- write/output paths must be relative to the project or under configured allowed
  roots;
- shell execution is not exposed through the JSON API; routes call registered
  Python entrypoints or controlled pipeline functions.

The gateway is designed as a local controlled execution layer, not as an
internet-facing production API.


## Manifest-first configuration

For new-machine deployment with a programming agent, use the skill-first
migration wrapper before editing local paths manually:

```bash
bash scripts/bootstrap_agent_deploy.sh
```

This creates local manifest copies, validates skills, runs read-only Doctor and
generates a `pipeline-full` dry-run plan. See `docs/SKILL_FIRST_MIGRATION.md`.

Controlled execution uses manifests as the stable boundary between source code,
external bioinformatics tools and large reference data. Large files such as VEP
cache, reference genome FASTA, annotation databases, HLA graphs and licensed
tool data are not stored in the repository. They are declared in manifests and
validated by Doctor/pipeline-full.

Manifest files:

- `configs/controlled_execution/sample_manifest.example.yaml`
- `configs/controlled_execution/reference_manifest.example.yaml`
- `configs/controlled_execution/tools_manifest.example.yaml`
- `configs/controlled_execution/gateway.example.json`

`sample_manifest` declares `sample_id`, disease/case metadata, tumor/normal BAM
inputs and optional analysis evidence such as somatic VCF, SV VCFs, fusion TSV,
expression TPM, HLA typing/LOH and purity/CNV tables.

`reference_manifest` declares GRCh38 FASTA, GENCODE GTF, VEP cache, protein and
normal safety references, HLA reference data, FACETS SNP VCF and optional tool
reference bundles.

`tools_manifest` accepts either a top-level tool mapping or a `tools:` block.
Each tool can declare `mode`, `executable`/`path`, `image` and
`license_required`. Container image digests written as `image@sha256:...` are
carried into run provenance.

For each `pipeline-full` run, the runner writes:

- `manifest_validation.tsv`
- `input_file_hashes.tsv`
- `reference_hashes.tsv`
- `tool_versions.tsv`
- `container_digests.tsv`
- `run_manifest.json` with `run_id`, `sample_id`, `pipeline_version`, `git_sha`,
  input/reference hashes, tool versions, container digests, start/end time,
  status and `output_manifest` path.

## Phase 2: pipeline-full runner

Phase 2 adds a manifest-driven pipeline-full planning runner.

```bash
neoag-pipeline-full \
  --sample-manifest configs/controlled_execution/sample_manifest.example.yaml \
  --tools-manifest configs/controlled_execution/tools_manifest.example.yaml \
  --reference-manifest configs/controlled_execution/reference_manifest.example.yaml \
  --outdir work/pipeline_full \
  --profile local
```

By default, this is dry-run planning. Use `--execute` to allow safe local
execution. Heavy external bioinformatics steps remain planned unless explicitly
implemented through a Gateway/HPC execution backend.

Outputs:

- `pipeline_plan.md`
- `pipeline_status.tsv`
- `run_manifest.json`
- `output_manifest.json`
- `audit_log.jsonl`
- `doctor/*`

## Design boundary

- Agent plans and explains.
- Skills provide stable SOPs and scripts.
- Gateway controls execution, approvals, and audit.
- Project B CLI/Nextflow/external tools perform computation.
- Missing evidence is reported as missing/unassessed, not as negative evidence.
- Candidate neoantigens are computational triage candidates, not confirmed
  therapeutic targets.
