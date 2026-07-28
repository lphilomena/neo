# Open-Neo public macro Skills

Open-Neo exposes three public task entrypoints while retaining the existing A/B/C/D fine-grained Skills as the internal implementation layer.

## 1. Installation and environment check

```bash
open-neo install-check \
  --project-root . \
  --deployment-tier prediction \
  --mode verify \
  --outdir work/install-check
```

`repair` and `install` delegate to the portable new-machine installer under
`neoag-remote-deploy`; they require `--approved`. Downloads require the
additional `--allow-download` flag, and licensed assets remain host-mounted or
machine-local. A supplied release tarball must include `--sha256`.

Release archives are safely staged under the Skill output directory; absolute,
traversal, link and device members are rejected. The staged tree must contain
one unambiguous `pyproject.toml` or `setup.py` root. Mutating modes record a
timeout-protected `deployment_checkpoint.json`; `--mode resume` reuses only a
matching successful checkpoint and otherwise reruns the idempotent installer.

For prediction/full tiers, READY requires the declared reference manifest,
FASTA sidecars and the required tier assets. Skill1 runs Doctor before and
after an approved installation and records `deployment_delta.tsv` plus the
authoritative `tier_requirements.tsv`.

Skill1 also writes `manifests/production_assets.local.tsv`, replacing generic
tool, licensed-tool and reference targets with the selected local roots. Use
`--asset-source-root` for a mounted asset tree or `--asset-source-host` for a
remote source. Required sources are checked before installation begins.

## 2. Input detection, Pipeline, and two rankings

Plan:

```bash
open-neo run \
  --sample-manifest configs/open_neo/sample_manifest.example.yaml \
  --mode plan \
  --outdir work/case-plan
```

Execute after review and approval:

```bash
open-neo run \
  --sample-manifest sample.yaml \
  --mode execute \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --project-root . \
  --outdir results/CASE001
```

Paired tumor RNA FASTQ can use the built-in fusion/splice profile:

```bash
open-neo run \
  --sample-manifest configs/local/RNA001.sample.yaml \
  --mode plan \
  --outdir work/RNA001-plan
```

With no explicit `production_manifest`, Skill2 generates
`manifests/capability_aware.production.toml`, `capability_decisions.tsv` and
`capability_plan.json`. The planner supports combined DNA/RNA inputs and selects
every compatible tool that has either a repository-owned sample runner or an
administrator-reviewed `command_template`. PATH presence alone is insufficient.
Execution requires Gateway approval; absent optional tools and references remain
`UNASSESSED/PARTIAL` rather than negative.

Use `--automatic-tool-policy all-available` (default) to select every validated
compatible runner. `balanced` and `minimal` are reserved for site policies that
prefer fewer cross-validation tools.

Reuse existing evidence:

```bash
open-neo run \
  --mode ranking-only \
  --comprehensive-evidence results/scoring/comprehensive_peptide_evidence.tsv \
  --weighted-baseline results/scoring/ranked_peptides.tsv \
  --outdir work/ranking-only
```

## 3. Event-level review, experiments, and reports

```bash
open-neo review \
  --result-dir results/CASE001 \
  --clinical-context configs/cases/CASE001.clinical.yaml \
  --disease-profile configs/diseases/sarcoma.yaml \
  --top-n 12 \
  --reports patient,technical \
  --outdir reviews/CASE001
```

Use `--reports none` for event review and experiment tables without generating
patient/technical documents or the one-page PPTX.

`result_dir` is the only required source entrypoint. Before reviewing any
candidate, Skill3 verifies that the run manifest, event- and peptide-level
evidence-consensus rankings, weighted baseline, canonical all-tool evidence,
and validation plan belong to a consistent run. Missing event-level consensus
returns `NEEDS_RANKING`; hash/run mismatches, promoted hard failures, and
incorrect missing-evidence semantics return `BLOCKED`.

Review is event-first. It preserves `pipeline_r_grade` and
`pipeline_event_rank`, then adds independent review fields such as
`review_status`, `review_reason`, `experiment_priority`, and
`recommended_validation`. One event contributes at most one or two
representative peptide-HLA pairs; phase and redundancy groups are deduplicated.

R1/R2 events with adequate RNA, presentation, safety, HLA/APPM, and clonality
evidence may enter the direct experimental-priority set. Eligible R3 events
may enter a clearly labelled evidence-completion set for targeted RNA,
fusion confirmation, or phasing. R4 drivers remain manual-review-only and are
never upgraded because of biological importance alone.

The deterministic first batch is an experimental planning set, not an
optimized vaccine design. It balances event type, HLA coverage, clonality,
RNA evidence, and safety while removing redundant windows. Outputs include:

- `candidate_review.tsv` and `first_batch_experiment_set.tsv`
- `experiment_candidates.tsv`, paired short-peptide, long-peptide, minigene,
  targeted-RNA, and manual-review tables
- weighted-versus-consensus ranking comparison
- APPM/HLA-LOH and CCF/clonality review summaries
- patient and technical reports, with optional DOCX/HTML/PPTX renderings when
  their document dependencies are installed
- integrity checks, hashes, reason codes, and a review run manifest

The review output directory must differ from the source result directory.

## Architecture boundary

The public macro Skills compose fine-grained Skills and production CLIs. They do not duplicate biological algorithms. `open-neo-run` never overwrites the weighted baseline, and `open-neo-review` never modifies Pipeline ranking outputs.

Skill3 can describe computational candidates, evidence gaps, experiment
priority, and research validation routes. It must not claim a confirmed
neoantigen, guaranteed benefit, clinical resistance, treatment failure, a drug
recommendation, or a finalized vaccine design.

Plan and dry-run use Input QC, Doctor and `pipeline-full`; approved
`open-neo-run` execute/resume requests are submitted to
NeoAg Gateway and dispatched to the production runner. Tool outputs are
cross-checked by evidence domain and written to `tool_run_status.tsv`,
`tool_consensus_summary.tsv`, domain consensus tables and conflict tables.

RNA evidence is first-class. Existing gene TPM, transcript TPM and RNA
alt/VAF tables are reused. With tumor RNA FASTQ and a declared Salmon index
plus tx2gene (or an RSEM reference), Skill2 generates both gene and transcript
TPM. With tumor RNA BAM plus somatic VCF, it generates RNA ref/alt reads, depth
and VAF. Fusion/splice junction-read fields remain separate event evidence.

## Gateway task contract

The existing compatible endpoints remain:

- `POST /open/install-check` (`task_type=open_install_check`)
- `POST /open/run` (`task_type=open_run`)
- `POST /open/review` (`task_type=open_review`)

Risk is assigned from the requested work: plan/dry-run and report-free review
are LOW; execution from precomputed evidence and document/PPT generation are
MEDIUM; install/repair and BAM/FASTQ or production-manifest execution are HIGH.
HIGH requests require `approved=true`. HPC submission is not part of these
public macro Skills.

Each accepted or approval-blocked request writes:

```text
work/neoag_gateway/jobs/<job_id>/
  request.json
  approval.json
  job_status.json
  audit_log.jsonl
```

The Gateway also retains its global audit log and legacy flat job-status JSON
for compatibility with existing clients.
