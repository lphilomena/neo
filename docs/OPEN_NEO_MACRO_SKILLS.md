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
`manifests/rna_fusion_splice.production.toml`. Its DAG covers FASTQ QC, STAR,
Salmon, EasyFuse, STAR-Fusion, Arriba, RegTools, optional SNAF/SpliceMutr,
candidate normalization, peptide generation and downstream evidence ranking.
Execution requires HLA, matching FASTA/GTF/STAR/CTAT/EasyFuse/Salmon assets and
Gateway approval. Missing optional cohort workflows remain `UNASSESSED`.

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
  --outdir reviews/CASE001
```

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
