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
  --top-n 12 \
  --outdir reviews/CASE001
```

The weighted baseline is optional for review. When present it enables ranking
comparison. R1/R2 events may enter the first experimental batch; R3 events are
written to a separate evidence-completion queue. The review output directory
must differ from the source result directory.

## Architecture boundary

The public macro Skills compose fine-grained Skills and production CLIs. They do not duplicate biological algorithms. `open-neo-run` never overwrites the weighted baseline, and `open-neo-review` never modifies Pipeline ranking outputs.

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
