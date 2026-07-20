# Multi-Entry Input Architecture (Project B)

Project B is an **event-level, multi-source neoantigen system**. It does not assume a single input type (e.g. pVACseq-only). Instead, diverse upstream callers are normalized into a **standard intermediate layer**, then scored with unified Layer-3 immunology ranking.

## Standard intermediate layer

| File | Role |
|------|------|
| `parsed/raw_events.tsv` | Layer-1/2 events (mutation source + peptide consequence) |
| `parsed/raw_peptides.tsv` | Candidate neoantigen peptides with HLA |
| `presentation/presentation_evidence.tsv` | HLA binding + presentation (Layer 3) |
| `parsed/expression_evidence.tsv` | Gene/event expression support |
| `parsed/rna_junction_evidence.tsv` | RNA junction, targeted fusion/splice validation, RNA alt reads, and RNA VAF support |
| `parsed/fusion_evidence.tsv` | EasyFuse / fusion caller structured evidence (optional) |
| `clonality/ccf_lite.tsv` | Clonality / CCF from VAF + purity + CNV |
| `safety/safety_evidence.tsv` | Normal-tissue / ligand safety pre-score |

## Multi-entry modes (`inputs.entry_mode`)

| Mode | Key | Typical inputs |
|------|-----|----------------|
| **A** | `snv_indel` | Annotated tumor/normal VCF → pVACseq; HLA; expression |
| **B** | `fusion` | EasyFuse `fusions.pass.csv` / STAR-Fusion / Arriba / AGFusion + pVACfuse (optional) |
| **C** | `splice_junction` | Annotated VCF + RegTools junction TSV + RNA support |
| **D** | `sv` | SV/BND VCF + GTF + reference (`sv-build-raw` output) |
| **E** | `peptide_only` | Peptide CSV/TSV + HLA (+ optional external evidence) |
| **F** | `e2e` | BAM/FASTQ end-to-end via `[tools.enabled]` upstream |
| — | `intermediates` | Pre-built `raw_events` + `raw_peptides` passthrough |

## CLI

```bash
# Build raw layer only (any entry mode)
neoag build-intermediates --config conf/run.example.toml --outdir results/P001/intermediates

# Build evidence sidecars from existing raw tables
neoag build-evidence-layer --outdir results/P001 --profile leukemia \
  --rna-vaf results/P001/parsed/rna_vaf.tsv \
  --rna-junction results/P001/parsed/rna_junctions.tsv

# Score from pre-built raw tables (no pVAC re-parse)
neoag run --outdir results/P001 \
  --raw-events results/P001/parsed/raw_events.tsv \
  --raw-peptides results/P001/parsed/raw_peptides.tsv \
  --netmhcpan ... --expression ...

# Full upstream + scoring (uses entry_mode from TOML)
neoag run-full --config conf/run.example.toml --outdir results/P001
```

## SV / SNV dedicated flows

- **SV**: `sv-run-full` / `sv-run-full-wes` → copies adapter output to standard `parsed/` layout
- **SNV WES**: `snv-run-full-wes` → Mutect2 + upstream + `run`

These remain first-class entry paths; set `entry_mode = "sv"` or use `inputs.sv_raw_events` when composing multi-source runs.

## EasyFuse (Mode B — RNA fusion discovery adapter)

EasyFuse is **not** a full neoantigen system. Use it as fusion discovery + RNA evidence input:

```text
WTS RNA-seq FASTQ → EasyFuse → fusions.pass.csv
  → EasyFuseAdapter → raw_events / raw_peptides / fusion_evidence
  → HLA binding + safety + score
```

Configure in TOML:

```toml
[inputs]
entry_mode = "fusion"
easyfuse_pass_csv = "/path/to/sample/fusions.pass.csv"
# or fusion_tsv = ".../fusions.pass.csv"   # auto-detected by BPID header
```

CLI:

```bash
neoag build-intermediates --entry-mode fusion \
  --easyfuse-tsv data/fixtures/easyfuse_fusions.pass.tsv \
  --outdir results/EF1/intermediates
```

Project B applies additional filters beyond `prediction_class=positive`:

- junction reads ≥ 3, anchor size ≥ 10
- `frame != no_frame`, non-empty `neo_peptide_sequence`
- excludes `cis_near` read-through fusions by default

Peptide rows from `neo_peptide_sequence` are stubs (HLA filled by upstream netMHCpan/MHCflurry). Prefer pVACfuse for full HLA-typed epitopes when available.

## Merging sources

`build-intermediates` can merge **pVAC outputs + fusion catalog + splice catalog** into one `raw_events.tsv`. Peptides come from pVAC/SV/peptide adapters; event-only catalogs (fusion/splice stubs) support event-level ranking before peptide calling completes.

## RNA and HLA LOH evidence additions

`build-evidence-layer` accepts `--rna-vaf` for RNA allele-count tables. The output appends `rna_alt_reads`, `rna_ref_reads`, `rna_depth`, `rna_vaf`, `rna_support_status`, and targeted fusion/splice validation fields to `parsed/rna_junction_evidence.tsv`.

For immune-escape inputs, convert LOHHLA and SpecHLA independently, then cross-check them before passing the consensus table downstream:

```bash
neoag convert-lohhla -i LOHHLA.HLAlossPrediction_CI.xls -o results/P001/tools/lohhla.hla_loh.tsv
neoag convert-spechla -i merge.hla.copy.txt -o results/P001/tools/spechla.hla_loh.tsv
neoag crosscheck-hla-loh \
  --lohhla-hla-loh results/P001/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/P001/tools/spechla.hla_loh.tsv \
  --out results/P001/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/P001/tools/hla_loh.consensus.tsv
```
