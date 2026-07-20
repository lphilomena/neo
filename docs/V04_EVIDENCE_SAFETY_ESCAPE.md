# Project B v0.4: Evidence, Safety, and Escape Layer

This release extends `neoag_event_pipeline` with three production-facing evidence layers:

1. **WES SV Phase 1.5**: capture-aware WES-limited SV/Junction candidate triage.
2. **Peptide safety gate**: reference proteome, normal ligandome, normal junction, WT/anchor-risk filters.
3. **Immune escape / HLA LOH layer**: peptide-level loss-of-restricting-HLA and global MHC/APM escape flags.

These modules provide computational risk stratification. They are not clinical safety certification and do not diagnose clinical drug resistance.

## 1. WES SV Phase 1.5

WES SV mode is explicitly treated as `EXOME_CAPTURE_LIMITED`; it is not equivalent to WGS SV discovery.

### Inputs

Required for WES mode:

```text
--sv-vcf <Manta/SvABA/etc VCF>
--reference-fasta GRCh38.fa
--gencode-gtf gencode.gtf
--hla hla.txt
--capture-bed exome_capture.bed
```

Recommended:

```text
--expression expression.tsv
--rna-junctions rna_junctions.tsv
--normal-expression normal_expression.tsv
--normal-hla-ligands normal_hla_ligands.tsv
--reference-proteome human_proteome.fa
--normal-junctions normal_junctions.tsv
--hla-loh hla_loh.tsv
--cnv facets_or_cnvkit_segments.tsv
```

### Capture-aware fields

WES SV output sidecars now include:

```text
breakend1_capture_status
breakend2_capture_status
breakend1_capture_distance_bp
breakend2_capture_distance_bp
capture_interpretability
wes_confidence_tier
priority_cap
filter_status
filter_reason
```

### Tiers and caps

```text
WES_Tier1: RNA junction supported and capture-interpretable; priority cap B
WES_Tier2: coding/capture evidence but no strong RNA junction; priority cap B_CAUTION
WES_Tier3: weak WES-only evidence; priority cap C
WES_UNINTERPRETABLE: off-target or insufficient evidence; priority cap D
```

## 2. Peptide safety gate

The new safety sidecar is written as:

```text
safety/peptide_safety.tsv
safety/event_safety.tsv
```

It evaluates:

```text
reference_proteome_exact_match
normal_hla_ligand_exact_match
normal_ligand_tissue
normal_junction_seen
WT peptide binding / anchor-only mutation risk
critical normal tissue expression
```

Hard safety failures are carried into `score` through `--peptide-safety`.

## 3. Immune escape / HLA LOH

The immune escape layer writes:

```text
immune_escape/immune_escape_events.tsv
immune_escape/immune_escape_summary.tsv
immune_escape/peptide_escape_flags.tsv
```

It evaluates:

```text
restricting HLA allele loss
B2M loss or damaging alteration
JAK1/JAK2 IFNγ-pathway risk
TAP1/TAP2 antigen processing risk
NLRC5 / CIITA transcriptional presentation risk
```

A peptide restricted by a lost HLA allele is marked with `escape_status=ESCAPE_REJECT` and is capped/rejected during scoring.

## 4. CLI examples

### Build WES SV raw inputs

```bash
PYTHONPATH=src python -m neoag.cli sv-build-raw-wes \
  --sample-id P001 \
  --profile sv_wes_phase1_5 \
  --sv-vcf wes_manta.somaticSV.vcf.gz \
  --callers Manta \
  --reference-fasta GRCh38.fa \
  --gencode-gtf gencode.gtf \
  --hla hla.txt \
  --capture-bed exome_capture.bed \
  --expression expression.tsv \
  --rna-junctions rna_junctions.tsv \
  --outdir results/P001_wes_sv_phase1_5
```

### Full WES SV scoring with safety and escape layers

```bash
PYTHONPATH=src python -m neoag.cli sv-run-full-wes \
  --sample-id P001 \
  --profile sv_wes_phase1_5 \
  --sv-vcf wes_manta.somaticSV.vcf.gz \
  --callers Manta \
  --reference-fasta GRCh38.fa \
  --gencode-gtf gencode.gtf \
  --hla hla.txt \
  --capture-bed exome_capture.bed \
  --expression expression.tsv \
  --rna-junctions rna_junctions.tsv \
  --normal-expression normal_expression.tsv \
  --normal-hla-ligands normal_hla_ligands.tsv \
  --reference-proteome human_proteome.fa \
  --normal-junctions normal_junctions.tsv \
  --hla-loh hla_loh.tsv \
  --cnv cnv_segments.tsv \
  --outdir results/P001_wes_sv_phase1_5
```

### Score from precomputed evidence sidecars

```bash
PYTHONPATH=src python -m neoag.cli score \
  --raw-events parsed/raw_events.tsv \
  --raw-peptides parsed/raw_peptides.tsv \
  --presentation presentation/presentation_evidence.tsv \
  --peptide-safety safety/peptide_safety.tsv \
  --peptide-escape-flags immune_escape/peptide_escape_flags.tsv \
  --out-events scoring/ranked_events.tsv \
  --out-peptides scoring/ranked_peptides.tsv
```

## 5. Interpretation boundaries

- WES SV Phase 1.5 is a capture-limited hypothesis-generation module.
- Peptide safety filtering reduces computational off-target risk but does not prove clinical safety.
- Immune escape output reports mechanism evidence and risk; it is not a clinical resistance diagnosis.
