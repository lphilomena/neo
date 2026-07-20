# Project B v0.4.1 — APPM 2.0 + CCF 2.0 + Immune Escape 2.0

This release upgrades the v0.4 evidence-safety-escape layer with three linked evidence modules.

## Scope

v0.4.1 remains a computational triage system. It provides antigen presentation machinery evidence, copy-number-aware clonality estimates, and immune-escape mechanism flags. It does not make a clinical resistance diagnosis.

## APPM 2.0

New module: `src/neoag/appm_v2.py`.

Outputs:

```text
appm/appm_gene_status.tsv
appm/appm_pathway_status.tsv
appm/peptide_appm_flags.tsv
appm/appm_summary.tsv
```

Key improvements:

```text
- Gene-level status for B2M, HLA-A/B/C, TAP1/2, NLRC5, CIITA, JAK/STAT/IFNGR and MHC-II regulators.
- Biallelic-status logic using damaging variant + CN/LOH/expression evidence.
- Pathway-level states: MHC_I_INTACT/CAUTION/DEFECTIVE/UNASSESSED, MHC_II_*, IFNG_RESPONSE_*.
- Peptide-level APPM flags and priority caps for B2M loss, lost restricting HLA, TAP defects and CIITA defects.
```

CLI:

```bash
neoag appm-2 \
  --sample-id P001 \
  --vep-tsv vep_appm.tsv \
  --expression expression.tsv \
  --cnv cnv_gene_status.tsv \
  --hla-loh hla_loh.tsv \
  --raw-peptides parsed/raw_peptides.tsv \
  --outdir appm
```

The legacy `appm-lite` command is backward compatible and now also writes APPM 2.0 sidecars.

## CCF 2.0

New module: `src/neoag/ccf_v2.py`.

Output:

```text
clonality/ccf_lite.tsv
```

The table keeps legacy fields (`ccf_estimate`, `ccf_status`, `clonality_multiplier`) and adds:

```text
ccf_best
ccf_min
ccf_max
ccf_ci_low
ccf_ci_high
multiplicity_best
total_cn
major_cn
minor_cn
loh_status
ccf_method
ccf_confidence
ccf_warning
```

Methods:

```text
SNV_INDEL_COPY_NUMBER_AWARE
SV_BREAKPOINT_APPROX
JUNCTION_APPROX
RNA_ONLY_UNRESOLVED
COPY_NUMBER_AWARE_APPROX
```

CLI:

```bash
neoag ccf-2 \
  --events parsed/raw_events.tsv \
  --purity purity.tsv \
  --cnv cnv_segments.tsv \
  --out clonality/ccf_lite.tsv
```

## Immune Escape 2.0

Updated module: `src/neoag/immune_escape.py`.

Outputs:

```text
immune_escape/immune_escape_events.tsv
immune_escape/immune_escape_summary.tsv
immune_escape/peptide_escape_flags.tsv
```

Key improvements:

```text
- Consumes APPM 2.0 sidecars when supplied.
- Supports treatment context: vaccine, tcr_target, immunomonitoring, discovery.
- Vaccine/TCR context: lost restricting HLA hard-caps the peptide at D.
- Immunomonitoring/discovery context: lost HLA is retained but flagged and capped for review.
- B2M biallelic loss hard-caps MHC-I peptides.
- TAP/NLRC5/CIITA/JAK-STAT defects are mapped to peptide-level multipliers and caps.
```

CLI:

```bash
neoag immune-escape \
  --sample-id P001 \
  --raw-peptides parsed/raw_peptides.tsv \
  --appm-gene-status appm/appm_gene_status.tsv \
  --appm-pathway-status appm/appm_pathway_status.tsv \
  --ccf clonality/ccf_lite.tsv \
  --hla-loh hla_loh.tsv \
  --therapy-context vaccine \
  --outdir immune_escape
```

## End-to-end integration

The existing `run`, `sv-run-full`, and `sv-run-full-wes` flows now keep compatibility with v0.4 output while producing richer APPM/CCF/immune-escape evidence where inputs are available.

Final scoring still consumes:

```text
peptide_safety.tsv
peptide_escape_flags.tsv
ccf_lite.tsv
appm_summary.tsv
```

Hard interpretation boundaries:

```text
APPM defect = antigen-presentation mechanism evidence, not absolute absence of immune response.
CCF = bulk sequencing inference, not single-cell truth.
Immune escape risk = mechanism evidence, not clinical resistance diagnosis.
```
