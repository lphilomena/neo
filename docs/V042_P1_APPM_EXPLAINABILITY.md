# Project B v0.4.2 P1: APPM explainability and mechanism-evidence upgrade

This release moves the APPM/CCF/immune-escape layer from “has a result” to “has a reviewable mechanism explanation”. It remains a computational triage system, not a clinical resistance or safety diagnosis.

## Main additions

### 1. APPM call confidence

APPM outputs now include:

- `appm_call_confidence`: `high`, `medium`, `low`, or `insufficient`
- `appm_call_confidence_score`: numeric 0–1 score
- `confidence_reason`: machine-readable reason codes
- `critical_missing_evidence`
- `evidence_conflict_impact`

These fields are appended to:

- `appm/appm_summary.tsv`
- `appm/appm_module_scores.tsv`
- `appm/appm_pathway_status.tsv`
- `appm/appm_peptide_modifiers.tsv`

Interpretation:

- DNA/CNV/LOH-supported biallelic loss is stronger evidence than RNA-only low expression.
- RNA-only low expression cannot become a high-confidence defect claim.
- Missing CNV/RNA/HLA-LOH evidence lowers confidence and is not treated as intact APPM.

### 2. MHC-I submodule scoring

New output:

```text
appm/appm_submodule_scores.tsv
```

The table currently includes:

- `MHC_I_CORE`: B2M and HLA-A/B/C
- `MHC_I_PROCESSING`: TAP1/2, TAPBP, ERAP1/2, PSMB8/9/10, chaperones
- `MHC_I_REGULATION`: NLRC5
- `MHC_I_HLA_LOH`: allele-level HLA LOH burden
- `MHC_II_CORE`: CIITA/RFX/HLA-II core genes
- `IFNG_SIGNALING`: IFNGR/JAK/STAT/IRF signaling

Hard rules are preserved:

- B2M biallelic loss strongly caps MHC-I interpretation.
- HLA LOH is allele-specific: peptide-level rejection/capping is handled by immune escape, not by globally zeroing all MHC-I peptides.

### 3. Immune escape burden and escape CCF context

`immune_escape_events.tsv` now includes:

- `escape_event_source`
- `escape_event_source_event_ids`
- `escape_event_ccf_best/min/max`
- `escape_event_ccf_confidence`
- `escape_event_clonality`
- `affected_candidate_count`
- `affected_top_candidate_count`
- `affected_hla_alleles`
- `affected_event_ids`
- `affected_peptide_ids`

`immune_escape_summary.tsv` now includes:

- `n_peptides_affected_by_hla_loh`
- `n_top_peptides_affected_by_hla_loh`
- `n_mhc_i_peptides_affected_by_b2m`
- `n_top_mhc_i_peptides_affected_by_b2m`
- `escape_burden_summary`

If an escape event cannot be mapped to a CCF-bearing source event, it is marked as `unresolved` / `not_mapped`; no pseudo-CCF is fabricated.

### 4. APPM-stratified ligandome/MS benchmark harness

`benchmark-system --mode ligandome-ms` now additionally writes:

- `appm_ms_stratified_validation.tsv`
- `appm_multiplier_delta.tsv`
- `hla_ligand_detection_by_appm.tsv`

If no external ligandome/MS data is provided, these outputs explicitly report `external_required`. MS-positive evidence is treated as strong presentation evidence; MS-negative is not treated as evidence of absence.

### 5. v0.4.2 HTML report

New CLI:

```bash
neoag report-v041 \
  --ranked-events scoring/ranked_events.tsv \
  --ranked-peptides scoring/ranked_peptides.tsv \
  --appm-summary appm/appm_summary.tsv \
  --appm-gene-status appm/appm_gene_status.tsv \
  --appm-module-scores appm/appm_module_scores.tsv \
  --appm-submodule-scores appm/appm_submodule_scores.tsv \
  --appm-peptide-modifiers appm/appm_peptide_modifiers.tsv \
  --immune-escape-summary immune_escape/immune_escape_summary.tsv \
  --peptide-escape-flags immune_escape/peptide_escape_flags.tsv \
  --peptide-safety safety/peptide_safety.tsv \
  --ccf clonality/ccf_lite.tsv \
  --out reports/evidence_report.v041.html
```

The report includes sample-level APPM cards, submodule scores, top APPM defects, immune escape burden, and peptide mechanism cards.

## CLI additions

### Immune escape affected-top counts

```bash
neoag immune-escape \
  --sample-id P001 \
  --raw-peptides parsed/raw_peptides.tsv \
  --ranked-peptides scoring/ranked_peptides.tsv \
  --top-priority-threshold B_CAUTION \
  --appm-gene-status appm/appm_gene_status.tsv \
  --appm-pathway-status appm/appm_pathway_status.tsv \
  --ccf clonality/ccf_lite.tsv \
  --hla-loh hla_loh.tsv \
  --therapy-context vaccine \
  --outdir immune_escape
```

### APPM/MS benchmark harness

```bash
neoag benchmark-system \
  --mode ligandome-ms \
  --ligandome-ms ligandome_ms.tsv \
  --ranked-peptides scoring/ranked_peptides.tsv \
  --appm-summary appm/appm_summary.tsv \
  --appm-module-scores appm/appm_module_scores.tsv \
  --appm-submodule-scores appm/appm_submodule_scores.tsv \
  --peptide-escape-flags immune_escape/peptide_escape_flags.tsv \
  --outdir benchmark
```

## Boundary

- APPM defect means antigen-presentation mechanism evidence, not a clinical resistance diagnosis.
- CCF is a bulk-sequencing inference, not single-cell truth.
- HLA LOH/B2M/JAK/APM evidence can cap, reject, or flag candidate peptides depending on therapy context, but cannot alone determine patient clinical resistance.
