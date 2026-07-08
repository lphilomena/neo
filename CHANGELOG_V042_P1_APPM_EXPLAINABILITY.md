# Changelog: v0.4.2 P1 APPM explainability

## Added

- `src/neoag_v03/appm_explainability.py`
  - APPM call confidence calculation
  - MHC-I/MHC-II/IFNG submodule scoring
  - APPM output post-processing for v0.4.2 P1 sidecars

- `src/neoag_v03/reports_v041.py`
  - APPM evidence card
  - MHC-I submodule card
  - Top APPM driver defect card
  - Peptide mechanism card combining APPM, immune escape, safety and CCF

- `appm/appm_submodule_scores.tsv`

- `neoag-v03 report-v041`

- Extended `neoag-v03 immune-escape`
  - `--ranked-peptides`
  - `--top-priority-threshold`
  - affected candidate burden and top-candidate counts

- Extended `neoag-v03 benchmark-system`
  - `--ranked-peptides`
  - `--appm-summary`
  - `--appm-module-scores`
  - `--appm-submodule-scores`
  - `--peptide-appm-flags`
  - `--peptide-escape-flags`

- Benchmark outputs:
  - `appm_ms_stratified_validation.tsv`
  - `appm_multiplier_delta.tsv`
  - `hla_ligand_detection_by_appm.tsv`

## Changed

- `appm_summary.tsv`, `appm_module_scores.tsv`, `appm_pathway_status.tsv`, and `appm_peptide_modifiers.tsv` now include APPM confidence fields.
- `immune_escape_events.tsv` now includes escape event CCF context and affected peptide counts.
- `immune_escape_summary.tsv` now includes escape burden summary fields.
- Project version updated to `0.4.2`.

## Tests

Added:

- `tests/test_v042_appm_explainability.py`
- `tests/test_v042_immune_escape_burden.py`
- `tests/test_v042_report_and_benchmark.py`

Validated subset:

```bash
PYTHONPATH=src pytest -q \
  tests/test_v041_appm_ccf_escape.py \
  tests/test_v042_appm_explainability.py \
  tests/test_v042_immune_escape_burden.py \
  tests/test_v042_report_and_benchmark.py \
  tests/test_v04_evidence_safety_escape.py
```

Result: `18 passed`.
