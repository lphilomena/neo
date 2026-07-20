# Changelog — v0.4.1 APPM 2.0 + CCF 2.0 + Immune Escape 2.0

## Added

- `src/neoag/appm_v2.py`
  - APPM gene status, pathway status and peptide-level APPM flags.
  - Biallelic status logic from damaging variant + CN/LOH/expression evidence.
- `src/neoag/ccf_v2.py`
  - Copy-number-aware CCF estimates with mutation multiplicity enumeration.
  - CCF min/max/best and confidence/method fields.
- Immune Escape 2.0 updates in `src/neoag/immune_escape.py`
  - APPM 2.0 and CCF 2.0 sidecar consumption.
  - Treatment context policy for `vaccine`, `tcr_target`, `immunomonitoring`, `discovery`.
- CLI commands:
  - `neoag appm-2`
  - `neoag ccf-2`
  - extended `neoag immune-escape` with APPM/CCF/context options.
- Tests:
  - `tests/test_v041_appm_ccf_escape.py`

## Changed

- `appm-lite` is now backed by APPM 2.0 while preserving the old return shape and output names.
- `ccf-lite` is now backed by CCF 2.0 while preserving legacy fields required by `score`.
- `run` and SV scoring flows pass CNV/raw peptide context into APPM and immune-escape evidence where available.

## Verified

```text
PYTHONPATH=src pytest -q
121 passed
```

## Interpretation boundary

The new modules provide computational evidence for antigen presentation capacity, clonality, and immune escape mechanisms. They do not provide clinical resistance or safety determinations.
