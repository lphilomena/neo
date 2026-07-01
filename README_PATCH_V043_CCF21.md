# Project B v0.4.3 CCF 2.1 P0-P1 patch

Copy this patch over a Project B v0.4.2 P1 tree. It updates CCF 2.0 to CCF 2.1 by adding:

- `ccf_input_qc.tsv`
- multiplicity confidence fields
- clonality confidence fields
- optional `--external-clonality` and `--svclone` inputs
- `ccf_conflicts.tsv`
- `ccf_cluster.tsv`
- event-type-aware handling of SNV/InDel, SV, WES SV, junction and RNA-only events

Targeted validation used during packaging:

```bash
PYTHONPATH=src pytest -q \
  tests/test_v041_appm_ccf_escape.py \
  tests/test_v042_appm_explainability.py \
  tests/test_v042_immune_escape_burden.py \
  tests/test_v043_ccf21.py \
  tests/test_v04_evidence_safety_escape.py
```

Result: `20 passed`.
