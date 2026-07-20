# v0.4 Evidence / Safety / Escape Patch

Implemented in this package:

- Capture-aware WES SV Phase 1.5:
  - `src/neoag/sv/wes_capture.py`
  - capture BED parsing, expanded BED sidecars, breakend capture status, WES confidence tiers, priority caps.
  - CLI options for `sv-build-raw-wes` and `sv-run-full-wes`: `--capture-bed`, `--capture-near-bp`, `--capture-slop-bp`.

- Peptide safety gate:
  - `src/neoag/peptide_safety_gate.py`
  - reference proteome exact match, normal ligandome match, normal junction match, anchor-only risk.
  - output sidecars: `safety/peptide_safety.tsv`, `safety/event_safety.tsv`.

- Immune escape / HLA LOH layer:
  - `src/neoag/immune_escape.py`
  - peptide-level lost restricting HLA flags, B2M/JAK/APM/CIITA risk summaries.
  - output sidecars: `immune_escape/immune_escape_events.tsv`, `immune_escape/immune_escape_summary.tsv`, `immune_escape/peptide_escape_flags.tsv`.

- Scoring integration:
  - `score` now accepts `--peptide-safety` and `--peptide-escape-flags`.
  - `sv-run-full` and `sv-run-full-wes` build and consume safety/escape evidence.

- Release fixture repair:
  - Added lightweight LOHHLA and BigMHC_IM fixtures so the lite package test suite is self-contained.

Validation:

```bash
PYTHONPATH=src pytest -q
# 118 passed
```

Boundaries:

- WES SV Phase 1.5 is capture-limited and should not be interpreted as WGS-equivalent SV discovery.
- Peptide safety gate reduces computational off-target risk but does not prove clinical safety.
- Immune escape output is mechanism/risk evidence and is not a clinical drug-resistance diagnosis.
