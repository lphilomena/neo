# Source-chain C1–C4 Test Report

## Scope

This release adds event-specific candidate source-chain confidence for SNV, InDel, Fusion and Splice candidates, while keeping C1–C4 independent from the downstream R1–R4 recommendation grade.

## Validation performed

- Python byte-code compilation: `python -m compileall -q src` — PASS.
- Focused source-chain, report-dimension, splice and report tests — 57 passed, 1 skipped.
- Complete default regression suite — 482 passed, 106 skipped.
- Regression case fixed: `NOT_CONFIRMED`, `NOT_PERFORMED` and `NOT_TESTED` are treated as `UNASSESSED`, not as positive orthogonal confirmation.

## Focused test command

```bash
PYTHONPATH=src pytest -q \
  tests/test_source_chain_confidence.py \
  tests/test_report_dimension_audit.py \
  tests/test_splice_v050_layer.py \
  tests/test_v042_report_and_benchmark.py
```

Result:

```text
57 passed, 1 skipped
```

Complete default suite:

```text
482 passed, 106 skipped
```

## Production boundary

The package does not include licensed tools, patient data, large reference assets, model weights or production HPC configuration. External-tool and large-reference validation must still be performed on the target deployment machine with `open-neo install-check` / NeoAg Doctor and the relevant smoke tests.

The integrated EC-v3 source-chain profile remains `PROVISIONAL_RESEARCH_ONLY`. The compatibility profile is recommended for historical-case review because it emits C1–C4 fields without changing the existing R1–R4/Pareto ranking.
