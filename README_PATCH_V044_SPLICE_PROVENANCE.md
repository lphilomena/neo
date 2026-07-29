# v0.4.4 Patch Guide

This patch can be applied over NeoAg v0.4.3. The full release archive is preferred because schema, adapters, production merge logic, tests, documentation, and version metadata change together.

Primary changed files:

```text
src/neoag/splice/coordinates.py
src/neoag/splice/registry.py
src/neoag/splice/normalization.py
src/neoag/provenance.py
src/neoag/adapters/event_catalog.py
src/neoag/adapters/pvactools_parser.py
src/neoag/adapters/splice_junction_adapter.py
src/neoag/open_neo/tool_consensus.py
src/neoag/open_neo/rna_fusion_splice_profile.py
src/neoag/production_runner.py
src/neoag/schemas.py
scripts/normalize_rna_fusion_splice.py
tests/test_splice_v044_provenance.py
```

Verification:

```bash
PYTHONPATH=src pytest -q tests/test_splice_v044_provenance.py
PYTHONPATH=src pytest -q tests/test_splice_junction.py --run-integration
```

The release archive also includes all compatibility tests and documentation.

Detailed validation results are recorded in `VALIDATION_V044_SPLICE_PROVENANCE.md`.
