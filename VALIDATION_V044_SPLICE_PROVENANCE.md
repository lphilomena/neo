# NeoAg v0.4.4 Validation Report

Validation date: 2026-07-29

## Scope

This validation covers the v0.4.4 exact-junction normalization, evidence non-leakage, provenance-preserving merge, exact splice consensus, compatibility projection, CLI, and package build.

## Implemented invariants

1. Canonical junction identity is `SJ|BUILD|CHR|INTRON_START_1BASED|INTRON_END_1BASED|STRAND`.
2. Verified `rna_junction_reads` can be transferred only by exact canonical identity, exact coordinates including build and strand, an exact unique source alias, or an explicit unique variant-to-junction relation.
3. Gene-level maximum reads, nearest-locus matching, approximate coordinate windows, and ambiguous aliases do not transfer verified reads.
4. Caller-provided unresolved counts remain in `provided_rna_junction_reads` and contribute zero to verified support.
5. Every merged event or peptide retains one long-form provenance row per source record.
6. Cross-domain splice confirmation requires the same canonical junction; same-gene but different-junction records are not confirmed.
7. A normal reference with insufficient locus-level coverage is not interpreted as definitive absence.

## Test results

### Python syntax

```text
python -m compileall -q src scripts tests
PASS
```

### Complete default test inventory

All 503 collected default tests were executed in one default test run on the target server:

```text
399 passed
104 skipped
0 failed
503 total
```

The skipped cases are tests intentionally guarded by optional integration, external-tool, licensed-resource, or environment conditions. No collected default test failed.

### Splice integration suite

```text
PYTHONPATH=src pytest -q --run-integration \
  tests/test_splice_junction.py \
  tests/test_splice_v044_provenance.py

21 passed
0 failed
```

### New v0.4.4 regression suite

The dedicated tests cover:

- RegTools annotated zero-based coordinate conversion;
- RegTools BED12 block-derived intron conversion;
- same-gene 5000-read versus 2-read junction non-leakage;
- ambiguous source alias rejection;
- separation of caller-provided and verified junction counts;
- multi-tool provenance retention;
- prevention of false same-gene cross-domain confirmation;
- exact pVACsplice-to-RegTools support assignment;
- production deduplication with full source provenance.

```text
8 passed
0 failed
```

### CLI and demo smoke tests

```text
neoag.__version__ = 0.4.4
normalize_rna_fusion_splice.py --help = PASS
neoag run-demo = PASS
```

The demo generated raw event/peptide tables, presentation, APPM, CCF, safety, immune-escape, weighted ranking, evidence-consensus ranking, validation plan, and technical/patient reports.

### Wheel build and installation

```text
neoag_event_pipeline-0.4.4-py3-none-any.whl
build = PASS
isolated target installation = PASS
import neoag = PASS
canonical junction API smoke = PASS
```

## Validation boundary

The regression tests use project fixtures and local stubs. They validate coordinate semantics, merge behavior, provenance, schemas, routing, and pipeline compatibility. They do not constitute biological or clinical validation of RegTools, SNAF, SpliceMutr, pVACtools, HLA predictors, or patient-specific neoantigen calls. External licensed binaries and large reference databases must still be validated in the deployment environment.
