# Controlled Execution Phase 0-2 Changelog Changelog

## Added

- `src/neoag/controlled_execution/`
  - Phase 0: `doctor.py`, `release_audit.py`
  - Phase 1: `gateway.py`
  - Phase 2: `pipeline_runner.py`
  - shared manifest, audit and IO helpers
- Console scripts:
  - `neoag-doctor`
  - `neoag-gateway`
  - `neoag-pipeline-full`
  - `neoag-release-audit`
- Example manifests under `configs/controlled_execution/`.
- Documentation: `docs/CONTROLLED_EXECUTION_PHASE0_2.md`.
- Tests: `tests/test_controlled_execution_phase0_2.py`.

## Safety

- Doctor is read-only.
- Gateway blocks high-risk execution unless `approved=true` is supplied.
- Pipeline-full defaults to dry-run planning.
- Release audit scans for cache/runtime artifacts and private/patient path hints.
