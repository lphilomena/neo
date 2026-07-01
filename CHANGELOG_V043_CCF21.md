# Changelog: v0.4.3 CCF 2.1 P0-P1

## Added

- CCF input QC sidecar: `ccf_input_qc.tsv`.
- Multiplicity candidates and confidence fields.
- Clonality confidence fields.
- Coarse CCF interval probability summaries.
- Event-type-aware methods including `WES_SV_CAPTURE_LIMITED_APPROX` and `RNA_ONLY_UNRESOLVED`.
- Optional external clonality input for PyClone-VI / PhylogicNDT-like outputs.
- Optional SVclone-like input for SV CCF evidence.
- CCF conflict sidecar: `ccf_conflicts.tsv`.
- CCF cluster sidecar: `ccf_cluster.tsv`.
- New tests: `tests/test_v043_ccf21.py`.

## Changed

- `ccf-2` CLI now exposes `--external-clonality`, `--svclone`, and sidecar output options.
- `build_ccf_2` remains backward-compatible with old call sites while writing sidecars by default.
- WES SV and RNA-only events now have explicit method/confidence behavior.

## Boundary

CCF 2.1 is a clonality evidence layer for computational triage. It does not replace dedicated clonal deconvolution tools and does not produce clinical clonality diagnoses.
