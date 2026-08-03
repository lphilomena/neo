# Changelog — Candidate Source-chain Confidence C1–C4

## v0.5.2 / source-chain-v1.0

### Added

- Event-specific SNV, InDel, Fusion and Splice source-chain evaluators.
- C1–C4 candidate source-chain confidence independent of R1–R4 recommendation grade.
- Explicit `NOT_APPLICABLE`, `UNASSESSED`, `INDETERMINATE_LOW_POWER`, `SUPPORTED`, `NEGATIVE` and `CONFLICT` statuses.
- Independent/cross-modal orthogonal confirmation logic; multiple callers on the same BAM are not treated as orthogonal confirmation.
- Compatibility profile that emits source-chain fields without changing legacy EC ranking.
- Integrated translational profile with C4 hard fail, C3 R3 cap, optional C2 R2 cap and source-chain Pareto dimension.
- Long-format requirement audit and source-chain summary outputs.
- `neoag source-chain` and `neoag report-dimension-audit` CLI commands.
- Machine-readable report-dimension map and technical report source-chain fields.
- Unit and regression tests for all four tracks and compatibility/integrated modes.
- Fixed orthogonal-status parsing so `NOT_CONFIRMED`, `NOT_PERFORMED` and `NOT_TESTED` remain `UNASSESSED` and cannot be misread as positive confirmation.
- Added explicit SNV base/MAP quality, FFPE, low-complexity and paralog QC.
- Added coverage-aware normal-junction interpretation; unqualified `NOT_DETECTED` is low-power rather than normal-negative evidence.
- Removed the automatic RNA-only Fusion R3 cap for complete C1/C2 source chains while retaining caps for incomplete chains.
- Added Fusion duplicate/junction-uniqueness evidence and Splice PSI, normal-isoform, translation-direction and NMD requirements.
- Added stable requirement-level applicability/value/source/conflict fields and separate C-tier/R-tier report presentation.

### Boundary

C1–C4 establish traceability from event to transcript/ORF and peptide. They do not establish natural HLA presentation, T-cell recognition or clinical benefit.
