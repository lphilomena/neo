# Open-Neo failure codes

Failure codes are stable machine-readable identifiers. Human-readable details
belong in the step `detail` field. Code definitions are owned by
`src/neoag/open_neo/errors.py`.

## Installation and release

- `CHECKSUM_REQUIRED`: a supplied release archive has no expected SHA256.
- `CHECKSUM_FAILED`: archive missing or observed SHA256 differs.
- `CORE_INSTALL_FAILED`: portable installer could not complete.
- `INSTALL_TIMEOUT`: portable installer exceeded the configured timeout.
- `RELEASE_STAGING_FAILED`: archive is unsafe, unsupported, or has no unique project root.
- `ASSET_SOURCE_UNCONFIGURED`: required assets cannot be reached from the declared local root or remote host.
- `PROJECT_ROOT_INVALID`: selected source/staged directory is not an Open-Neo project root.
- `TOOL_MISSING` / `TOOL_SMOKE_FAILED`: tool entrypoint or mini-smoke failed.
- `REFERENCE_MISSING` / `REFERENCE_HASH_MISMATCH`: required reference absent or changed.
- `LICENSE_BLOCKED`: licensed software/assets cannot be installed or redistributed automatically.
- `PRIVATE_PATH_DETECTED` / `PATIENT_DATA_IN_RELEASE`: release boundary is unsafe.

## Run and Gateway

- `AMBIGUOUS_INPUT`, `ROUTE_FAILED`, `HLA_MISSING`, `GENOME_BUILD_MISMATCH`.
- `APPROVAL_REQUIRED`: high-risk execution lacks explicit approval.
- `GATEWAY_REQUIRED`: direct heavy execution is forbidden.
- `DOCTOR_BLOCKED`: required production dependency failed Doctor.
- `PRODUCTION_MANIFEST_REQUIRED`: raw production inputs have no controlled profile.
- `RNA_FUSION_SPLICE_MISSING:<field>`: automatic RNA profile lacks a required asset.
- `PIPELINE_STAGE_FAILED`, `CONSENSUS_RANKING_FAILED`.

## Review

- `NEEDS_RANKING`: event-level evidence-consensus output is absent.
- `REVIEW_INTEGRITY_BLOCKED`: run/hash/evidence integrity checks failed.
- `EVENT_MAPPING_FAILED`: event representatives cannot map to peptide rows.
- `REPORT_BOUNDARY_VIOLATION`: review output would overwrite source results.
- `INVALID_REPORT_SELECTION`: use only `patient`, `technical`, `onepage`, or `none`.

`APPROVAL_REQUIRED` uses CLI exit code 3; `NEEDS_RANKING` uses 4; other
terminal failures use 2. Gateway may map approval to HTTP 403 and state or
integrity conflicts to 409/422.
