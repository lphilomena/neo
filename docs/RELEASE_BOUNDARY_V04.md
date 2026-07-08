# V04 Release Boundary

This document defines the lightweight v04 release boundary for the current server checkout.

## Bundled in v04 lightweight release

- Python package source under `src/neoag_v03/`
- Project entry points under `bin/` that are small wrappers/scripts
- Profiles under `profiles/`
- Example configs under `conf/`
- Nextflow modules/workflows under `modules/` and `workflows/`
- Tests under `tests/`
- Small fixtures under `data/fixtures`, `data/fixtures_snv`, `data/fixtures_sv`, `assets`, and `resources`
- Public documentation under `docs/`
- Release notes/changelogs in the repository root
- Packaging/reproducibility scripts under `scripts/`

## Excluded from lightweight release

- Local virtual environments: `.venv/`, `.venv.local/`
- Nextflow/runtime caches: `.nextflow/`, `.nextflow.log*`, `work/`, `.pytest_cache/`
- Heavy external tools: `tools/`
- Conda packs and downloaded references: `conda_packs/`, `data/ref/`, `data/vep/`, large external archives
- Prior run outputs and benchmark outputs: `results/`
- Historical release tarballs under `dist/`
- Migration logs and ad hoc local logs

## External or optional requirements

Some tool integrations are intentionally optional for the lightweight release. The code should record missing optional tools without failing release tests. Current optional/may-be-external tools include NetMHCpan, PRIME, STAR-Fusion, and FusionCatcher, depending on workflow mode.

## Clinical-use boundary

The v04 release remains a computational prototype. It can rank and annotate neoantigen candidates, safety flags, and immune-escape hypotheses, but it does not produce clinical treatment recommendations.
