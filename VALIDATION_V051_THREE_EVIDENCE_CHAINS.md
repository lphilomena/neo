# NeoAg v0.5.1 validation report

## Release scope

v0.5.1 extends the formal v0.5.0 Splice Provenance Layer with:

- RNA-driven dual-generator evidence: ImmunoPepper + moPepGen;
- DNA-causal evidence: splice2neo + EasyQuant + pVACsplice;
- normal-background evidence: coverage-aware normal sources + k4neo;
- exact query-ID and exact variant/junction/event provenance;
- multi-pass production orchestration and auditable external-tool manifests.

## Verification results

| Check | Result |
|---|---:|
| Full Python test suite | 450 passed, 105 skipped, 0 failed |
| v0.4.4 + v0.5.0 + v0.5.1 splice regression | 48 passed, 13 skipped, 0 failed |
| v0.5.1 focused tests | 11 passed, 0 failed |
| Python compileall | PASS |
| Bash syntax for v0.5.1 runners/installers | PASS |
| Synthetic three-chain Python end-to-end build | PASS |
| Synthetic `run_splice_provenance_v051.sh` smoke | PASS |
| EasyQuant wrapper with contract-compatible stub | PASS |
| k4neo wrapper with contract-compatible stub and license gate | PASS |
| pVACsplice wrapper with contract-compatible stub | PASS |
| Wheel build with local non-isolated build environment | PASS |
| Wheel installation in a fresh virtual environment | PASS |
| `neoag.__version__` and schema version | 0.5.1 |
| `neoag-splice-layer --help` from installed wheel | PASS |
| Incremental patch applied to a clean v0.5.0 source tree | PASS |
| Patch-installed focused regression suite | PASS |
| Full-source ZIP/TAR extraction and checksum verification | PASS |
| Re-extracted source focused regression suite | PASS |

Skipped tests are pre-existing opt-in tests for external tools, large references, licensed resources, workflow engines, and site-specific installations. They are not counted as passes.

## Provenance invariants covered by tests

- no same-gene or nearest-locus junction-read fallback;
- canonical strand-aware junction IDs;
- canonical variant IDs;
- exact splice2neo variant-junction linking;
- EasyQuant import only by project-generated query ID;
- k4neo import only by project-generated cts_id/query ID;
- pVACsplice import only by exact variant and strand-aware junction/event;
- dual-generator peptide consensus only for the same event and exact peptide;
- peptide-level moPepGen output does not falsely establish full-ORF consensus;
- k4neo-negative evidence remains distinct from coverage-aware locus negativity;
- unknown or ambiguous external IDs are retained as conflicts and do not contribute positive evidence;
- foreign-key and manifest hash validation remains intact.

## External-tool validation boundary

The release environment did not execute production installations of moPepGen, splice2neo, EasyQuant, pVACsplice, or k4neo against real patient data. External adapters and wrappers were tested with synthetic, contract-compatible fixtures and stub executables. Therefore this release establishes software contracts and provenance behavior, not biological sensitivity/specificity or clinical validity.

Production deployment must independently verify:

- exact external-tool version/commit or immutable container digest;
- parser compatibility with the local version’s real outputs;
- reference FASTA/GTF/build consistency;
- sample identity and RNA/DNA pairing;
- HLA caller and pVACtools version;
- normal-reference composition and coverage;
- k4neo database/index provenance and licensing;
- licensed prediction software and model assets;
- performance on public truth sets and local validation samples.

## License boundary

k4neo is not bundled. The production CLI and runner require explicit `--k4neo-license-accepted` before importing or executing k4neo results. This flag records operator acknowledgement only; it does not grant a commercial license.

## Scientific boundary

A resolved ORF remains a computational translation hypothesis unless supported by long-read, proteomic, ligandomic, or functional evidence. HLA prediction does not prove endogenous presentation, and endogenous presentation does not prove T-cell recognition. The final evidence tiers are research prioritization outputs, not clinical treatment decisions.
