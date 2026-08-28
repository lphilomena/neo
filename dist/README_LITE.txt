# Lite release package (v0.3.0rc2)

Upload file:

```
dist/neoag_event_pipeline_lite_<DATE>.tar.gz
```

SHA256: see matching `.sha256` file alongside the tarball.

## Extract

```bash
tar -xzf neoag_event_pipeline_lite_<DATE>.tar.gz
cd neoag_event_pipeline
pip install -e ".[test]"
source conf/tools.env.sh
bash scripts/install_immunogenicity_tools.sh   # PRIME + BigMHC (~5 GB)
bash scripts/install_deepimmuno.sh             # optional DeepImmuno-CNN
neoag check-tools
```

## Included

- Pipeline source, tests, configs, docs, install scripts
- HCC1395 demo inputs (`data/examples/HCC1395/`, ~27 MB)
- Verified outputs (reports/metrics/scoring; tool caches excluded):
  - `results/HCC1395/` — real PRIME + BigMHC_IM rescored peptides
  - `results/bench_phase1_complete/`, `bench_phase2/`, `bench_phase2_deepimmuno/` — CEDAR benchmark
- CEDAR benchmark labels (`data/improve/`)

## Excluded (install separately)

| Component | Restore |
|-----------|---------|
| BigMHC models (~4.8 GB) | `bash scripts/install_immunogenicity_tools.sh` |
| PRIME + MixMHCpred | same script |
| DeepImmuno-CNN (~50 MB + TensorFlow) | `bash scripts/install_deepimmuno.sh` |
| NetMHCpan binary | `bash scripts/install_netmhcpan.sh` |
| LOHHLA / FACETS / pctGCdata | `bash scripts/install_lohhla_facets_nextflow.sh` |
| `.git` history | use full tarball if needed |
| `results/*/tools/` caches | re-run benchmark or copy from full install |

## Quick verify

```bash
pytest -q

neoag benchmark-improve --dataset cedar --outdir results/bench_phase2 \
  --skip-netmhcpan --skip-mhcflurry --skip-stabpan \
  --reuse-tools-dir results/benchmark_cedar_v2/tools
```

See `RELEASE.md` in project root.
