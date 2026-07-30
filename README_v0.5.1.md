# NeoAg Event Pipeline v0.5.1

This release adds a dual-generator RNA branch, an exact mutation-caused-splicing branch, and a separately audited normal-background branch to the formal v0.5.0 Splice Provenance Layer.

Start here:

```bash
python -m pip install -e .
neoag-splice-layer --help
bash scripts/run_splice_provenance_v051.sh --help
PYTHONPATH=src pytest -q tests/test_splice_v051_three_chains.py
```

Read `docs/V051_THREE_EVIDENCE_CHAINS.md` before production deployment. External tools and databases are not bundled. Their exact versions, licenses, references and local assets remain site-controlled.
