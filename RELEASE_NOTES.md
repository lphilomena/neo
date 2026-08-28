# Open-Neo Main Release 2026-08-28

This release merges `na0707_upload_release` into `main` and delivers a production-ready three-skill workflow:

- Improves `open-neo-install-check`, `open-neo-run`, and `open-neo-review` for installation readiness, complete analysis, and report generation.
- Supports resumable public-asset synchronization from Hugging Face, configurable mirrors, and portable deployment checks.
- Hardens SV, fusion, and splice evidence provenance, QC, peptide-origin validation, and candidate deduplication.
- Expands integrated ranking with RNA expression/VAF, purity and CNV consensus, HLA-LOH, CCF, MT/WT specificity, and quantitative presentation evidence.
- Adds controlled CPU/memory scheduling, interruption-safe resume behavior, and clearer patient and technical reports.

Validation before the merge: **745 passed, 115 skipped**.

> Licensed predictors such as NetMHCpan, NetMHCstabpan, and NetChop still require legally obtained distributions or license files.
