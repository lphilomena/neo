# Release packages (2026-06-08)

## Lite (recommended for upload)

```
../neoag_event_pipeline_v03_rc_lite_20260608.tar.gz   (27 MB)
```

SHA256: `dist/neoag_event_pipeline_v03_rc_lite_20260608.tar.gz.sha256`  
Details: `dist/README_LITE.txt`

## Full (offline / all tools included)

```
../neoag_event_pipeline_v03_rc_full_with_git_20260608.tar.gz   (~5.3 GB)
```

SHA256: `dist/neoag_event_pipeline_v03_rc_full_with_git_20260608.tar.gz.sha256`  
Git commit inside full archive: `65f67ec` (branch `main`)

## Extract (lite)

```bash
tar -xzf neoag_event_pipeline_v03_rc_lite_20260608.tar.gz
cd neoag_event_pipeline_v03_rc
pip install -e ".[test]"
source conf/tools.env.sh
bash scripts/install_immunogenicity_tools.sh
neoag-v03 check-tools
```

See `RELEASE.md` in project root.
