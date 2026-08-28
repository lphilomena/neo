# Deployment layout

Use a four-way separation on target machines:

```text
/opt/neoag/
  code/
    neoag_event_pipeline/
  conf/
    tools_manifest.yaml
    reference_manifest.yaml
    profiles/
  envs/
    conda/
    containers/
      docker/
      apptainer/
      image_manifest.yaml
  refs/
    GRCh38/
    hla/
    mhc/
    cnv/
    fusion/
    safety/
    validated/
    testdata/
  runs/
    SAMPLE001/
      run_YYYYMMDD_NNN/
        sample_manifest.yaml
        run_manifest.json
        outputs/
        logs/
        audit_log.jsonl
        provenance.json
```

Rules:

- `code/` contains source, small fixtures, schema, docs and tests only.
- `envs/` contains conda environment definitions, container manifests and local
  image locations, not patient data.
- `refs/` is central, versioned and read-only during pipeline execution.
- `runs/` contains sample/run-specific manifests, logs, outputs and provenance.
- New derived references should be staged under `refs/_staging/` and promoted
  only after review.
