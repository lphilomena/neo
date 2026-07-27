# Open-Neo macro Skills validation report

## Validation

```bash
PYTHONPATH=src python -m compileall -q src
PYTHONPATH=src python -m pytest -q
PYTHONPATH=src python -m neoag.skill_taxonomy.cli validate --root . --outdir /tmp/openneo_skill_validate
PYTHONPATH=src python -m neoag.open_neo.cli --help
```

Result:

```text
350 passed, 104 skipped
Skill directory validation: PASS
Open-Neo CLI: PASS
```

Focused macro-Skill and Gateway tests: `22 passed`.

## External boundary

Skipped tests include site-specific external tools, licensed assets, large
references, and network-dependent workflows. Doctor and deployment-tier
outputs report these capabilities as missing or unassessed; they are not
interpreted as biological negative evidence.
