#!/usr/bin/env python3
"""Patch PSAURON/TD2 for CPU fallback and explicit subprocess validation."""

from __future__ import annotations

import argparse
from pathlib import Path


def patch_file(path: Path, replacements: list[tuple[str, str]]) -> bool:
    original = path.read_text()
    updated = original
    for old, new in replacements:
        updated = updated.replace(old, new)
    if updated == original:
        return False
    backup = path.with_suffix(path.suffix + ".neoag.bak")
    if not backup.exists():
        backup.write_text(original)
    path.write_text(updated)
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--env", required=True, help="Conda environment root")
    args = parser.parse_args()

    site = Path(args.env) / "lib/python3.11/site-packages"
    psauron = site / "psauron/psauron.py"
    td2 = site / "TD2/Predict.py"
    if not psauron.is_file() or not td2.is_file():
        raise SystemExit(f"PSAURON/TD2 package files missing under {site}")

    psauron_changed = patch_file(
        psauron,
        [("torch.cuda.device_count() > 0", "torch.cuda.is_available()")],
    )

    marker = "        result_psauron = subprocess.run(command_psauron, capture_output=True, text=True, env=env)\n"
    validation = marker + "\n    if result_psauron.returncode != 0:\n        if result_psauron.stdout:\n            print(result_psauron.stdout, file=sys.stderr)\n        if result_psauron.stderr:\n            print(result_psauron.stderr, file=sys.stderr)\n        raise RuntimeError(f\"PSAURON failed with exit code {result_psauron.returncode}: {' '.join(command_psauron)}\")\n    if not os.path.isfile(p_score) or os.path.getsize(p_score) == 0:\n        raise RuntimeError(f\"PSAURON did not create a non-empty score file: {p_score}\")\n"
    td2_changed = False
    td2_text = td2.read_text()
    if "PSAURON did not create a non-empty score file" not in td2_text:
        if marker not in td2_text:
            raise SystemExit("TD2 Predict.py layout is unsupported; validation marker not found")
        td2_changed = patch_file(td2, [(marker, validation)])

    final = psauron.read_text()
    if "torch.cuda.device_count() > 0" in final:
        raise SystemExit("PSAURON patch verification failed")
    print(f"psauron_cpu_fallback={'patched' if psauron_changed else 'already_ok'}")
    print(f"td2_subprocess_validation={'patched' if td2_changed else 'already_ok'}")


if __name__ == "__main__":
    main()
