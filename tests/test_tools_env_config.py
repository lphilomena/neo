from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_tools_env_respects_external_tools_root(tmp_path):
    external_root = tmp_path / "tool_bundle"
    external_root.mkdir()
    tools_env = ROOT / "conf" / "tools.env.sh"
    cmd = "\n".join([
        f"export NEOAG_TOOLS_ROOT='{external_root}'",
        f"source '{tools_env}' >/dev/null 2>&1",
        "printf '%s\\n%s\\n' \"$NEOAG_PROJECT_ROOT\" \"$NEOAG_TOOLS_ROOT\"",
    ])
    out = subprocess.check_output(["bash", "-lc", cmd], text=True).splitlines()
    assert out[0] == str(ROOT)
    assert out[1] == str(external_root)
