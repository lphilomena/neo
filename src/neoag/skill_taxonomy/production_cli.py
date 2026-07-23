from __future__ import annotations

import contextlib
import io
import shlex
from typing import Sequence


def invoke_production_cli(command: str, arguments: Sequence[str]) -> dict[str, object]:
    """Invoke a registered production CLI in-process without duplicating its logic."""
    tokens = shlex.split(command)
    if tokens[:1] == ["neoag"]:
        tokens = tokens[1:]
    if tokens != ["evidence-rank"]:
        raise ValueError(f"Unsupported production CLI command: {command}")

    from neoag.cli import main as production_main

    stdout = io.StringIO()
    stderr = io.StringIO()
    exit_code = 0
    try:
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            result = production_main([*tokens, *map(str, arguments)])
            if isinstance(result, int):
                exit_code = result
    except SystemExit as exc:
        exit_code = int(exc.code or 0) if isinstance(exc.code, int) else 1
    if exit_code:
        raise RuntimeError(
            f"Production CLI failed ({exit_code}): {command}\n"
            f"stdout:\n{stdout.getvalue()}\nstderr:\n{stderr.getvalue()}"
        )
    return {
        "command": command,
        "arguments": list(map(str, arguments)),
        "exit_code": exit_code,
        "stdout": stdout.getvalue(),
        "stderr": stderr.getvalue(),
    }
