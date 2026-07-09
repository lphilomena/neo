from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from .schemas import CaseState


def new_case_id(prefix: str = "neoag_case") -> str:
    return f"{prefix}_{time.strftime('%Y%m%d_%H%M%S')}"


def write_case_state(path: str | Path, state: CaseState) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(state.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def append_audit(path: str | Path, event: dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(event, ensure_ascii=False, sort_keys=True) + "\n")
