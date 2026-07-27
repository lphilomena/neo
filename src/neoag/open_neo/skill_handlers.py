from __future__ import annotations

from typing import Any

from .install_check import run_install_check
from .run import run_open_neo
from .review import run_review


def run_open_neo_install_check(args: dict[str, Any]) -> dict[str, Any]:
    return run_install_check(args)


def run_open_neo_run(args: dict[str, Any]) -> dict[str, Any]:
    return run_open_neo(args)


def run_open_neo_review(args: dict[str, Any]) -> dict[str, Any]:
    return run_review(args)
