"""Load and check conf/tools.toml — the single, TOML-sectioned source of truth for
which external tool each entry_mode needs, and how to verify it's installed.

Replaces ad-hoc reading of dozens of flat NEOAG_* environment variables scattered
across conf/tools.env.sh: each tool gets one section here, with a declared check
type (bin / dir / env) and the list of entry_mode values that need it.
"""
from __future__ import annotations

import os
import re
import shutil
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TOOLS_TOML = ROOT / "conf" / "tools.toml"

_VAR_RE = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(:-([^}]*))?\}")


def _expand(value: str) -> str:
    def repl(m: re.Match) -> str:
        name, _, default = m.groups()
        return os.environ.get(name, default or "")

    return _VAR_RE.sub(repl, value)


def load_tools_toml(path: str | Path | None = None) -> dict:
    """Load conf/tools.toml, expanding ${VAR}/${VAR:-default} in string values."""
    p = Path(path or os.environ.get("NEOAG_TOOLS_TOML") or DEFAULT_TOOLS_TOML)
    if not p.is_file():
        return {}
    with p.open("rb") as fh:
        data = tomllib.load(fh)
    for section in data.values():
        if not isinstance(section, dict):
            continue
        for key, val in list(section.items()):
            if isinstance(val, str):
                section[key] = _expand(val)
    return data


def check_tool(name: str, section: dict) -> tuple[bool, str]:
    """Check a single tool section. Returns (available, message)."""
    check = section.get("check", "bin")
    optional = bool(section.get("optional", False))
    tag = "optional" if optional else "required"

    if check == "bin":
        bin_value = (section.get("bin") or "").strip()
        if not bin_value:
            return False, f"{name} ({tag}): MISSING (no bin configured)"
        found = shutil.which(bin_value) or (Path(bin_value).is_file() and bin_value)
        if found:
            return True, f"{name} ({tag}): OK ({found})"
        return False, f"{name} ({tag}): MISSING (bin '{bin_value}' not found on PATH)"

    if check == "dir":
        dir_value = (section.get("dir") or "").strip()
        if dir_value and Path(dir_value).is_dir():
            return True, f"{name} ({tag}): OK ({dir_value})"
        return False, f"{name} ({tag}): MISSING (dir '{dir_value or '<unset>'}' not found)"

    if check == "env":
        env_var = section.get("env_var", "")
        val = os.environ.get(env_var, "") if env_var else ""
        if val:
            return True, f"{name} ({tag}): OK ({env_var}={val})"
        return False, f"{name} ({tag}): MISSING ({env_var or '<unset>'} not set)"

    return False, f"{name} ({tag}): MISSING (unknown check type {check!r})"


def tools_for_entry(entry_mode: str, data: dict | None = None) -> dict[str, dict]:
    data = data if data is not None else load_tools_toml()
    return {
        name: section
        for name, section in data.items()
        if isinstance(section, dict) and entry_mode in (section.get("entries") or [])
    }


def check_entry_tools(entry_mode: str, data: dict | None = None) -> dict[str, tuple[bool, str]]:
    """Check every tool section relevant to entry_mode. Returns {name: (ok, message)}."""
    sections = tools_for_entry(entry_mode, data)
    return {name: check_tool(name, section) for name, section in sections.items()}
