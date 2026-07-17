from __future__ import annotations
import tomllib
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
PROFILE_DIR = ROOT / "profiles"
DEFAULT_PROFILE_NAME = "default"


def _resolve_profile_path(profile: str | Path) -> Path:
    p = Path(profile)
    if not p.exists():
        p = PROFILE_DIR / (str(profile) if str(profile).endswith(".toml") else f"{profile}.toml")
    if not p.exists():
        raise FileNotFoundError(f"Profile not found: {profile}")
    return p


def _read_toml(p: Path) -> dict[str, Any]:
    with p.open("rb") as fh:
        return tomllib.load(fh)


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge ``override`` on top of ``base``.

    - dict values are merged key-by-key (so e.g. supplying only
      ``[l3_weights] hla_binding = 0.3`` in a profile no longer silently
      drops every other l3_weights entry -- the rest are inherited from the
      base profile instead of falling back to hardcoded Python defaults).
    - Any other value type (str/int/float/bool/list) in ``override``
      replaces the base value outright.
    """
    merged: dict[str, Any] = dict(base)
    for key, value in override.items():
        base_value = merged.get(key)
        if isinstance(base_value, Mapping) and isinstance(value, Mapping):
            merged[key] = _deep_merge(base_value, value)
        else:
            merged[key] = value
    return merged


def load_profile(profile: str | Path) -> dict:
    """Load a profile TOML, merged on top of ``profiles/default.toml``.

    ``default.toml`` is always loaded first as the base configuration; the
    requested profile's fields are then deep-merged on top of it, so a
    profile only needs to declare the settings it wants to *override* --
    anything it omits (e.g. ``leukemia.toml`` not declaring ``[l3_weights]``
    at all) is inherited from ``default.toml`` explicitly, rather than
    silently depending on hardcoded Python-side fallback constants matching
    default.toml's numbers by coincidence (scoring audit fix #4).

    Requesting the default profile itself (by name "default", or by path)
    is a no-op merge (base merged with itself) and returns the same result
    as before this change.
    """
    p = _resolve_profile_path(profile)
    data = _read_toml(p)

    default_path = PROFILE_DIR / f"{DEFAULT_PROFILE_NAME}.toml"
    if p.resolve() != default_path.resolve() and default_path.exists():
        base = _read_toml(default_path)
        data = _deep_merge(base, data)

    data["_profile_name"] = p.stem
    data["_profile_path"] = str(p)
    return data

def source_priority(profile: dict, event_type: str, mutation_source: str | None = None) -> float:
    src = profile.get("source_priority", {})
    ms = (mutation_source or "").strip()
    if ms and ms in src:
        return float(src[ms])
    if ms and ms.lower() in {k.lower(): k for k in src}:
        key = next(k for k in src if k.lower() == ms.lower())
        return float(src[key])
    return float(src.get(event_type, src.get(event_type.lower(), src.get("Other", 1.0))))
