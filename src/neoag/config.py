from __future__ import annotations
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROFILE_DIR = ROOT / "profiles"

def load_profile(profile: str | Path) -> dict:
    p = Path(profile)
    if not p.exists():
        p = PROFILE_DIR / (str(profile) if str(profile).endswith(".toml") else f"{profile}.toml")
    if not p.exists():
        raise FileNotFoundError(f"Profile not found: {profile}")
    with p.open("rb") as fh:
        data = tomllib.load(fh)
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
