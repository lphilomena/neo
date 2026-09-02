from __future__ import annotations

import hashlib
import tomllib
from pathlib import Path
from typing import Any


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_cohort_rule_contract(path: str | Path) -> dict[str, Any]:
    contract_path = Path(path).expanduser().resolve()
    payload = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    metadata = dict(payload.get("metadata") or {})
    compatibility = dict(payload.get("compatibility") or {})
    knowledge = dict(payload.get("knowledge") or {})
    if not metadata.get("id") or not metadata.get("version"):
        raise ValueError(f"cohort rule contract is missing metadata.id/version: {contract_path}")
    if not compatibility.get("ranking_profile") or not compatibility.get("evidence_consensus_rules"):
        raise ValueError(f"cohort rule contract is missing the profile/rules pair: {contract_path}")
    project_root = contract_path.parents[2]
    profile = (project_root / str(compatibility["ranking_profile"])).resolve()
    rules = (project_root / str(compatibility["evidence_consensus_rules"])).resolve()
    if not profile.is_file() or not rules.is_file():
        raise ValueError(f"cohort rule contract references missing files: profile={profile}; rules={rules}")
    knowledge_file = str(knowledge.get("disease_knowledge") or "").strip()
    knowledge_path = (project_root / knowledge_file).resolve() if knowledge_file else None
    if knowledge_file and not knowledge_path.is_file():
        raise ValueError(
            f"cohort rule contract references missing disease knowledge file: {knowledge_path}"
        )
    return {
        "path": str(contract_path),
        "id": str(metadata["id"]),
        "version": str(metadata["version"]),
        "status": str(metadata.get("status") or "RESEARCH_ONLY"),
        "disease": str(metadata.get("disease") or ""),
        "ranking_profile": str(profile),
        "ranking_profile_version": str(compatibility.get("ranking_profile_version") or ""),
        "evidence_consensus_rules": str(rules),
        "evidence_consensus_rules_version": str(compatibility.get("evidence_consensus_rules_version") or ""),
        "report_contract_version": str(compatibility.get("report_contract_version") or ""),
        "release_audit_policy": str(compatibility.get("release_audit_policy") or ""),
        "contract_sha256": sha256_file(contract_path),
        "disease_knowledge_file": str(knowledge_path) if knowledge_path else "",
        "ranking_profile_sha256": sha256_file(profile),
        "evidence_consensus_rules_sha256": sha256_file(rules),
    }


def validate_cohort_rule_pair(
    contract: dict[str, Any],
    *,
    ranking_profile: str | Path,
    evidence_consensus_rules: str | Path,
) -> list[str]:
    mismatches: list[str] = []
    actual_profile = Path(ranking_profile).expanduser().resolve()
    actual_rules = Path(evidence_consensus_rules).expanduser().resolve()
    if actual_profile != Path(str(contract["ranking_profile"])):
        mismatches.append(f"ranking_profile={actual_profile}; expected={contract['ranking_profile']}")
    if actual_rules != Path(str(contract["evidence_consensus_rules"])):
        mismatches.append(f"evidence_consensus_rules={actual_rules}; expected={contract['evidence_consensus_rules']}")
    if actual_profile.is_file() and sha256_file(actual_profile) != contract["ranking_profile_sha256"]:
        mismatches.append("ranking_profile_sha256 changed after contract resolution")
    if actual_rules.is_file() and sha256_file(actual_rules) != contract["evidence_consensus_rules_sha256"]:
        mismatches.append("evidence_consensus_rules_sha256 changed after contract resolution")
    return mismatches


def discover_matching_cohort_contract(
    ranking_profile: str | Path,
    evidence_consensus_rules: str | Path,
    *,
    explicit: str | Path | None = None,
) -> dict[str, Any] | None:
    profile = Path(ranking_profile).expanduser().resolve()
    rules = Path(evidence_consensus_rules).expanduser().resolve()
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    for source in (profile, rules):
        for parent in source.parents:
            directory = parent / "configs" / "cohorts"
            if directory.is_dir():
                candidates.extend(
                    path for path in sorted(directory.glob("*.toml"))
                    if not path.name.startswith(".")
                )
                break
    for candidate in dict.fromkeys(path.resolve() for path in candidates if path.is_file()):
        try:
            contract = load_cohort_rule_contract(candidate)
        except (OSError, UnicodeDecodeError, ValueError, tomllib.TOMLDecodeError):
            continue
        if not validate_cohort_rule_pair(
            contract,
            ranking_profile=profile,
            evidence_consensus_rules=rules,
        ):
            return contract
    return None
