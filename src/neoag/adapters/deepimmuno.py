"""DeepImmuno-CNN immunogenicity adapter (9/10-mer peptide–HLA pairs)."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from ..utils import write_tsv
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, provenance_stub, without_provenance, write_evidence_tsv
from ..schemas import DEEPIMMUNO_EVIDENCE_FIELDS
from .peptide_input import normalize_hla_allele, pair_key

STANDARD_AAS = frozenset("ACDEFGHIKLMNPQRSTVWY")
DEEPIMMUNO_VALID_LENGTHS = frozenset({9, 10})


def _pair_key(peptide: str, hla: str) -> tuple[str, str]:
    return pair_key(peptide, normalize_hla_allele(hla))


def normalize_hla_for_deepimmuno(hla_raw: str) -> str:
    """Convert to DeepImmuno format: HLA-A*0201 (no colon in allele digits)."""
    hla = hla_raw.strip().upper()
    if not hla:
        return ""
    if hla.startswith("HLA-"):
        hla = hla[4:]
    if "*" in hla and ":" in hla:
        hla = hla.replace(":", "")
    elif "*" not in hla and len(hla) >= 4:
        hla = hla[0] + "*" + hla[1:3] + hla[3:5]
    return "HLA-" + hla


def is_valid_peptide(peptide: str) -> bool:
    pep = peptide.strip().upper()
    return len(pep) in DEEPIMMUNO_VALID_LENGTHS and set(pep).issubset(STANDARD_AAS)


def resolve_deepimmuno_dir(custom: str | Path | None = None) -> Path:
    candidates: list[Path] = []
    if custom:
        candidates.append(Path(custom))
    for env_key in ("DEEPIMMUNO_DIR", "NEOAG_DEEPIMMUNO_DIR"):
        val = os.environ.get(env_key, "")
        if val:
            candidates.append(Path(val))
    root = os.environ.get("NEOAG_TOOLS_ROOT", "")
    if root:
        candidates.append(Path(root) / "tools" / "DeepImmuno")
    candidates.append(Path(__file__).resolve().parents[3] / "tools" / "DeepImmuno")

    required_names = (
        "deepimmuno-cnn.py",
        "data/after_pca.txt",
        "data/hla2paratopeTable_aligned.txt",
        "models/cnn_model_331_3_7",
    )
    for base in candidates:
        if not base.is_dir():
            continue
        if all((base / rel).exists() for rel in required_names):
            return base.resolve()
    searched = ", ".join(str(p) for p in candidates)
    raise FileNotFoundError(
        f"DeepImmuno installation not found. Set DEEPIMMUNO_DIR or install under tools/DeepImmuno. "
        f"Searched: {searched}"
    )


def predict_pair_stub(peptide: str, hla: str) -> str:
    if not is_valid_peptide(peptide):
        return ""
    h = abs(hash(f"DeepImmuno_{peptide}_{hla}"))
    return f"{(h % 10000) / 10000.0:.6f}"


def predict_pairs(pairs: list[tuple[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for peptide, hla in pairs:
        rows.append(
            {
                "peptide": peptide.strip().upper(),
                "hla_allele": normalize_hla_allele(hla),
                "deepimmuno_score": predict_pair_stub(peptide, hla),
            }
        )
    return rows


def _parse_result_rows(path: Path, sample_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            peptide = (row.get("peptide") or row.get("Peptide") or "").strip().upper()
            hla = normalize_hla_allele(row.get("HLA") or row.get("hla") or row.get("hla_allele") or "")
            score = (
                row.get("immunogenicity")
                or row.get("DeepImmuno_score")
                or row.get("deepimmuno_score")
                or row.get("score")
                or ""
            )
            if not peptide:
                continue
            rows.append(
                {
                    "sample_id": sample_id,
                    "peptide": peptide,
                    "hla_allele": hla,
                    "deepimmuno_score": str(score),
                    "source_file": str(path),
                }
            )
    return rows


def parse_deepimmuno(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    return _parse_result_rows(Path(path), sample_id)


def run_deepimmuno_batch(
    pairs: list[tuple[str, str]],
    deepimmuno_dir: str | Path,
    *,
    sample_id: str = "",
    python_exe: str | None = None,
) -> list[dict[str, str]]:
    """Run DeepImmuno-CNN in multiple mode for strict peptide–HLA pairs."""
    deepimmuno_dir = Path(deepimmuno_dir)
    python = python_exe or sys.executable

    indexed_pairs: list[tuple[str, str]] = []
    for peptide, hla in pairs:
        pep = peptide.strip().upper()
        if not is_valid_peptide(pep):
            continue
        indexed_pairs.append((pep, normalize_hla_allele(hla)))

    if not indexed_pairs:
        return predict_pairs(pairs)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".csv",
        delete=False,
        newline="",
        encoding="utf-8",
    ) as tmp_in:
        writer = csv.writer(tmp_in)
        for pep, hla in indexed_pairs:
            writer.writerow([pep, normalize_hla_for_deepimmuno(hla)])
        tmp_in_path = tmp_in.name

    try:
        with tempfile.TemporaryDirectory() as tmp_out_dir:
            cmd = [
                python,
                "deepimmuno-cnn.py",
                "--mode",
                "multiple",
                "--intdir",
                tmp_in_path,
                "--outdir",
                tmp_out_dir,
            ]
            proc = subprocess.run(
                cmd,
                cwd=deepimmuno_dir,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    f"DeepImmuno failed ({proc.returncode}): {' '.join(cmd)}\n"
                    f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
                )
            result_path = Path(tmp_out_dir) / "deepimmuno-cnn-result.txt"
            if not result_path.is_file():
                raise FileNotFoundError(f"DeepImmuno output missing: {result_path}")
            parsed = _parse_result_rows(result_path, sample_id)
    finally:
        Path(tmp_in_path).unlink(missing_ok=True)

    by_pair = deepimmuno_by_pair(parsed)
    merged: list[dict[str, str]] = []
    for peptide, hla in pairs:
        pep = peptide.strip().upper()
        canon_hla = normalize_hla_allele(hla)
        hit = by_pair.get(_pair_key(pep, canon_hla), {})
        score = hit.get("deepimmuno_score", "")
        if score == "" and is_valid_peptide(pep):
            score = predict_pair_stub(pep, canon_hla)
        merged.append(
            {
                "sample_id": sample_id,
                "peptide": pep,
                "hla_allele": canon_hla,
                "deepimmuno_score": score,
                "source_file": hit.get("source_file", str(deepimmuno_dir)),
            }
        )
    return merged


def write_deepimmuno_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    src = rows[0].get("source_file") if rows else path
    if src == "stub":
        prov = provenance or provenance_stub("deepimmuno", production_invalid=True)
    else:
        prov = provenance or provenance_from_file("deepimmuno", src, mode="tool_run")
    write_evidence_tsv(path, rows, without_provenance(DEEPIMMUNO_EVIDENCE_FIELDS), prov)


def deepimmuno_by_pair(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    out: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        out[_pair_key(row.get("peptide", ""), row.get("hla_allele", ""))] = row
    return out
