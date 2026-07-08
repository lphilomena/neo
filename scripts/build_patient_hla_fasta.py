#!/usr/bin/env python3
"""Build patient-specific HLA FASTA from Polysolver winners.hla.txt."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_winners(path: Path) -> list[str]:
    alleles: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3 and parts[0].upper().startswith("HLA-"):
            alleles.extend([parts[1].strip(), parts[2].strip()])
        elif len(parts) == 1 and parts[0]:
            alleles.append(parts[0].strip())
    out = sorted({a for a in alleles if a and a.lower() not in {"na", "none", "-"}})
    if not out:
        raise SystemExit(f"No alleles parsed from {path}")
    return out


def _faidx_extract(samtools: str, ref_fasta: Path, allele: str) -> str | None:
    proc = subprocess.run(
        [samtools, "faidx", str(ref_fasta), allele],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    lines = proc.stdout.splitlines()
    if len(lines) < 2:
        return None
    return "\n".join(lines) + "\n"


def extract_allele_record(
    allele: str,
    ref_fasta: Path,
    complete_dir: Path,
    samtools: str,
) -> str:
    record = _faidx_extract(samtools, ref_fasta, allele)
    if record:
        return record

    per_allele = complete_dir / f"{allele}.fasta"
    if per_allele.is_file():
        return per_allele.read_text(encoding="utf-8", errors="replace")

    # Polysolver allele IDs occasionally differ by trailing suffix; try prefix match.
    if complete_dir.is_dir():
        matches = sorted(complete_dir.glob(f"{allele}*.fasta"))
        if matches:
            return matches[0].read_text(encoding="utf-8", errors="replace")

    raise SystemExit(f"Could not extract sequence for allele '{allele}' from {ref_fasta} or {complete_dir}")


def build_patient_fasta(
    winners: Path,
    ref_fasta: Path,
    complete_dir: Path,
    out_fasta: Path,
    samtools: str,
) -> list[str]:
    alleles = parse_winners(winners)
    records: list[str] = []
    for allele in alleles:
        records.append(
            extract_allele_record(allele, ref_fasta, complete_dir, samtools).rstrip() + "\n"
        )

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_fasta.write_text("".join(records), encoding="utf-8")
    return alleles


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--winners", required=True, type=Path, help="Polysolver winners.hla.txt")
    parser.add_argument(
        "--ref-fasta",
        type=Path,
        default=Path(os.environ.get("POLYSOLVER_HOME", "tools/polysolver")) / "data" / "abc_complete.fasta",
        help="Polysolver HLA reference FASTA (default: $POLYSOLVER_HOME/data/abc_complete.fasta or tools/polysolver/data/abc_complete.fasta)",
    )
    parser.add_argument(
        "--complete-dir",
        type=Path,
        default=Path(os.environ.get("POLYSOLVER_HOME", "tools/polysolver")) / "data" / "complete",
        help="Per-allele FASTA fallback directory (default: $POLYSOLVER_HOME/data/complete or tools/polysolver/data/complete)",
    )
    parser.add_argument("--out-fasta", required=True, type=Path, help="Patient HLA FASTA output")
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    parser.add_argument(
        "--novoindex",
        default="",
        help="novoindex executable; if set, build companion .nix next to out-fasta",
    )
    args = parser.parse_args(argv)

    if not args.ref_fasta.is_file():
        raise SystemExit(f"Reference FASTA missing: {args.ref_fasta}")
    if not args.winners.is_file():
        raise SystemExit(f"winners.hla.txt missing: {args.winners}")

    alleles = build_patient_fasta(
        args.winners,
        args.ref_fasta,
        args.complete_dir,
        args.out_fasta,
        args.samtools,
    )
    print(f"Wrote {args.out_fasta} ({len(alleles)} alleles)")

    if args.novoindex:
        nix = args.out_fasta.with_suffix(".nix")
        subprocess.run([args.novoindex, str(args.out_fasta)], check=True)
        if nix.is_file() or args.out_fasta.with_name(args.out_fasta.name + ".nix").is_file():
            print(f"Built Novoalign index for {args.out_fasta}")
        else:
            # novoindex names index from fasta basename
            print(f"Ran novoindex on {args.out_fasta}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
