from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VcfPreflight:
    path: str
    exists: bool
    has_csq_header: bool
    has_gt_format: bool
    variant_rows: int
    sample_columns: int

    @property
    def pvacseq_ready(self) -> bool:
        return self.exists and self.has_csq_header and self.has_gt_format

    def require_pvacseq_ready(self) -> None:
        missing: list[str] = []
        if not self.exists:
            missing.append("file_exists")
        if not self.has_csq_header:
            missing.append("VEP CSQ header")
        if not self.has_gt_format:
            missing.append("sample GT format")
        if missing:
            raise ValueError(
                f"VCF preflight failed for {self.path}: missing {', '.join(missing)}. "
                "Run neoag vep-annotate before pVACseq enrich, or set inputs.pvacseq_enrich = false for fixture/stub runs."
            )


def vcf_preflight(path: str | Path | None) -> VcfPreflight:
    if path is None:
        return VcfPreflight("", False, False, False, 0, 0)
    p = Path(path)
    if not p.is_file():
        return VcfPreflight(str(p), False, False, False, 0, 0)
    opener = gzip.open if str(p).endswith(".gz") else open
    has_csq = False
    has_gt = False
    variants = 0
    samples = 0
    with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                has_csq = True
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = max(0, len(parts) - 9)
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 8:
                variants += 1
            if len(parts) >= 9 and "GT" in parts[8].split(":"):
                has_gt = True
    return VcfPreflight(str(p), True, has_csq, has_gt, variants, samples)
