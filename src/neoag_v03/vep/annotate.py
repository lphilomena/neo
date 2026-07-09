from __future__ import annotations

import os
import subprocess
import tomllib
from pathlib import Path
from typing import Any


def load_vep_annotate_config(path: str | Path) -> dict[str, Any]:
    return tomllib.loads(Path(path).read_text(encoding="utf-8"))


def resolve_vep_annotate_from_config(cfg: dict[str, Any], *, root: str | Path | None = None) -> dict[str, Any]:
    base = Path(root or ".").resolve()
    vep = cfg.get("vep") or cfg.get("inputs") or cfg

    def path_value(*keys: str) -> str | None:
        for key in keys:
            val = vep.get(key)
            if val:
                p = Path(str(val))
                return str(p if p.is_absolute() else (base / p))
        return None

    return {
        "input_vcf": path_value("input_vcf", "vcf", "tumor_vcf"),
        "output_vcf": path_value("output_vcf", "out_vcf", "annotated_vcf"),
        "reference_fasta": path_value("reference_fasta", "fasta"),
        "workdir": path_value("workdir", "work_dir"),
        "cache_dir": path_value("cache_dir"),
        "plugins_dir": path_value("plugins_dir", "plugin_dir"),
        "online": vep.get("online"),
        "fork": int(vep.get("fork", 4) or 4),
        "pick": bool(vep.get("pick", True)),
        "expression_custom": path_value("expression_custom"),
        "index_vcf": bool(vep.get("index_vcf", True)),
        "vep_bin": vep.get("vep_bin") or os.environ.get("NEOAG_VEP_BIN") or "vep",
        "cache_version": vep.get("cache_version") or os.environ.get("NEOAG_VEP_CACHE_VERSION"),
    }


def build_vep_pvacseq_command(
    *,
    input_vcf: str | Path,
    output_vcf: str | Path,
    reference_fasta: str | Path,
    cache_dir: str | Path | None = None,
    plugins_dir: str | Path | None = None,
    online: bool | None = None,
    fork: int = 4,
    pick: bool = True,
    expression_custom: str | Path | None = None,
    vep_bin: str | Path | None = None,
    cache_version: str | int | None = None,
) -> list[str]:
    cmd = [
        str(vep_bin or os.environ.get("NEOAG_VEP_BIN") or "vep"),
        "--input_file", str(input_vcf),
        "--output_file", str(output_vcf),
        "--format", "vcf",
        "--vcf",
        "--force_overwrite",
        "--symbol",
        "--terms", "SO",
        "--tsl",
        "--biotype",
        "--hgvs",
        "--fasta", str(reference_fasta),
        "--fork", str(fork),
    ]
    if pick:
        cmd.append("--pick")
    if online is True or (online is None and os.environ.get("NEOAG_VEP_ONLINE", "").lower() in {"1", "true", "yes"}):
        # VEP online mode should use Ensembl public databases explicitly.
        cmd.extend(["--database", "--species", "homo_sapiens"])
    else:
        cmd.extend(["--cache", "--offline"])
        if cache_dir:
            cmd.extend(["--dir_cache", str(cache_dir)])
        ver = cache_version or os.environ.get("NEOAG_VEP_CACHE_VERSION")
        if ver:
            cmd.extend(["--cache_version", str(ver)])
    if plugins_dir:
        cmd.extend(["--dir_plugins", str(plugins_dir)])
        cmd.extend(["--plugin", "Wildtype"])
        cmd.extend(["--plugin", "Frameshift"])
    if expression_custom:
        cmd.extend(["--custom", f"{expression_custom},Expression,tsv,exact,0,TPM"])
    return cmd


def run_vep_pvacseq_annotate(
    *,
    input_vcf: str | Path,
    output_vcf: str | Path,
    reference_fasta: str | Path,
    workdir: str | Path | None = None,
    cache_dir: str | Path | None = None,
    plugins_dir: str | Path | None = None,
    online: bool | None = None,
    fork: int = 4,
    pick: bool = True,
    expression_custom: str | Path | None = None,
    index_vcf: bool = True,
    vep_bin: str | Path | None = None,
    cache_version: str | int | None = None,
) -> dict[str, str]:
    out = Path(output_vcf)
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = build_vep_pvacseq_command(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        reference_fasta=reference_fasta,
        cache_dir=cache_dir,
        plugins_dir=plugins_dir,
        online=online,
        fork=fork,
        pick=pick,
        expression_custom=expression_custom,
        vep_bin=vep_bin,
        cache_version=cache_version,
    )
    subprocess.run(cmd, cwd=workdir or None, check=True)
    if index_vcf and str(out).endswith(".gz"):
        subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=False)
    return {"annotated_vcf": str(out), "command": " ".join(cmd)}
