"""Standalone peptide–HLA prediction workflow (neoantigen2-style input adapter)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from .adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from .adapters.netmhcstabpan import parse_netmhcstabpan, write_netmhcstabpan_evidence
from .adapters.peptide_input import PeptideInputSummary, convert_peptide_input
from .config import load_profile
from .immunogenicity_composite import apply_immunogenicity_evidence, run_immunogenicity_predictors
from .presentation import build_presentation_evidence
from .evidence_provenance import ProvenanceRegistry
from .schemas import PRESENTATION_FIELDS
from .tools.registry import RunContext
from .tools.runner import run_bigmhc_im, run_mhcflurry, run_netmhcpan, run_netmhcstabpan
from .utils import read_tsv, write_json, write_tsv


def run_peptide_predict(
    input_path: str | Path,
    outdir: str | Path,
    *,
    sample_id: str = "SAMPLE001",
    profile_name: str = "default",
    stub: bool = False,
    skip_netmhcpan: bool = False,
    skip_mhcflurry: bool = False,
    skip_prime: bool = False,
    skip_bigmhc_im: bool = False,
    skip_deepimmuno: bool = False,
    skip_stabpan: bool = False,
    tool_executables: dict[str, str] | None = None,
) -> dict[str, Any]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    profile = load_profile(profile_name)

    summary = convert_peptide_input(input_path, outdir, sample_id=sample_id)
    raw_peptides = Path(summary.raw_peptides_tsv)
    pres_dir = outdir / "presentation"
    pres_dir.mkdir(exist_ok=True)

    ctx = RunContext(
        sample_id=sample_id,
        outdir=outdir,
        stub=stub,
        raw_peptides=raw_peptides,
        executables=tool_executables or {},
        hla_alleles=read_hla_list(summary.hla_alleles_txt),
    )

    provenance_registry = ProvenanceRegistry()
    net_path = mhc_path = stab_path = None
    if not skip_netmhcpan:
        net_raw = outdir / "tools" / "netmhcpan.xls"
        run_netmhcpan(ctx, net_raw)
        net_path = pres_dir / "netmhcpan_evidence.tsv"
        net_prov = ctx.tool_provenance.get("netmhcpan")
        if net_prov:
            provenance_registry.set(net_prov)
        write_netmhcpan_evidence(net_path, parse_netmhcpan(net_raw, sample_id), net_prov)
    if not skip_mhcflurry:
        mhc_raw = outdir / "tools" / "mhcflurry.csv"
        run_mhcflurry(ctx, mhc_raw)
        mhc_path = pres_dir / "mhcflurry_evidence.tsv"
        if not stub:
            provenance_registry.register_real("mhcflurry", mhc_raw)
        else:
            provenance_registry.register_stub("mhcflurry")
        write_mhcflurry_evidence(mhc_path, parse_mhcflurry(mhc_raw, sample_id))
    if not skip_stabpan:
        stab_raw = outdir / "tools" / "netmhcstabpan.tsv"
        run_netmhcstabpan(ctx, stab_raw)
        stab_path = pres_dir / "netmhcstabpan_evidence.tsv"
        provenance_registry.register_real("netmhcstabpan", stab_raw) if not stub else provenance_registry.register_stub("netmhcstabpan")
        write_netmhcstabpan_evidence(stab_path, parse_netmhcstabpan(stab_raw, sample_id))

    pres_path = pres_dir / "presentation_evidence.tsv"
    build_presentation_evidence(
        raw_peptides, net_path, mhc_path, profile, pres_path, stab_path,
        provenance_registry=provenance_registry,
    )
    pres_rows = read_tsv(pres_path)

    immuno_paths: dict[str, Path | None] = {}
    if not (skip_prime and skip_bigmhc_im and skip_deepimmuno):
        skip_sources = set()
        if skip_prime:
            skip_sources.add("prime")
        if skip_bigmhc_im:
            skip_sources.add("bigmhc_im")
        if skip_deepimmuno:
            skip_sources.add("deepimmuno")
        immuno_paths = run_immunogenicity_predictors(
            raw_peptides,
            outdir,
            profile,
            ctx,
            skip=stub,
            skip_sources=skip_sources,
            provenance_registry=provenance_registry,
        )
    apply_immunogenicity_evidence(pres_rows, immuno_paths, profile)
    write_tsv(pres_path, pres_rows, PRESENTATION_FIELDS)

    predictions = outdir / "peptide_predictions.tsv"
    write_tsv(predictions, pres_rows, PRESENTATION_FIELDS)

    manifest = outdir / "peptide_predict_summary.json"
    write_json(
        manifest,
        {
            "sample_id": sample_id,
            "profile": profile["_profile_name"],
            "input": summary.__dict__,
            "outputs": {
                "raw_peptides": str(raw_peptides),
                "pairs_tsv": summary.pairs_tsv,
                "presentation_evidence": str(pres_path),
                "peptide_predictions": str(predictions),
                **{k: str(v) for k, v in immuno_paths.items() if v},
            },
            "skipped": {
                "netmhcpan": skip_netmhcpan,
                "mhcflurry": skip_mhcflurry,
                "prime": skip_prime,
                "bigmhc_im": skip_bigmhc_im,
                "deepimmuno": skip_deepimmuno,
                "netmhcstabpan": skip_stabpan,
            },
        },
    )
    return {
        "summary": summary,
        "raw_peptides": str(raw_peptides),
        "presentation_evidence": str(pres_path),
        "peptide_predictions": str(predictions),
        "manifest": str(manifest),
    }


def read_hla_list(path: str | Path) -> list[str]:
    alleles = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            alleles.append(line)
    return alleles
