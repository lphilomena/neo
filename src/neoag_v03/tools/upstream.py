from __future__ import annotations
from pathlib import Path
import gzip
import os
import tomllib

from ..adapters.pvacseq_enrich import enrich_pvacseq_outputs, pvacseq_enrich_enabled
from ..adapters.pvactools_parser import parse_pvactools_outputs
from ..adapters.splice_junction_adapter import run_splice_junction_upstream
from ..adapters.variant_peptide_adapter import (
    run_variant_peptide_upstream,
    variant_peptide_extraction_enabled,
)
from ..config import load_profile
from ..input_router import resolve_entry_mode
from ..preflight import vcf_preflight
from .registry import RunContext, ROOT
from .runner import run_tool
from .postprocess import lohhla_to_hla_loh_tsv, spechla_to_hla_loh_tsv, facets_to_cnv_tsv


def load_run_config(path: str | Path) -> dict:
    p = Path(path)
    return tomllib.loads(p.read_text(encoding="utf-8"))


def _cfg_get(cfg: dict, *keys, default=None):
    cur = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def _path_or_none(val) -> Path | None:
    if not val:
        return None
    p = Path(val)
    if not p.is_absolute():
        p = (ROOT / p).resolve()
    return p


def _open_vcf_text(path: Path):
    """Open VCF text regardless of whether a .gz suffix is actually compressed."""
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def vcf_has_csq_annotations(path: str | Path) -> bool:
    """Return true when the VCF carries a VEP CSQ INFO definition."""
    p = Path(path)
    if not p.is_file():
        return False
    with _open_vcf_text(p) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                return True
            if line.startswith("#CHROM"):
                return False
    return False


def _auto_annotate_variants_vcf(
    cfg: dict,
    *,
    variants_vcf: Path,
    tools_dir: Path,
    sample_id: str,
) -> tuple[Path, dict[str, str]]:
    """Annotate an unannotated somatic VCF with VEP before peptide extraction."""
    inputs = _cfg_get(cfg, "inputs", default={}) or {}
    if vcf_has_csq_annotations(variants_vcf):
        return variants_vcf, {}
    if inputs.get("auto_vep_annotate", True) is False:
        return variants_vcf, {}

    from ..vep.annotate import run_vep_pvacseq_annotate

    ref_fasta = _path_or_none(inputs.get("reference_fasta")) or _path_or_none(
        os.environ.get("NEOAG_REFERENCE_FASTA")
    )
    if not ref_fasta or not ref_fasta.is_file():
        for candidate in (
            str(ROOT / "data" / "ref" / "hg38" / "Homo_sapiens_assembly38.fasta"),
        ):
            p = Path(candidate)
            if p.is_file():
                ref_fasta = p
                break
    if not ref_fasta or not ref_fasta.is_file():
        raise ValueError(
            "variant_peptide_extraction received a VCF without VEP CSQ annotations, "
            "but automatic VEP annotation needs inputs.reference_fasta or "
            "NEOAG_REFERENCE_FASTA to point to a valid FASTA."
        )

    out_vcf = _path_or_none(inputs.get("vep_annotated_vcf") or inputs.get("annotated_vcf"))
    if not out_vcf:
        out_vcf = tools_dir / f"{sample_id}.vep.annotated.vcf"
    if out_vcf.is_file() and vcf_has_csq_annotations(out_vcf):
        return out_vcf, {"vep_annotated_vcf": str(out_vcf)}

    cache_dir = _path_or_none(inputs.get("vep_cache") or inputs.get("cache_dir")) or _path_or_none(
        os.environ.get("NEOAG_VEP_CACHE")
    )
    plugins_dir = _path_or_none(inputs.get("vep_plugins") or inputs.get("plugins_dir")) or _path_or_none(
        os.environ.get("NEOAG_VEP_PLUGINS")
    )
    workdir = tools_dir / "vep_annotate_work"
    workdir.mkdir(parents=True, exist_ok=True)

    result = run_vep_pvacseq_annotate(
        input_vcf=variants_vcf,
        output_vcf=out_vcf,
        reference_fasta=ref_fasta,
        workdir=workdir,
        cache_dir=cache_dir if cache_dir and cache_dir.exists() else None,
        plugins_dir=plugins_dir if plugins_dir and plugins_dir.exists() else None,
        online=True if inputs.get("vep_online") else None,
        fork=int(inputs.get("vep_fork", 4) or 4),
        pick=not bool(inputs.get("vep_no_pick", False)),
        expression_custom=_path_or_none(inputs.get("expression_custom")),
        index_vcf=False,
        vep_bin=inputs.get("vep_bin") or os.environ.get("NEOAG_VEP_BIN"),
        cache_version=inputs.get("vep_cache_version") or os.environ.get("NEOAG_VEP_CACHE_VERSION"),
    )
    if not vcf_has_csq_annotations(out_vcf):
        raise ValueError(f"Automatic VEP annotation did not produce CSQ annotations in {out_vcf}")
    return out_vcf, {"vep_annotated_vcf": str(out_vcf), "vep_annotate_command": result.get("command", "")}


def build_context(cfg: dict, outdir: Path) -> RunContext:
    tools = _cfg_get(cfg, "tools", default={}) or {}
    inputs = _cfg_get(cfg, "inputs", default={}) or {}
    ex = _cfg_get(tools, "executables", default={}) or {}
    hla = _cfg_get(inputs, "hla_alleles", default=[]) or []
    sample_id = _cfg_get(cfg, "sample", "id", default="SAMPLE001")
    parsed = outdir / "parsed"
    raw_pep = parsed / "raw_peptides.tsv"
    return RunContext(
        sample_id=sample_id,
        outdir=outdir,
        stub=bool(tools.get("stub", False)),
        executables={k: str(v) for k, v in ex.items()},
        hla_alleles=list(hla),
        raw_peptides=raw_pep if raw_pep.is_file() else _path_or_none(inputs.get("raw_peptides")),
        tumor_vcf=_path_or_none(inputs.get("tumor_vcf")),
        normal_vcf=_path_or_none(inputs.get("normal_vcf")),
        tumor_sample_name=inputs.get("tumor_sample_name"),
        normal_sample_name=inputs.get("normal_sample_name"),
        phased_vcf=_path_or_none(inputs.get("phased_vcf")),
        prediction_algorithms=inputs.get("prediction_algorithms") or "NetMHCpan",
        fusion_tsv=_path_or_none(inputs.get("fusion_tsv")),
        expression_tsv=_path_or_none(inputs.get("expression_tsv")),
        variants_vcf=_path_or_none(inputs.get("variants_vcf")),
        splice_junction_tsv=_path_or_none(inputs.get("splice_junction_tsv") or inputs.get("regtools_tsv")),
        reference_fasta=_path_or_none(inputs.get("reference_fasta")),
        gencode_gtf=_path_or_none(inputs.get("gencode_gtf")),
        facets_rds=_path_or_none(inputs.get("facets_rds")),
        lohhla_prediction=_path_or_none(inputs.get("lohhla_prediction")),
        pass_only=bool(inputs.get("pass_only", True)),
    )


def run_upstream(config_path: str | Path, outdir: str | Path | None = None) -> dict[str, str]:
    """Run enabled upstream bioinformatics tools; returns paths for downstream v0.3."""
    cfg = load_run_config(config_path)
    out = Path(outdir or _cfg_get(cfg, "sample", "outdir", default="results/upstream"))
    out.mkdir(parents=True, exist_ok=True)
    tools_dir = out / "tools"
    tools_dir.mkdir(exist_ok=True)
    parsed = out / "parsed"
    parsed.mkdir(exist_ok=True)

    ctx = build_context(cfg, out)
    enabled = _cfg_get(cfg, "tools", "enabled", default=[]) or []
    outputs: dict[str, str] = {}

    profile = _cfg_get(cfg, "sample", "profile", default="default")
    profile_name = load_profile(profile)["_profile_name"]
    sample_id = ctx.sample_id
    inputs_cfg = _cfg_get(cfg, "inputs", default={}) or {}
    entry_mode = resolve_entry_mode(cfg)

    variants_vcf = (
        _path_or_none(inputs_cfg.get("variants_vcf"))
        or ctx.variants_vcf
        or ctx.tumor_vcf
    )
    use_variant_peptides = variant_peptide_extraction_enabled(cfg, variants_vcf)
    splice_tsv = _path_or_none(inputs_cfg.get("splice_junction_tsv") or inputs_cfg.get("regtools_tsv"))

    if use_variant_peptides and variants_vcf and variants_vcf.is_file():
        variants_vcf, vep_outputs = _auto_annotate_variants_vcf(
            cfg,
            variants_vcf=variants_vcf,
            tools_dir=tools_dir,
            sample_id=sample_id,
        )
        outputs.update(vep_outputs)
        ctx.variants_vcf = variants_vcf

    if entry_mode == "splice_junction" and splice_tsv:
        if not ctx.hla_alleles:
            raise ValueError("splice_junction mode requires inputs.hla_alleles")
        pvacsplice_out = None
        if "pvacsplice" in enabled:
            pvacsplice_out = tools_dir / "pvacsplice_aggregated.tsv"
            try:
                run_tool("pvacsplice", ctx, pvacsplice_out)
                outputs["pvacsplice"] = str(pvacsplice_out)
            except (RuntimeError, FileNotFoundError, ValueError) as exc:
                if not variants_vcf:
                    raise
                pvacsplice_out = None
                outputs["pvacsplice_fallback"] = str(exc)
        splice_outputs = run_splice_junction_upstream(
            cfg,
            splice_path=splice_tsv,
            variants_vcf=variants_vcf,
            parsed_dir=parsed,
            tools_dir=tools_dir,
            sample_id=sample_id,
            profile_name=profile_name,
            hla_alleles=list(ctx.hla_alleles),
            pvacsplice_tsv=pvacsplice_out,
        )
        outputs.update(splice_outputs)
        ctx.raw_peptides = Path(splice_outputs["raw_peptides"])
        use_variant_peptides = False

    if use_variant_peptides:
        if not ctx.hla_alleles:
            raise ValueError("variant_peptide_extraction requires inputs.hla_alleles")
        vp_outputs = run_variant_peptide_upstream(
            cfg,
            variants_vcf=variants_vcf,
            parsed_dir=parsed,
            tools_dir=tools_dir,
            sample_id=sample_id,
            profile_name=profile_name,
            hla_alleles=list(ctx.hla_alleles),
        )
        outputs.update(vp_outputs)
        ctx.raw_peptides = Path(vp_outputs["raw_peptides"])
        outputs["peptide_source"] = "extract-variant-peptides"

    if "pvacseq" in enabled and not use_variant_peptides and entry_mode != "splice_junction":
        pvac_out = tools_dir / "pvacseq_aggregated.tsv"
        run_tool("pvacseq", ctx, pvac_out)
        outputs["pvacseq"] = str(pvac_out)

    if "pvacfuse" in enabled:
        fuse_out = tools_dir / "pvacfuse_aggregated.tsv"
        run_tool("pvacfuse", ctx, fuse_out)
        outputs["pvacfuse"] = str(fuse_out)

    pvac_paths = []
    if "pvacseq" in outputs:
        pvac_paths.append(outputs["pvacseq"])
    if "pvacfuse" in outputs:
        pvac_paths.append(outputs["pvacfuse"])
    pre_pvac = _cfg_get(cfg, "inputs", "pvac_files", default=[]) or []
    pvac_paths.extend([str(_path_or_none(p)) for p in pre_pvac if _path_or_none(p)])

    if pvac_paths and not use_variant_peptides:
        parse_pvactools_outputs(
            pvac_paths,
            sample_id,
            profile_name,
            parsed / "raw_events.tsv",
            parsed / "raw_peptides.tsv",
        )
        ctx.raw_peptides = parsed / "raw_peptides.tsv"
        outputs["peptide_source"] = "pvactools"

    if "netmhcpan" in enabled:
        p = tools_dir / "netmhcpan.xls"
        run_tool("netmhcpan", ctx, p)
        outputs["netmhcpan"] = str(p)

    if "mhcflurry" in enabled:
        p = tools_dir / "mhcflurry.csv"
        run_tool("mhcflurry", ctx, p)
        outputs["mhcflurry"] = str(p)

    if "netmhcstabpan" in enabled:
        p = tools_dir / "netmhcstabpan.tsv"
        run_tool("netmhcstabpan", ctx, p)
        outputs["netmhcstabpan"] = str(p)

    for tool_key, input_key in (
        ("netmhcpan", "netmhcpan"),
        ("mhcflurry", "mhcflurry"),
        ("netmhcstabpan", "netmhcstabpan"),
    ):
        if tool_key in outputs:
            continue
        cached = _path_or_none(_cfg_get(cfg, "inputs", input_key))
        if cached and cached.is_file():
            outputs[tool_key] = str(cached)

    pvacseq_agg_path = outputs.get("pvacseq")
    if (
        pvacseq_agg_path
        and variants_vcf
        and variants_vcf.is_file()
        and pvacseq_enrich_enabled(cfg, variants_vcf, has_pvacseq_output=True)
    ):
        vcf_preflight(variants_vcf).require_pvacseq_ready()
        enrich_inputs = dict(inputs_cfg)
        if outputs.get("netmhcpan"):
            enrich_inputs["netmhcpan"] = outputs["netmhcpan"]
        if outputs.get("mhcflurry"):
            enrich_inputs["mhcflurry"] = outputs["mhcflurry"]
        enrich_outputs = enrich_pvacseq_outputs(
            {**cfg, "inputs": enrich_inputs},
            aggregated_tsv=pvacseq_agg_path,
            variants_vcf=variants_vcf,
            raw_peptides_tsv=parsed / "raw_peptides.tsv",
            raw_events_tsv=parsed / "raw_events.tsv",
            out_enriched_tsv=tools_dir / "pvacseq_enriched.tsv",
            sample_id=sample_id,
            profile_name=profile_name,
        )
        outputs.update(enrich_outputs)

    if ctx.hla_alleles and outputs.get("variant_peptides"):
        from ..adapters.peptide_netmhcpan import (
            annotate_raw_peptides_tsv,
            annotate_variant_peptide_tsv,
        )

        if use_variant_peptides and ctx.raw_peptides and Path(ctx.raw_peptides).is_file():
            if outputs.get("netmhcpan") or outputs.get("mhcflurry"):
                ann_raw = annotate_raw_peptides_tsv(
                    ctx.raw_peptides,
                    outputs.get("netmhcpan"),
                    outputs.get("mhcflurry"),
                )
                outputs["raw_peptides_binding_annotated"] = str(ann_raw.get("rows", 0))

            ann = annotate_variant_peptide_tsv(
                outputs["variant_peptides"],
                list(ctx.hla_alleles),
                netmhcpan_xls=outputs.get("netmhcpan"),
                mhcflurry_csv=outputs.get("mhcflurry"),
                netmhcstabpan_tsv=outputs.get("netmhcstabpan"),
                output_tsv=tools_dir / "variant_peptides.annotated.tsv",
            )
            outputs["variant_peptides_annotated"] = ann["output_tsv"]

    if "vep" in enabled:
        p = tools_dir / "vep_appm.tsv"
        run_tool("vep", ctx, p)
        outputs["vep_appm"] = str(p)

    if "lohhla" in enabled:
        p = tools_dir / "hla_loh.tsv"
        run_tool("lohhla", ctx, p)
        outputs["hla_loh"] = str(p)
    elif "hla_loh" not in outputs:
        cached = _path_or_none(inputs_cfg.get("hla_loh"))
        if cached and cached.is_file():
            outputs["hla_loh"] = str(cached)
        else:
            pred = _path_or_none(inputs_cfg.get("lohhla_prediction")) or ctx.lohhla_prediction
            if pred and pred.is_file():
                p = tools_dir / "hla_loh.tsv"
                lohhla_to_hla_loh_tsv(pred, p)
                outputs["hla_loh"] = str(p)
            else:
                spechla_merge = _path_or_none(
                    inputs_cfg.get("spechla_merge")
                    or inputs_cfg.get("spechla_loh")
                    or inputs_cfg.get("spechla_loh_merge")
                )
                if spechla_merge and spechla_merge.is_file():
                    p = tools_dir / "hla_loh.tsv"
                    spechla_to_hla_loh_tsv(spechla_merge, p)
                    outputs["hla_loh"] = str(p)

    if "facets" in enabled:
        p = tools_dir / "purity.tsv"
        run_tool("facets", ctx, p)
        outputs["purity"] = str(p)
        work_cncf = ctx.outdir / "tools" / "facets" / "facets_cncf.tsv"
        if work_cncf.is_file():
            cnv_out = tools_dir / "cnv_segments.tsv"
            facets_to_cnv_tsv(work_cncf, cnv_out)
            outputs["cnv"] = str(cnv_out)
    elif "purity" not in outputs:
        cached = _path_or_none(inputs_cfg.get("purity"))
        if cached and cached.is_file():
            outputs["purity"] = str(cached)
        elif ctx.facets_rds and ctx.facets_rds.is_file():
            p = tools_dir / "purity.tsv"
            run_tool("facets", ctx, p)
            outputs["purity"] = str(p)
            work_cncf = ctx.outdir / "tools" / "facets" / "facets_cncf.tsv"
            if work_cncf.is_file() and "cnv" not in outputs:
                cnv_out = tools_dir / "cnv_segments.tsv"
                facets_to_cnv_tsv(work_cncf, cnv_out)
                outputs["cnv"] = str(cnv_out)

    if variants_vcf and variants_vcf.is_file() and inputs_cfg.get("extract_appm_from_vcf", True):
        from ..adapters.vcf_appm import (
            extract_gx_expression,
            write_expression_tsv,
            write_vep_appm_from_vcf,
        )

        if "vep_appm" not in outputs:
            vep_cached = _path_or_none(inputs_cfg.get("vep_appm"))
            if vep_cached and vep_cached.is_file():
                outputs["vep_appm"] = str(vep_cached)
            else:
                vep_path = tools_dir / "vep_appm.tsv"
                write_vep_appm_from_vcf(variants_vcf, vep_path)
                outputs["vep_appm"] = str(vep_path)

        if "expression" not in outputs:
            expr_cached = _path_or_none(inputs_cfg.get("expression"))
            if expr_cached and expr_cached.is_file():
                outputs["expression"] = str(expr_cached)
            else:
                expr_path = tools_dir / "gene_expression.tsv"
                tumor_sample = inputs_cfg.get("tumor_sample_name") or ctx.tumor_sample_name
                write_expression_tsv(
                    extract_gx_expression(variants_vcf, tumor_sample=tumor_sample),
                    expr_path,
                )
                outputs["expression"] = str(expr_path)

    for key in ("purity", "hla_loh", "expression", "cnv", "normal_expression", "normal_hla_ligands"):
        if key in outputs:
            continue
        val = inputs_cfg.get(key) or inputs_cfg.get(f"{key}_tsv")
        p = _path_or_none(val)
        if p and p.is_file():
            outputs[key] = str(p)

    cached_cncf = _path_or_none(inputs_cfg.get("facets_cncf"))
    if cached_cncf and cached_cncf.is_file() and "cnv" not in outputs:
        outputs["cnv"] = str(cached_cncf)

    outputs["outdir"] = str(out)
    if ctx.raw_peptides and ctx.raw_peptides.is_file():
        outputs["raw_peptides"] = str(ctx.raw_peptides)
    if (parsed / "raw_events.tsv").is_file():
        outputs["raw_events"] = str(parsed / "raw_events.tsv")

    return outputs
