#!/usr/bin/env python3
"""Dependency-free validation for the hardened SV adapter."""

from __future__ import annotations

import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from neoag.source_chain import derive_source_chain_confidence, source_chain_track  # noqa: E402
from neoag.sv.exact_evidence import load_exact_junction_evidence  # noqa: E402
from neoag.sv.identity import canonical_breakpoint_key  # noqa: E402
from neoag.sv.phase1 import build_sv_phase1_raw  # noqa: E402
from neoag.sv.sv_callset import parse_vcf_records  # noqa: E402
from neoag.sv.sv_merge import cluster_sv_records  # noqa: E402
from neoag.utils import read_tsv  # noqa: E402


FX = ROOT / "data" / "fixtures_sv"


def validate() -> None:
    try:
        parse_vcf_records(FX / "mini_sv.vcf", "GRIDSS2")
    except ValueError:
        pass
    else:
        raise AssertionError("multi-sample VCF sample-order inference was not blocked")

    key1 = canonical_breakpoint_key("GRCh38", "chr1", 60, "+", "chr2", 31, "-")
    key2 = canonical_breakpoint_key("grch38", "2", 31, "-", "1", 60, "+")
    assert key1 == key2
    assert load_exact_junction_evidence(FX / "rna_junctions_exact.tsv", default_build="GRCh38")[key1].qc_pass

    records = parse_vcf_records(
        FX / "mini_sv.vcf",
        "GRIDSS2",
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
    )
    assert len(records) == 1
    assert len(cluster_sv_records(records)) == 1

    with TemporaryDirectory() as tmp:
        out = build_sv_phase1_raw(
            sample_id="SVMINI",
            sv_vcfs=[FX / "mini_sv.vcf"],
            callers=["GRIDSS2"],
            tumor_sample_name="TUMOR",
            normal_sample_name="NORMAL",
            reference_fasta=FX / "mini_ref.fa",
            gencode_gtf=FX / "mini.gtf",
            hla=FX / "hla.txt",
            outdir=tmp,
            rna_junction_tsv=FX / "rna_junctions_exact.tsv",
            expressed_products_tsv=FX / "expressed_products.tsv",
        )
        events = read_tsv(out["sv_events_full"])
        peptides = read_tsv(out["raw_peptides"])
        assert events[0]["rna_evidence_match"] == "EXACT_BREAKPOINT"
        assert events[0]["rna_evidence_qc"] == "PASS"
        assert peptides

    with TemporaryDirectory() as tmp:
        out = build_sv_phase1_raw(
            sample_id="SVMINI",
            sv_vcfs=[FX / "mini_sv.vcf"],
            callers=["GRIDSS2"],
            tumor_sample_name="TUMOR",
            normal_sample_name="NORMAL",
            reference_fasta=FX / "mini_ref.fa",
            gencode_gtf=FX / "mini.gtf",
            hla=FX / "hla.txt",
            outdir=tmp,
        )
        assert read_tsv(out["raw_peptides"]) == []
        assert read_tsv(out["sv_events_full"])[0]["reconstruction_status"] == "hypothesis_only"

    assert source_chain_track({"mutation_source": "SV", "event_type": "SV_Fusion"}) == "DNA_SV"
    assert derive_source_chain_confidence({"source_chain_track": "DNA_SV"}).track == "DNA_SV"


if __name__ == "__main__":
    validate()
    print("SV hardening validation: PASS")
