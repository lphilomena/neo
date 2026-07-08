from __future__ import annotations

import pytest

from neoag_v03.preflight import vcf_preflight


def test_vcf_preflight_detects_csq_and_gt(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele|Consequence|SYMBOL\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|TP53\tGT:AD\t0/1:10,5\n",
        encoding="utf-8",
    )
    status = vcf_preflight(vcf)
    assert status.pvacseq_ready
    assert status.variant_rows == 1
    assert status.sample_columns == 1


def test_vcf_preflight_reports_missing_csq(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    status = vcf_preflight(vcf)
    assert not status.pvacseq_ready
    with pytest.raises(ValueError, match="VEP CSQ header"):
        status.require_pvacseq_ready()
