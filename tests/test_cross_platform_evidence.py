from neoag.cross_platform import (
    annotate_event_cross_platform,
    canonical_variant_key,
    propagate_cross_platform_to_peptide,
)


def test_canonical_variant_key_normalizes_chr_and_simple_indel():
    assert canonical_variant_key("chr12", "50350977", "GC", "G") == "12:50350977:GC:G"
    assert canonical_variant_key("chr1", 10, "AC", "AT") == "1:11:C:T"


def test_common_event_is_annotated_without_penalty():
    event = {"mutation_source": "SNV", "chrom": "chr12", "pos": "25245350", "ref": "C", "alt": "T"}
    evidence = {
        "12:25245350:C:T": {
            "variant_key": "12:25245350:C:T",
            "comparison_status": "COMMON",
            "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
        }
    }
    result = annotate_event_cross_platform(event, evidence, {"cross_platform": {}})
    assert result["cross_platform_confidence"] == "HIGH"
    assert result["cross_platform_multiplier"] == "1.0000"
    assert result["cross_platform_review_required"] == "no"


def test_normal_support_caps_and_propagates_to_peptide():
    event = {"mutation_source": "InDel", "chrom": "chr1", "pos": 20, "ref": "GA", "alt": "G"}
    evidence = {
        "1:20:GA:G": {
            "variant_key": "1:20:GA:G",
            "comparison_status": "WES_ONLY",
            "cross_platform_status": "NORMAL_SUPPORT_REVIEW",
            "normal_depth": "100",
            "normal_alt_count": "12",
        }
    }
    profile = {"cross_platform": {"normal_support_multiplier": 0.2}}
    result = annotate_event_cross_platform(event, evidence, profile)
    peptide = propagate_cross_platform_to_peptide({}, result)
    assert result["cross_platform_priority_cap"] == "D"
    assert result["cross_platform_multiplier"] == "0.2000"
    assert peptide["normal_alt_count"] == "12"
    assert peptide["cross_platform_review_required"] == "yes"


def test_power_limited_absence_is_not_penalized():
    event = {"mutation_source": "SNV", "chrom": "chr3", "pos": 30, "ref": "A", "alt": "G"}
    evidence = {
        "3:30:A:G": {
            "variant_key": "3:30:A:G",
            "comparison_status": "WES_ONLY",
            "cross_platform_status": "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF",
        }
    }
    result = annotate_event_cross_platform(event, evidence, {"cross_platform": {}})
    assert result["cross_platform_multiplier"] == "1.0000"
    assert result["cross_platform_confidence"] == "LOW_UNASSESSED"


def test_phased_event_inherits_concordant_component_evidence():
    event = {
        "mutation_source": "SNV",
        "chrom": "chr2",
        "pos": 100,
        "ref": "C;A",
        "alt": "A;C",
        "component_event_ids": "TBR1|chr2:100C>A;TBR1|chr2:102A>C",
    }
    evidence = {
        "2:100:C:A": {"variant_key": "2:100:C:A", "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT"},
        "2:102:A:C": {"variant_key": "2:102:A:C", "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT"},
    }
    result = annotate_event_cross_platform(event, evidence, {"cross_platform": {}})
    assert result["cross_platform_status"] == "CROSS_PLATFORM_PASS_CONCORDANT"
    assert result["cross_platform_variant_key"] == "2:100:C:A;2:102:A:C"
