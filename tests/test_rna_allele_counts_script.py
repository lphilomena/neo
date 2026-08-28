from __future__ import annotations

import importlib.util
from pathlib import Path


def _module():
    path = Path.cwd() / "scripts/rna_allele_counts_pysam.py"
    spec = importlib.util.spec_from_file_location("rna_allele_counts_pysam", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_classifies_snv_and_indel_observations():
    module = _module()
    assert module.classify_observation(query_sequence="AT", query_position=0, indel=0, ref="A", alt="T") == "ref"
    assert module.classify_observation(query_sequence="TT", query_position=0, indel=0, ref="A", alt="T") == "alt"
    assert module.classify_observation(query_sequence="ATG", query_position=0, indel=2, ref="A", alt="ATG") == "alt"
    assert module.classify_observation(query_sequence="A", query_position=0, indel=-2, ref="ATG", alt="A") == "alt"
    assert module.classify_observation(query_sequence="A", query_position=0, indel=0, ref="ATG", alt="A") == "ref"


def test_complex_variant_is_explicitly_unassessed():
    module = _module()
    assert module.variant_type("AT", "GCAT") == "COMPLEX"
    assert module.classify_observation(query_sequence="AT", query_position=0, indel=0, ref="AT", alt="GCAT") == "other"
