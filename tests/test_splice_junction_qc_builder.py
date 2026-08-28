import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_splice_junction_qc_from_star_bam.py"
SPEC = importlib.util.spec_from_file_location("splice_junction_qc_builder", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_snaf_consensus_and_splicemutr_origin_form_exact_crossvalidation(tmp_path):
    consensus = tmp_path / "splice_consensus.tsv"
    consensus.write_text(
        "event_id\tcanonical_junction_id\tsupport_tools\tstatus\n"
        "SJ|GRCh38|chr1|101|200|+\tSJ|GRCh38|chr1|101|200|+\tRegTools;SNAF\tCROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION\n",
        encoding="utf-8",
    )

    snaf, splicemutr = MODULE.load_crossvalidated_keys(consensus)
    splicemutr.update(MODULE.splicemutr_origin_events([{
        "event_id": "SJ|GRCh38|chr1|101|200|+",
        "origin_peptide_id": "POR|fixture",
    }]))

    event_id = "SJ|GRCh38|chr1|101|200|+"
    assert event_id in snaf
    assert event_id in splicemutr


def test_event_junction_prefers_canonical_and_accepts_structured_coordinates():
    assert MODULE.event_junction({
        "canonical_junction_id": "SJ|GRCh38|chr2|101|200|-",
        "event_name": "JUNC0001",
    }) == ("2", 101, 200, "-")
    assert MODULE.event_junction({
        "junction_chrom": "chr3", "junction_start": "301",
        "junction_end": "400", "junction_strand": "+",
    }) == ("3", 301, 400, "+")
