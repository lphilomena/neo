from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from neoag.open_neo.review import _fusion_report_data
from neoag.open_neo.tool_consensus import _write_fusion


ROOT = Path(__file__).resolve().parents[1]


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_union(tmp_path: Path, *args: str) -> Path:
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    outdir = tmp_path / "union"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/build_fusion_caller_union.py"),
            "--sample-id", "S1", "--hla-file", str(hla), *args,
            "--no-targeted-fusion-rescue", "--outdir", str(outdir),
        ],
        cwd=ROOT,
        check=True,
    )
    return outdir


def test_easyfuse_pass_and_embedded_callers_have_separate_roles(tmp_path: Path) -> None:
    header = (
        "BPID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;prediction_class;prediction_prob;"
        "starfusion_detected;fusioncatcher_detected;arriba_detected;ft_junc_cnt;ft_anch_cnt;"
        "frame;type;neo_peptide_sequence;neo_peptide_sequence_bp\n"
    )
    easyfuse = write(
        tmp_path / "easyfuse/fusions.pass.csv",
        header
        + "bp1;EWSR1_WT1;chr22:100:+;chr11:200:-;tx1;positive;0.95;1;0;1;12;20;"
        "in_frame;trans;ACDEFGHIKLMNPQRSTVWY;10\n",
    )
    outdir = run_union(tmp_path, "--easyfuse", str(easyfuse))
    consensus = read_tsv(outdir / "fusion_consensus.tsv")
    assert len(consensus) == 1
    assert consensus[0]["fixed_panel_version"] == "OPEN_NEO_SHORT_READ_FUSION_V1"
    assert consensus[0]["easyfuse_pass"] == "yes"
    assert consensus[0]["fixed_panel_support_callers"] == "Arriba,STAR-Fusion"
    assert consensus[0]["n_fixed_panel_support"] == "2"
    assert consensus[0]["status"] == "SHORT_READ_MULTI_CALLER"
    availability = {row["caller"]: row for row in read_tsv(outdir / "fusion_caller_availability.tsv")}
    assert availability["FusionCatcher"]["availability_status"] == "AVAILABLE_EASYFUSE_EMBEDDED_COLUMNS"


def test_same_gene_pair_with_different_adjacencies_is_not_cross_validated(tmp_path: Path) -> None:
    star = write(
        tmp_path / "star.tsv",
        "#FusionName\tLeftBreakpoint\tRightBreakpoint\tJunctionReadCount\n"
        "EWSR1--WT1\tchr22:100:+\tchr11:200:-\t10\n",
    )
    arriba = write(
        tmp_path / "arriba.tsv",
        "gene1\tgene2\tbreakpoint1\tbreakpoint2\tsplit_reads1\n"
        "EWSR1\tWT1\tchr22:150:+\tchr11:250:-\t12\n",
    )
    outdir = run_union(tmp_path, "--star-fusion", str(star), "--arriba", str(arriba))
    consensus = read_tsv(outdir / "fusion_consensus.tsv")
    assert len(consensus) == 2
    assert {row["status"] for row in consensus} == {"SHORT_READ_SINGLE_CALLER"}
    assert all(row["n_fixed_panel_support"] == "1" for row in consensus)


def test_same_oriented_adjacency_is_cross_validated_across_fixed_callers(tmp_path: Path) -> None:
    star = write(
        tmp_path / "star.tsv",
        "#FusionName\tLeftBreakpoint\tRightBreakpoint\tJunctionReadCount\n"
        "EWSR1--WT1\tchr22:100:+\tchr11:200:-\t10\n",
    )
    arriba = write(
        tmp_path / "arriba.tsv",
        "gene1\tgene2\tbreakpoint1\tbreakpoint2\tsplit_reads1\n"
        "EWSR1\tWT1\tchr22:100:+\tchr11:200:-\t12\n",
    )
    outdir = run_union(tmp_path, "--star-fusion", str(star), "--arriba", str(arriba))
    consensus = read_tsv(outdir / "fusion_consensus.tsv")
    assert len(consensus) == 1
    assert consensus[0]["status"] == "SHORT_READ_MULTI_CALLER"
    assert consensus[0]["fixed_panel_support_callers"] == "Arriba,STAR-Fusion"


def test_dna_sv_support_is_written_back_to_authoritative_consensus(tmp_path: Path) -> None:
    events = write(
        tmp_path / "events.tsv",
        "event_id\tgene\nF1\tEWSR1::WT1\n",
    )
    peptides = write(tmp_path / "peptides.tsv", "event_id\tpeptide\nF1\tACDEFGHIK\n")
    union = write(
        tmp_path / "caller_union.tsv",
        "event_id\tgene_pair\tleft_breakpoint\tright_breakpoint\n"
        "F1\tEWSR1::WT1\tchr22:100:+\tchr11:200:-\n",
    )
    consensus = write(
        tmp_path / "consensus.tsv",
        "event_id\tmember_event_ids\tadjacency_key\tfusion\tdna_sv_support\n"
        "F1\tF1\tGRCh38|EWSR1::WT1|chr22:100|chr11:200|.\tEWSR1::WT1\tUNASSESSED\n",
    )
    sv = write(
        tmp_path / "sv.tsv",
        "event_id\tgene1\tgene2\tchrom1\tpos1\tstrand1\tchrom2\tpos2\tstrand2\tcallers\tcaller_count\n"
        "SV1\tEWSR1\tWT1\tchr22\t100\t+\tchr11\t200\t-\tManta,SvABA\t2\n",
    )
    outdir = tmp_path / "linked"
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/link_dna_sv_rna_fusions.py"),
         "--fusion-events", str(events), "--fusion-peptides", str(peptides),
         "--fusion-union", str(union), "--fusion-consensus", str(consensus),
         "--sv-events", str(sv), "--outdir", str(outdir)],
        cwd=ROOT,
        check=True,
    )
    row = read_tsv(outdir / "fusion_consensus.tsv")[0]
    assert row["dna_sv_support"] == "SUPPORTED"
    assert row["dna_sv_callers"] == "Manta,SvABA"


def test_skill3_summary_uses_canonical_union_not_legacy_directory_scan(tmp_path: Path) -> None:
    fusion = tmp_path / "pipeline/production/branches/fusion"
    write(
        fusion / "intermediates/fusion_caller_union.tsv",
        "event_id\tsource_tool\tadmission_policy\n"
        "F1\tEasyFuse\tCALLER_PASS\nF1\tArriba\tCALLER_PASS\nF1\tSTAR-Fusion\tCALLER_PASS\n",
    )
    write(
        fusion / "consensus/fusion_caller_availability.tsv",
        "caller\tavailability_status\tsource_files\n"
        "Arriba\tAVAILABLE_EASYFUSE_EMBEDDED_COLUMNS\ta.csv\n"
        "STAR-Fusion\tAVAILABLE_EASYFUSE_EMBEDDED_COLUMNS\ts.csv\n"
        "FusionCatcher\tAVAILABLE_EASYFUSE_EMBEDDED_COLUMNS\tf.csv\n",
    )
    write(
        fusion / "consensus/fusion_consensus.tsv",
        "event_id\tfusion\tstatus\tfixed_panel_support_callers\tlong_read_support\tdna_sv_support\trescue_status\n"
        "F1\tEWSR1::WT1\tSHORT_READ_MULTI_CALLER\tArriba,STAR-Fusion\t\tSUPPORTED\tNONE\n",
    )
    write(fusion / "dna_sv_linked/raw_events.tsv", "event_id\tgene\nF1\tEWSR1::WT1\n")
    write(fusion / "dna_sv_linked/raw_peptides.tsv", "event_id\tpeptide\nF1\tACDEFGHIK\n")
    rows, _events, _note = _fusion_report_data(tmp_path)
    by_layer = {row["tool_or_layer"]: row for row in rows}
    assert by_layer["EasyFuse PASS (aggregator)"]["records"] == "1"
    assert by_layer["Fixed-panel multi-caller consensus"]["records"] == "1"
    assert by_layer["DNA-SV orthogonal support"]["records"] == "1"


def test_macro_tool_consensus_reuses_authoritative_adjacency_file(tmp_path: Path) -> None:
    canonical = write(
        tmp_path / "fusion_consensus.tsv",
        "event_id\tadjacency_key\tfusion\tfixed_panel_version\tstatus\n"
        "F1\tGRCh38|EWSR1::WT1|chr22:100|chr11:200|.\tEWSR1::WT1\t"
        "OPEN_NEO_SHORT_READ_FUSION_V1\tSHORT_READ_MULTI_CALLER\n",
    )
    overall, rows = _write_fusion(
        {"fusion_consensus_tsv": str(canonical)}, {"fusion": {}}, tmp_path / "out",
    )
    assert overall == "MULTI_CALLER_STRONG"
    assert rows[0]["adjacency_key"].startswith("GRCh38|")


def test_macro_tool_consensus_refuses_gene_pair_only_inference(tmp_path: Path) -> None:
    raw = write(tmp_path / "arriba.tsv", "gene1\tgene2\nEWSR1\tWT1\n")
    overall, rows = _write_fusion({}, {"fusion": {"arriba": str(raw)}}, tmp_path / "out")
    assert overall == "UNASSESSED"
    assert rows[0]["status"] == "UNASSESSED"
    assert "gene-pair-only consensus is forbidden" in rows[0]["reason"]
