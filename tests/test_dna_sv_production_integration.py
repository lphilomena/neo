from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from neoag.controlled_execution.doctor import CheckRow
from neoag.open_neo.install_check import _assess_tier


ROOT = Path(__file__).resolve().parents[1]
LINKER = ROOT / "scripts/link_dna_sv_rna_fusions.py"
GENERATOR = ROOT / "scripts/generate_production_from_results_manifest.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_linker(tmp_path: Path, fusion_rows: list[dict[str, str]], sv_rows: list[dict[str, str]]):
    events = tmp_path / "fusion.events.tsv"
    peptides = tmp_path / "fusion.peptides.tsv"
    union = tmp_path / "fusion.union.tsv"
    sv = tmp_path / "sv.events.tsv"
    out = tmp_path / "out"
    write_tsv(events, [{"event_id": row["event_id"], "gene": row["gene_pair"]} for row in fusion_rows])
    write_tsv(peptides, [{"event_id": row["event_id"], "peptide_id": f"P_{row['event_id']}", "peptide": "AAAAAAAAA"} for row in fusion_rows])
    write_tsv(union, fusion_rows)
    write_tsv(sv, sv_rows)
    subprocess.run(
        [sys.executable, str(LINKER), "--fusion-events", str(events), "--fusion-peptides", str(peptides), "--fusion-union", str(union), "--sv-events", str(sv), "--outdir", str(out)],
        check=True,
    )
    return read_tsv(out / "raw_events.tsv"), read_tsv(out / "raw_peptides.tsv"), read_tsv(out / "dna_sv_rna_fusion_links.tsv")


def sv(event_id: str, pos1: int, pos2: int, gene1: str = "EWSR1", gene2: str = "WT1") -> dict[str, str]:
    return {
        "sv_event_id": event_id,
        "adjacency_key": f"GRCh38|chr22:{pos1}:+|chr11:{pos2}:-",
        "gene1": gene1,
        "gene2": gene2,
        "chrom1": "chr22",
        "pos1": str(pos1),
        "strand1": "+",
        "chrom2": "chr11",
        "pos2": str(pos2),
        "strand2": "-",
        "callers": "Manta;GRIDSS2",
        "caller_count": "2",
        "tumor_sr": "7",
        "tumor_pe": "11",
        "tumor_alt_support": "18",
    }


def test_exact_adjacency_is_attached_to_events_and_peptides(tmp_path: Path):
    fusion = [{"event_id": "F1", "gene_pair": "EWSR1::WT1", "left_breakpoint": "chr22:100:+", "right_breakpoint": "chr11:200:-"}]
    events, peptides, links = run_linker(tmp_path, fusion, [sv("SV1", 101, 198)])
    assert events[0]["dna_sv_confirmation_status"] == "SUPPORTED"
    assert events[0]["dna_sv_match_method"] == "EXACT_ADJACENCY"
    assert events[0]["dna_sv_support_reads"] == "18"
    assert peptides[0]["dna_sv_event_id"] == "SV1"
    assert links[0]["dna_sv_callers"] == "Manta;GRIDSS2"


def test_unique_gene_chromosome_projection_does_not_equate_intron_and_exon_coordinates(tmp_path: Path):
    fusion = [{"event_id": "F1", "gene_pair": "EWSR1::WT1", "left_breakpoint": "chr22:100000:+", "right_breakpoint": "chr11:200000:-"}]
    events, _, _ = run_linker(tmp_path, fusion, [sv("SV1", 150000, 250000)])
    assert events[0]["dna_sv_confirmation_status"] == "SUPPORTED"
    assert events[0]["dna_sv_match_method"] == "TRANSCRIPT_PROJECTION_UNIQUE"


def test_multiple_same_gene_pair_adjacencies_are_ambiguous(tmp_path: Path):
    fusion = [{"event_id": "F1", "gene_pair": "EWSR1::WT1", "left_breakpoint": "chr22:100000:+", "right_breakpoint": "chr11:200000:-"}]
    events, peptides, _ = run_linker(tmp_path, fusion, [sv("SV1", 150000, 250000), sv("SV2", 160000, 260000)])
    assert events[0]["dna_sv_confirmation_status"] == "AMBIGUOUS"
    assert not events[0]["dna_sv_event_id"]
    assert peptides[0]["dna_sv_match_method"] == "AMBIGUOUS"


def test_gene_pair_only_never_borrows_incompatible_dna_event(tmp_path: Path):
    fusion = [{"event_id": "F1", "gene_pair": "EWSR1::WT1", "left_breakpoint": "chr22:100:+", "right_breakpoint": "chr11:200:-"}]
    events, _, _ = run_linker(tmp_path, fusion, [sv("SV_OTHER", 100, 200, gene1="EWSR1", gene2="FLI1")])
    assert events[0]["dna_sv_confirmation_status"] == "NOT_DETECTED"
    assert not events[0]["dna_sv_event_id"]


def test_incompatible_orientation_is_not_projected(tmp_path: Path):
    fusion = [{"event_id": "F1", "gene_pair": "EWSR1::WT1", "left_breakpoint": "chr22:100:+", "right_breakpoint": "chr11:200:+"}]
    events, _, _ = run_linker(tmp_path, fusion, [sv("SV1", 500, 600)])
    assert events[0]["dna_sv_confirmation_status"] == "NOT_DETECTED"
    assert not events[0]["dna_sv_event_id"]


def test_production_entrypoints_expose_required_contracts():
    wrapper = (ROOT / "scripts/run_production_case.sh").read_text(encoding="utf-8")
    generator = (ROOT / "scripts/generate_production_from_results_manifest.py").read_text(encoding="utf-8")
    workflow = (ROOT / "workflows/sv_phase1_wgs.nf").read_text(encoding="utf-8")
    schemas = (ROOT / "src/neoag/schemas.py").read_text(encoding="utf-8")
    for token in ("--tumor-dna-bam", "--normal-dna-bam", "--assay-type", "--capture-bed", "--sv-vcf"):
        assert token in wrapper
        assert token in generator
    assert "sample_identity_bam_matcher" in generator
    assert "fusion_dna_sv_link" in generator
    assert "tumor_sample_name" in workflow and "normal_sample_name" in workflow
    assert "dna_sv_confirmation_status" in schemas
    assert "sample_identity_status" in schemas


def generator_inputs(tmp_path: Path) -> list[str]:
    paths = {}
    for name in ("profile.toml", "rules.toml", "ref.fa", "genes.gtf", "hla.txt", "somatic.vcf", "sv.vcf", "purity.tsv", "cnv.tsv", "lohhla.tsv", "spechla.tsv"):
        path = tmp_path / name
        path.write_text("HLA-A*02:01\n" if name == "hla.txt" else "x\n", encoding="utf-8")
        paths[name] = path
    return [
        sys.executable, str(GENERATOR), "--project-root", str(ROOT), "--sample-id", "S1",
        "--outdir", str(tmp_path / "run"), "--output", str(tmp_path / "production.toml"),
        "--profile", str(paths["profile.toml"]), "--evidence-consensus-rules", str(paths["rules.toml"]),
        "--reference-fasta", str(paths["ref.fa"]), "--gencode-gtf", str(paths["genes.gtf"]),
        "--hla-file", str(paths["hla.txt"]), "--somatic-vcf", str(paths["somatic.vcf"]),
        "--purity", str(paths["purity.tsv"]), "--cnv", str(paths["cnv.tsv"]),
        "--lohhla", str(paths["lohhla.tsv"]), "--spechla-loh", str(paths["spechla.tsv"]),
        "--sv-vcf", str(paths["sv.vcf"]), "--sv-caller", "Manta",
        "--tumor-sample-name", "TUMOR", "--normal-sample-name", "NORMAL", "--assay-type", "WGS",
    ]


def test_manifest_generates_dna_sv_stage_for_existing_vcf(tmp_path: Path):
    command = generator_inputs(tmp_path)
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    manifest = (tmp_path / "production.toml").read_text(encoding="utf-8")
    assert "[stages.dna_sv_discovery]" in manifest
    assert 'source = "DNA_SV"' in manifest
    assert "sv-build-raw" in manifest
    assert "--tumor-sample-name" in manifest and "--normal-sample-name" in manifest


def test_wes_manifest_rejects_missing_capture_bed(tmp_path: Path):
    command = generator_inputs(tmp_path)
    command[command.index("WGS")] = "WES"
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode != 0
    assert "requires --capture-bed" in result.stderr


def test_full_readiness_requires_every_dna_sv_group_member():
    rows = [CheckRow("tool", name, "OK") for name in ("manta", "svaba")]
    _, requirements = _assess_tier("full", rows, None)
    dna_sv = next(row for row in requirements if row["requirement"] == "dna_sv")
    assert dna_sv["kind"] == "tool_group_all"
    assert dna_sv["status"] == "MISSING"
    assert "missing=gridss" in dna_sv["evidence"]

    rows.append(CheckRow("tool", "gridss", "OK"))
    _, requirements = _assess_tier("full", rows, None)
    dna_sv = next(row for row in requirements if row["requirement"] == "dna_sv")
    assert dna_sv["status"] == "OK"
    assert "missing=" in dna_sv["evidence"]


def test_skill1_installs_complete_dna_sv_group():
    installer = (ROOT / "scripts/install_dna_sv_tools.sh").read_text(encoding="utf-8")
    deployer = (ROOT / "scripts/deploy_external_tools.sh").read_text(encoding="utf-8")
    env_yml = (ROOT / "conda/env.neoag-sv.yml").read_text(encoding="utf-8")
    gridss_yml = (ROOT / "conda/env.neoag-gridss.yml").read_text(encoding="utf-8")
    for tool in ("manta", "svaba"):
        assert tool in env_yml.lower()
    assert "gridss=2.13.2" in gridss_yml.lower()
    assert "configManta.py" in installer
    assert "DNA-SV capability group READY" in installer
    assert "with_bioc_data_cache.sh" in installer
    assert "genomeinfodbdata-1.2.13" in installer
    assert "Install DNA-SV group (Manta / SvABA / GRIDSS)" in deployer
    assert "SKIP_DNA_SV" in deployer
