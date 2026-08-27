from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
import importlib.util
from pathlib import Path

from neoag.splice.consensus import build_consensus
from neoag.splice.normal_background import parse_normal_junctions
from neoag.splice.normalization import normalize_splice_sources
from neoag.utils import read_tsv


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


class SpliceMutrNormalP0Tests(unittest.TestCase):
    def test_catalog_nonmembership_does_not_increase_composite_score(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "prepare_snaf_splice_branch", SCRIPTS / "prepare_snaf_splice_branch.py"
        )
        assert spec and spec.loader
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        base = {
            "junction_count": "20", "tumor_specificity_mean": "0.1",
            "tumor_specificity_mle": "0.05", "binding_affinity": "0.4",
            "immunogenicity": "0.7", "frame_evidence_status": "TRANSCRIPT_EVIDENCE_PRESENT",
        }
        not_listed = module.score_candidate({
            **base, "normal_junction_status": "NOT_LISTED_IN_NORMAL_CATALOG"
        })
        unassessed = module.score_candidate({
            **base, "normal_junction_status": "UNASSESSED"
        })
        self.assertEqual(not_listed, unassessed)

    def test_recount3_catalog_positive_is_not_read_as_zero(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            source = write(
                Path(td) / "normal.tsv",
                "junction_id\tchromosome\tstart\tend\tstrand\tgenome_build\t"
                "normal_samples\tnormal_reads\tnormal_total_reads\ttotal_samples\t"
                "sample_prevalence\tnormal_tissue_count\ttissue\tsource\tdataset\treference_release\n"
                "SJ|GRCh38|chr1|101|200|+\tchr1\t101\t200\t+\tGRCh38\t"
                "7\t19\t83\t251\t0.0278884462\t1\tLiver\t"
                "recount3_GTEx_v8\tGTEx_v8\trecount3_GTEx_v8_GRCh38\n",
            )
            bundle = parse_normal_junctions(
                source,
                sample_id="S1",
                allowed_junction_ids={"SJ|GRCh38|chr1|101|200|+"},
                strict=True,
            )
            row = bundle["normal_background"][0]
            self.assertEqual(row["detection_status"], "DETECTED")
            self.assertEqual(row["assessment_status"], "DETECTED_BROAD_NORMAL")
            self.assertEqual(row["junction_reads"], "19")
            self.assertEqual(row["normal_total_junction_reads"], "83")
            self.assertEqual(row["sample_count"], "7")
            self.assertEqual(row["total_samples"], "251")
            self.assertEqual(row["normal_tissues"], "Liver")
            self.assertEqual(row["source_dataset"], "GTEx_v8")
            self.assertEqual(row["normal_source_type"], "EXTERNAL_NORMAL_JUNCTION_CATALOG")

    def test_catalog_nonmembership_remains_n1_and_caps_r3(self) -> None:
        jid = "SJ|GRCh38|chr1|101|200|+"
        tables = {
            "events": [{"splice_event_id": "SEV|1", "alternative_junction_ids": jid}],
            "event_junction_links": [{"splice_event_id": "SEV|1", "junction_id": jid}],
            "orfs": [{
                "orf_id": "ORF|1", "splice_event_id": "SEV|1",
                "protein_sequence_sha256": "abc", "frame_status": "IN_FRAME",
                "orf_validity_status": "VALID_FULL_LENGTH", "source_generator": "SpliceMutr",
            }],
            "peptide_origins": [{
                "origin_peptide_id": "POR|1", "peptide_id": "PEP|1", "orf_id": "ORF|1",
                "splice_event_id": "SEV|1", "required_junction_ids": jid,
                "crosses_junction": "true", "contains_novel_aa": "true",
                "tumor_specificity_status": "UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT",
            }],
            "junction_read_qc": [{"junction_id": jid, "qc_status": "PASS", "resolution_status": "RESOLVED_EXACT"}],
            "tool_evidence": [{
                "entity_id": jid, "evidence_group": "RNA_JUNCTION", "source_assay_id": "RNA1",
                "verified_value": "10", "resolution_status": "RESOLVED_EXACT",
            }],
            "presentation": [{"origin_peptide_id": "POR|1"}],
            "normal_background": [],
        }
        row = build_consensus(tables, sample_id="S1")[0]
        self.assertEqual(row["normal_safety_grade"], "N1")
        self.assertEqual(row["normal_background_status"], "NORMAL_BACKGROUND_INCOMPLETE")
        self.assertEqual(row["final_evidence_tier"], "R3")
        self.assertIn("CAP_NORMAL_BACKGROUND_INCOMPLETE_R3", row["cap_codes"])

    def test_main_normalization_path_is_neutral_and_r3_capped(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            jid = "SJ|GRCh38|chr1|101|200|+"
            primary = write(
                tmp / "primary.tsv",
                "junction_id\tunique_split_reads\n" + jid + "\t10\n",
            )
            snaf = write(
                tmp / "snaf.tsv",
                "junction_id\tpeptide\thla_allele\tcrosses_junction\tcontains_novel_aa\tbinding_rank\n"
                + jid + "\tSYFPEITHI\tHLA-A*02:01\tyes\tyes\t0.2\n",
            )
            normal = write(
                tmp / "normal.tsv",
                "junction_id\tnormal_samples\tnormal_reads\n"
                "SJ|GRCh38|chr2|301|400|+\t2\t5\n",
            )
            paths = normalize_splice_sources(
                sample_id="S1", junctions=primary, outdir=tmp / "out",
                snaf=snaf, normal_junctions=normal, candidate_only=True,
            )
            event = read_tsv(Path(paths["raw_events"]))[0]
            peptide = read_tsv(Path(paths["raw_peptides"]))[0]
            self.assertEqual(event["normal_junction_status"], "NOT_LISTED_IN_NORMAL_CATALOG")
            self.assertEqual(event["tumor_specificity"], "0.5")
            self.assertEqual(event["cohort_analysis_status"], "UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT")
            self.assertEqual(event["priority_cap"], "R3")
            self.assertEqual(peptide["mutant_specificity_status"], "UNASSESSED")
            self.assertEqual(peptide["mutant_specificity_gate_status"], "REVIEW_REQUIRED")
            self.assertEqual(peptide["structural_novelty_status"], "ALTERED_JUNCTION_SPANNING_SEQUENCE")
            self.assertEqual(peptide["splice_consensus_tier"], "R3")

    def test_splicemutr_altered_path_is_not_called_tumor_specific(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            candidates = write(
                tmp / "candidates.tsv",
                "canonical_junction_id\tpeptide\tgene\n"
                "SJ|GRCh38|chr1|101|199|+\tCDE\tGENE1\n",
            )
            formed = tmp / "formed_transcripts"
            formed.mkdir()
            write(
                formed / "x_data_splicemutr_cp_corrected.txt",
                "chr\tstart\tend\tstrand\tpeptide\tpep_junc_loc\tmodified\ttx_id\tgene\tcluster\n"
                "chr1\t100\t200\t+\tABCDEFGH\t4\tTUMOR\tENST1\tGENE1\tC1\n",
            )
            out = tmp / "out"
            env = dict(os.environ)
            env["PYTHONPATH"] = str(ROOT / "src")
            subprocess.run(
                [
                    sys.executable, str(SCRIPTS / "rebuild_splice_origins_from_splicemutr.py"),
                    "--sample-id", "S1", "--candidates", str(candidates),
                    "--splicemutr-glob", str(formed / "*.txt"), "--outdir", str(out),
                ],
                check=True, env=env, capture_output=True, text=True,
            )
            origin = read_tsv(out / "splice_peptide_origins.tsv")[0]
            peptide = read_tsv(out / "raw_peptides.formal_origins.tsv")[0]
            self.assertEqual(origin["structural_novelty_status"], "ALTERED_JUNCTION_SPANNING_SEQUENCE")
            self.assertEqual(origin["tumor_specificity_status"], "UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT")
            self.assertEqual(peptide["mutant_specificity_status"], "UNASSESSED")
            self.assertEqual(peptide["mutant_specificity_gate_status"], "REVIEW_REQUIRED")
            self.assertEqual(peptide["mutant_specificity_priority_cap"], "R3")
            self.assertNotEqual(peptide["mutant_specificity_status"], "MT_SPECIFIC")

    def test_recount3_builder_records_release_and_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            rr = write(
                tmp / "junctions.tsv",
                "chromosome\tstart\tend\tstrand\tannotated\n1\t101\t200\t+\t1\n",
            )
            mm = write(
                tmp / "counts.mtx",
                "%%MatrixMarket matrix coordinate integer general\n1 3 2\n1 1 4\n1 3 2\n",
            )
            output = tmp / "liver.tsv.gz"
            subprocess.run(
                [
                    sys.executable, str(SCRIPTS / "build_recount3_normal_junctions.py"),
                    "--rr", str(rr), "--mm", str(mm), "--tissue", "Liver",
                    "--output", str(output),
                ],
                check=True, capture_output=True, text=True,
            )
            import gzip
            with gzip.open(output, "rt", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["normal_samples"], "2")
            self.assertEqual(row["total_samples"], "3")
            self.assertEqual(row["normal_reads"], "4")
            self.assertEqual(row["normal_total_reads"], "6")
            self.assertEqual(row["reference_release"], "recount3_GTEx_v8_GRCh38")

    def test_migration_downgrades_historical_mt_specific_splice_row(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            source = write(
                tmp / "old.tsv",
                "event_type\tmutation_source\tcontains_novel_aa\tmutant_specificity_status\t"
                "normal_junction_status\trna_evidence_completeness\n"
                "Splice\tSpliceMutr\tyes\tMT_SPECIFIC\tABSENT_GTEX_V11\tCOMPLETE\n",
            )
            output = tmp / "new.tsv"
            subprocess.run(
                [
                    sys.executable, str(SCRIPTS / "migrate_splicemutr_normal_p0.py"),
                    "--input", str(source), "--output", str(output), "--kind", "peptides",
                ],
                check=True, capture_output=True, text=True,
            )
            row = read_tsv(output)[0]
            self.assertEqual(row["normal_junction_status"], "NOT_LISTED_IN_NORMAL_CATALOG")
            self.assertEqual(row["mutant_specificity_status"], "UNASSESSED")
            self.assertEqual(row["cohort_analysis_status"], "UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT")
            self.assertEqual(row["normal_safety_grade"], "N1")
            self.assertEqual(row["splice_consensus_tier"], "R3")

    def test_pan_tissue_merge_preserves_denominator_semantics(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            header = (
                "junction_id\tchromosome\tstart\tend\tstrand\tgenome_build\tannotated\t"
                "normal_samples\tnormal_reads\tnormal_total_reads\ttotal_samples\t"
                "sample_prevalence\tnormal_tissue_count\ttissue\tsource\tdataset\treference_release\n"
            )
            liver = write(
                tmp / "liver.tsv",
                header + "SJ|GRCh38|chr1|101|200|+\tchr1\t101\t200\t+\tGRCh38\t1\t"
                "2\t5\t7\t10\t0.2\t1\tLiver\trecount3_GTEx_v8\tGTEx_v8\trecount3_GTEx_v8_GRCh38\n",
            )
            lung = write(
                tmp / "lung.tsv",
                header + "SJ|GRCh38|chr1|101|200|+\tchr1\t101\t200\t+\tGRCh38\t1\t"
                "3\t8\t12\t20\t0.15\t1\tLung\trecount3_GTEx_v8\tGTEx_v8\trecount3_GTEx_v8_GRCh38\n",
            )
            output = tmp / "merged.tsv"
            subprocess.run(
                [
                    sys.executable, str(SCRIPTS / "merge_normal_junction_tissues.py"),
                    "--inputs", str(liver), str(lung), "--output", str(output),
                ],
                check=True, capture_output=True, text=True,
            )
            row = read_tsv(output)[0]
            self.assertEqual(row["normal_samples"], "5")
            self.assertEqual(row["total_samples"], "30")
            self.assertEqual(row["normal_reads"], "8")
            self.assertEqual(row["normal_total_reads"], "19")
            self.assertEqual(row["normal_tissue_count"], "2")
            self.assertEqual(row["tissue"], "Liver,Lung")
            metadata = json.loads(Path(str(output) + ".meta.json").read_text(encoding="utf-8"))
            self.assertIn("sample-tissue records", metadata["denominator_semantics"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
