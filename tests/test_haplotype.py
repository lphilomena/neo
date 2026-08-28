from neoag.haplotype import annotate_nearby_variant_groups, combined_haplotype_peptides, read_bases_at_positions
from neoag.validation import make_validation_plan


def test_cigar_aware_base_extraction():
    fields = ["r1", "0", "chr2", "100", "60", "2S5M1I5M", "*", "0", "0", "NNACGTAACGTAA", "I" * 13]
    bases = read_bases_at_positions(fields, [100, 104, 105])
    assert bases == {100: "A", 104: "A", 105: "C"}


def test_nearby_snvs_require_phasing():
    events = [
        {"event_id": "e1", "chrom": "chr2", "pos": "100", "ref": "C", "alt": "A", "transcript_id": "tx"},
        {"event_id": "e2", "chrom": "chr2", "pos": "102", "ref": "A", "alt": "C", "transcript_id": "tx"},
    ]
    annotate_nearby_variant_groups(events)
    assert {row["haplotype_status"] for row in events} == {"PHASING_REQUIRED"}
    assert events[0]["phase_group_id"] == events[1]["phase_group_id"]


def test_combined_peptides_contain_both_adjacent_changes():
    rows = [
        {"variant_key": "e1", "minigene": "AAAA|Q|HBBBB", "protein_position": "5", "amino_acids": "H/Q"},
        {"variant_key": "e2", "minigene": "AAAAH|P|BBBB", "protein_position": "6", "amino_acids": "H/P"},
    ]
    ids, change, peptides = combined_haplotype_peptides(rows, peptide_lengths=(8,))
    assert ids == "e1;e2"
    assert "QP" in change.split(">", 1)[1]
    assert peptides
    assert all("QP" in row["peptide"] for row in peptides)


def test_validation_shortlist_keeps_only_two_per_phase_group():
    peptides = [
        {
            "peptide_id": f"p{i}", "event_id": "combined", "gene": "TBR1",
            "peptide": f"PEPTIDE{i}", "wildtype_peptide": f"WILDTYP{i}",
            "hla_allele": "HLA-A*02:06", "event_type": "SNV",
            "peptide_consequence": "phased_multi_substitution", "final_priority": "B",
            "safety_status": "PASS", "phase_group_id": "phase_tbr1",
            "redundancy_group": "phase_tbr1", "haplotype_status": "PHASED_CIS_COMBINED",
        }
        for i in range(3)
    ]
    rows = make_validation_plan(peptides)
    assert [row["shortlist_status"] for row in rows] == [
        "SHORTLISTED", "SHORTLISTED", "REDUNDANT_NOT_SHORTLISTED"
    ]
