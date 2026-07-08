from pathlib import Path
import csv

from neoag_v03.vep.extract_peptides import (
    MULTI_AA_COMPLEX,
    MULTI_AA_FRAMESHIFT,
    MULTI_AA_INFRAME,
    MULTI_AA_SINGLE,
    MULTI_AA_SUBSTITUTION,
    _parse_vaf_from_format,
    build_minigene,
    build_mutant_protein,
    build_proteome_kmer_index,
    classify_multi_aa_flag,
    extract_variant_peptides_from_vcf,
    parse_peptide_lengths,
    peptide_in_normal_proteome,
    sliding_full_mutant_mode,
    sliding_variant_peptides,
)

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
    "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
    "Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|"
    'TSL|HGVS_OFFSET|FrameshiftSequence|WildtypeProtein">'
)


def test_sliding_missense_excludes_wt_identical():
    wt = "ACDEFGHIKLMN"
    mut = "ACDENGHIKLMN"
    peps = list(sliding_variant_peptides(mut, wt, anchor_start=5, anchor_end=5, lengths=(8,)))
    assert any(p["mutant_peptide"] == "CDENGHIK" for p in peps)
    assert all(p["mutant_peptide"] != p["wildtype_peptide"] for p in peps if p["wildtype_peptide"])


def test_build_mutant_protein_missense():
    mut, a0, a1 = build_mutant_protein(
        "missense_variant",
        "MKLLVV",
        protein_position_raw="3",
        amino_acids="L/V",
        frameshift_sequence="",
    )
    assert mut == "MKVLVV"
    assert a0 == 3 and a1 == 3


def test_build_mutant_protein_multi_aa_substitution():
    wt = "ACDEFGHIKLMN"
    mut, a0, a1 = build_mutant_protein(
        "missense_variant",
        wt,
        protein_position_raw="4-5",
        amino_acids="EF/IS",
        frameshift_sequence="",
    )
    assert mut == "ACDISGHIKLMN"
    assert a0 == 4 and a1 >= 5


def test_classify_multi_aa_flag():
    assert classify_multi_aa_flag("missense_variant", "K/N", "10") == MULTI_AA_SINGLE
    assert classify_multi_aa_flag("missense_variant", "MG/IS", "1388-1389") == MULTI_AA_SUBSTITUTION


def test_parse_vaf_from_format():
    fmt = "GT:AD:AF:DP"
    tumor = "0/1:31,4:0.103:35"
    assert _parse_vaf_from_format(fmt, tumor) == "0.1030"
    assert _parse_vaf_from_format(fmt, "0/0:293,0:3.409e-03:293") == "0.0034"
    assert _parse_vaf_from_format("GT:AD", "0/1:20,5") == "0.2000"


def test_parse_allele_metrics_from_format():
    from neoag_v03.vep.extract_peptides import _parse_allele_metrics_from_format

    vaf, depth, alt = _parse_allele_metrics_from_format("GT:AD:AF:DP", "0/1:20,5:0.25:25")
    assert vaf == "0.2500"
    assert depth == "25"
    assert alt == "5"

    vaf2, depth2, alt2 = _parse_allele_metrics_from_format("GT:AD", "0/1:30,10")
    assert vaf2 == "0.2500"
    assert depth2 == "40"
    assert alt2 == "10"


def test_extract_variant_peptides_includes_dna_and_rna_metrics(tmp_path):
    vcf = tmp_path / "rna.vcf"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR_DNA\tTUMOR_RNA",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            "NM_0001.1|1|||ACDEFGHIKLMNPQRSTVWYACDEFG\tGT:AD:AF\t0/1:20,5:0.25\t0/1:40,12:0.2308",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "peptides.tsv"
    extract_variant_peptides_from_vcf(
        vcf,
        out,
        sample_id="T1",
        lengths=(8,),
        tumor_sample_name="TUMOR_DNA",
        rna_sample_name="TUMOR_RNA",
    )
    header, first_row = out.read_text(encoding="utf-8").strip().split("\n")[:2]
    cols = header.split("\t")
    cells = first_row.split("\t")
    assert cells[cols.index("vaf")] == "0.2500"
    assert cells[cols.index("tumor_depth")] == "25"
    assert cells[cols.index("tumor_alt_count")] == "5"
    assert cells[cols.index("rna_vaf")] == "0.2308"
    assert cells[cols.index("rna_alt_reads")] == "12"
    assert cells[cols.index("rna_depth")] == "52"


def test_extract_variant_peptides_auto_detects_rna_column(tmp_path):
    vcf = tmp_path / "rna_auto.vcf"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHCC1395_TUMOR_DNA\tHCC1395_TUMOR_RNA",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            "NM_0001.1|1|||ACDEFGHIKLMNPQRSTVWYACDEFG\tGT:AD\t0/1:20,5\t0/1:40,12",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "peptides.tsv"
    extract_variant_peptides_from_vcf(
        vcf,
        out,
        sample_id="T1",
        lengths=(8,),
        tumor_sample_name="HCC1395_TUMOR_DNA",
    )
    header, first_row = out.read_text(encoding="utf-8").strip().split("\n")[:2]
    cols = header.split("\t")
    cells = first_row.split("\t")
    assert cells[cols.index("rna_vaf")] == "0.2308"
    assert cells[cols.index("rna_alt_reads")] == "12"


def test_extract_variant_peptides_includes_vaf(tmp_path):
    vcf = tmp_path / "vaf.vcf"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            "NM_0001.1|1|||ACDEFGHIKLMNPQRSTVWYACDEFG\tGT:AD:AF\t0/1:20,5:0.25\t0/0:30,0:0",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "peptides.tsv"
    extract_variant_peptides_from_vcf(
        vcf, out, sample_id="T1", lengths=(8,), tumor_sample_name="TUMOR",
    )
    header, first_row = out.read_text(encoding="utf-8").strip().split("\n")[:2]
    assert "vaf" in header.split("\t")
    cols = header.split("\t")
    assert first_row.split("\t")[cols.index("vaf")] == "0.2500"


def test_extract_variant_peptides_from_mini_vcf(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            "NM_0001.1|1|||ACDEFGHIKLMNPQRSTVWYACDEFG",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "peptides.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,))
    assert summary["peptide_rows"] > 0
    text = out.read_text(encoding="utf-8")
    assert "GENE1" in text
    assert "multi_aa_flag" in text
    assert "peptide_source" in text
    assert MULTI_AA_SINGLE in text
    assert "\tsnv\t" in text


def test_build_minigene_missense():
    mg, mg_nt = build_minigene(
        "MKLLVV",
        "MKVLVV",
        anchor_start=3,
        anchor_end=3,
        amino_acids="L/V",
        mini_len=2,
    )
    assert mg == "MK|V|LV"
    assert mg_nt.count("|") == 2
    assert "GTA" in mg_nt


def test_parse_peptide_lengths_range():
    assert parse_peptide_lengths(length_min=8, length_max=11) == (8, 9, 10, 11)
    assert parse_peptide_lengths("9,11") == (9, 11)


def test_open_text_maybe_gz_plain_vcf_with_gz_suffix(tmp_path):
    from neoag_v03.utils import is_gzip_file, open_text_maybe_gz

    plain = tmp_path / "annot.vcf.gz"
    plain.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    assert not is_gzip_file(plain)
    with open_text_maybe_gz(plain) as fh:
        assert fh.readline().startswith("##fileformat")


def test_normal_proteome_filter(tmp_path):
    vcf = tmp_path / "mini.vcf"
    wt = "ACDEFGHIKLMNPQRSTVWYACDEFG"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            f"NM_0001.1|1|||{wt}",
        ]) + "\n",
        encoding="utf-8",
    )
    proteome = tmp_path / "normal.fa"
    proteome.write_text(">ref\nACDEFGHI\n", encoding="utf-8")
    out = tmp_path / "filtered.tsv"
    summary = extract_variant_peptides_from_vcf(
        vcf,
        out,
        sample_id="T1",
        lengths=(8,),
        normal_proteome_fasta=proteome,
        filter_normal_proteome=True,
    )
    index = build_proteome_kmer_index(["ACDEFGHI"], (8,))
    assert peptide_in_normal_proteome("ACDEFGHI", index)
    assert summary["peptides_filtered_normal_proteome"] >= 0


def test_exclude_multi_aa_filter(tmp_path):
    vcf = tmp_path / "multi.vcf"
    wt = "ACDEFGHIKLMNPQRSTVWYACDEFG"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t200\t.\tC\tG\t.\tPASS\tCSQ=G|missense_variant|MODERATE|GENE2|ENSG0002|"
            "Transcript|ENST00002.1|protein_coding|1/1||ENST00002.1:c.20C>G|"
            "ENSP00002.1:p.Glu6Gly|20|20|6-7|EF/IS|||1||HGNC|HGNC:2|YES||MANE_Select|"
            f"NM_0002.1|1|||{wt}",
        ]) + "\n",
        encoding="utf-8",
    )
    out_all = tmp_path / "all.tsv"
    out_single = tmp_path / "single.tsv"
    all_summary = extract_variant_peptides_from_vcf(vcf, out_all, sample_id="T2", lengths=(8,))
    single_summary = extract_variant_peptides_from_vcf(
        vcf, out_single, sample_id="T2", lengths=(8,), exclude_multi_aa=True,
    )
    assert all_summary["peptide_rows"] > 0
    assert single_summary["peptide_rows"] == 0
    assert single_summary["variants_skipped_multi_aa"] == 1


def test_sliding_full_mutant_mode_flags():
    assert sliding_full_mutant_mode("frameshift_variant", MULTI_AA_FRAMESHIFT)
    assert sliding_full_mutant_mode(
        "stop_gained&inframe_insertion",
        MULTI_AA_COMPLEX,
    )
    assert not sliding_full_mutant_mode("missense_variant", MULTI_AA_SINGLE)
    assert not sliding_full_mutant_mode(
        "missense_variant&splice_region_variant",
        MULTI_AA_COMPLEX,
    )


def test_rrm1_like_stop_inframe_insertion_includes_pvac_best_peptide(tmp_path):
    """RRM1 regression: pVAC Best Peptide YVTQDLNEV lies in the inserted novel tail."""
    wt = "M" + "L" * 177 + "KE" + "DIDAAIETYN"
    alt_insert = "-/HCSYVTQDLNEVVWE*EKGASRWTKS"
    consequence = "stop_gained&inframe_insertion"
    multi_aa = classify_multi_aa_flag(consequence, alt_insert, "180-181")
    assert multi_aa == MULTI_AA_COMPLEX
    assert sliding_full_mutant_mode(consequence, multi_aa)

    mut, anchor_start, anchor_end = build_mutant_protein(
        consequence,
        wt,
        protein_position_raw="180-181",
        amino_acids=alt_insert,
        frameshift_sequence="",
    )
    assert "YVTQDLNEV" in mut
    assert "EKGASRWTKS" not in mut

    narrow = {
        p["mutant_peptide"]
        for p in sliding_variant_peptides(
            mut,
            wt,
            anchor_start=anchor_start,
            anchor_end=anchor_end,
            lengths=(9,),
            frameshift_mode=False,
        )
    }
    full = {
        p["mutant_peptide"]
        for p in sliding_variant_peptides(
            mut,
            wt,
            anchor_start=anchor_start,
            anchor_end=anchor_end,
            lengths=(9,),
            frameshift_mode=True,
        )
    }
    assert "YVTQDLNEV" not in narrow
    assert "YVTQDLNEV" in full

    vcf = tmp_path / "rrm1.vcf"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr11\t4111952\t.\tA\tACACTGCTCATATGTGACCCAAGATCTAAATGAAGTAGTTTGGGAGTGAGAGAAGGGAGCGAGCAGATGGACCAAAAGC\t.\tPASS\t"
            f"CSQ=ACACTGCTCATATGTGACCCAAGATCTAAATGAAGTAGTTTGGGAGTGAGAGAAGGGAGCGAGCAGATGGACCAAAAGC|"
            f"{consequence}|HIGH|RRM1|ENSG00000167325|Transcript|ENST00000300738.10|protein_coding|"
            "12/13||ENST00000300738.10:c.540_541insCACTGCTCATATGTGACCCAAGATCTAAATGAAGTAGTTTGGGAGTGAGAGAAGGGAGCGAGCAGATGGACCAAAAGC|"
            "ENSP00000300738.5:p.Lys180_Glu181insHisCysSerTyrValThrGlnAspLeuAsnGluValValTrpGluTer|"
            "754|540-541|180-181|"
            f"{alt_insert}|"
            "-/CACTGCTCATATGTGACCCAAGATCTAAATGAAGTAGTTTGGGAGTGAGAGAAGGGAGCGAGCAGATGGACCAAAAGC|"
            "|||1||HGNC:10451|YES||MANE_Select|1|||"
            f"{wt}",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "rrm1_peptides.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="RRM1", lengths=(9,))
    with out.open(encoding="utf-8") as fh:
        peptides = {row["mutant_peptide"] for row in csv.DictReader(fh, delimiter="\t")}
    assert summary["peptide_rows"] > 0
    assert "YVTQDLNEV" in peptides



def test_build_minigene_indel_can_center_on_short_peptide():
    mg, mg_nt = build_minigene(
        "ACDEFGHIKLMNPQRST",
        "ACDEYWVFGHIKLMNPQRST",
        anchor_start=5,
        anchor_end=6,
        amino_acids="-/YWV",
        mini_len=2,
        peptide_start_aa=3,
        peptide_end_aa=10,
        peptide_centered=True,
    )
    assert mg == "AC|DEYWVFGH|IK"
    assert mg.split("|")[1] == "DEYWVFGH"
    assert mg_nt.count("|") == 2


def test_extract_indel_minigene_is_centered_on_each_mutant_peptide(tmp_path):
    wt = "ACDEFGHIKLMNPQRST"
    vcf = tmp_path / "indel.vcf"
    vcf.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tATGGTACTGG\t.\tPASS\t"
            "CSQ=ATGGTACTGG|inframe_insertion|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.15_16insTGGTACTGG|"
            "ENSP00001.1:p.Glu5_Phe6insTyrTrpVal|15|15-16|5-6|-/YWV|"
            "-/TGGTACTGG|||1||HGNC|HGNC:1|YES||MANE_Select|NM_0001.1|||"
            f"{wt}",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "indel_peptides.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,), mini_len=2)
    assert summary["peptide_rows"] > 0
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["multi_aa_flag"] for r in rows} == {MULTI_AA_INFRAME}
    assert all(r["minigene"].split("|")[1] == r["mutant_peptide"] for r in rows)
    assert any(r["minigene"] == "AC|DEYWVFGH|IK" for r in rows)
