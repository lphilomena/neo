from pathlib import Path

from neoag_v03.adapters.vcf_appm import (
    extract_appm_variants,
    extract_gx_expression,
    extract_appm_inputs_from_vcf,
)

ROOT = Path(__file__).resolve().parents[1]
HCC_VCF = ROOT / "data/examples/HCC1395/HCC1395_inputs/annotated.expression.vcf.gz"


def test_extract_gx_expression_hcc1395():
    if not HCC_VCF.is_file():
        return
    expr = extract_gx_expression(HCC_VCF, tumor_sample="HCC1395_TUMOR_DNA")
    assert expr["HLA-A"] > 100
    assert expr["NLRC5"] > 10


def test_extract_appm_variants_hcc1395():
    if not HCC_VCF.is_file():
        return
    variants = extract_appm_variants(HCC_VCF)
    assert "STAT1" in variants
    assert "missense_variant" in variants["STAT1"]


def test_extract_appm_inputs_from_vcf(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(
        "\n".join([
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR",
            'chr1\t100\t.\tC\tT\t.\tPASS\tCSQ=T|stop_gained|HIGH|B2M|ENSG00000166710|||||||||||||||||SNV||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\tGT:GX\t0/1:B2M|25.5',
        ]) + "\n",
        encoding="utf-8",
    )
    outs = extract_appm_inputs_from_vcf(vcf, tmp_path / "derived", tumor_sample="TUMOR")
    assert Path(outs["expression"]).is_file()
    assert Path(outs["vep_appm"]).is_file()
    text = Path(outs["expression"]).read_text(encoding="utf-8")
    assert "B2M" in text and "25.5" in text
    vep = Path(outs["vep_appm"]).read_text(encoding="utf-8")
    assert "B2M" in vep and "stop_gained" in vep


def test_extract_appm_inputs_include_damaging_missense_columns(tmp_path):
    vcf = tmp_path / "damaging_missense.vcf"
    vcf.write_text(
        "\n".join([
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SOURCE|GENE_PHENO|SIFT|PolyPhen">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR",
            'chr9\t100\t.\tC\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|JAK1|ENSG00000162434|Transcript|ENST00000342505||||c.1C>T|p.Ala1Val|||1|A/V|gCc/gTc||||SNV|||YES||||||||deleterious(0)|probably_damaging(0.99)\tGT:GX\t0/1:JAK1|12.0',
        ]) + "\n",
        encoding="utf-8",
    )
    outs = extract_appm_inputs_from_vcf(vcf, tmp_path / "derived", tumor_sample="TUMOR")
    text = Path(outs["vep_appm"]).read_text(encoding="utf-8")
    assert "is_damaging_missense" in text
    assert "JAK1" in text
    assert "damaging_missense" in text


def test_extract_appm_inputs_from_vep_tabular(tmp_path):
    vep = tmp_path / "vep_raw.tsv"
    vep.write_text(
        "\n".join([
            "## ENSEMBL VARIANT EFFECT PREDICTOR",
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\tSYMBOL\tHGVSc\tHGVSp\tSIFT\tPolyPhen",
            "var1\tchr15:450000\tA\tENSG00000138642\tENST0001\tTranscript\tsplice_donor_variant\tHIGH\tB2M\tc.1+1G>A\t-\t-\t-",
            "var2\tchr1:1\tT\tENSG0\tENST0\tTranscript\tmissense_variant\tMODERATE\tNOTAPPM\t-\t-\tdeleterious(0)\tprobably_damaging(1)",
        ]) + "\n",
        encoding="utf-8",
    )
    outs = extract_appm_inputs_from_vcf(vep, tmp_path / "derived")
    text = Path(outs["vep_appm"]).read_text(encoding="utf-8")
    assert "B2M" in text
    assert "splice_disrupting_variant" in text
    assert "NOTAPPM" not in text


def test_extract_appm_inputs_from_maf(tmp_path):
    maf = tmp_path / "mini.maf"
    maf.write_text(
        "\n".join([
            "Hugo_Symbol\tVariant_Classification\tIMPACT\tSIFT\tPolyPhen\tHGVSc\tHGVSp\tTumor_Seq_Allele2",
            "NLRC5\tFrame_Shift_Del\tHIGH\t\t\tc.10del\tp.G4fs\t-",
            "CIITA\tMissense_Mutation\tMODERATE\tdeleterious(0.01)\tprobably_damaging(0.9)\tc.3A>T\tp.K1N\tT",
        ]) + "\n",
        encoding="utf-8",
    )
    outs = extract_appm_inputs_from_vcf(maf, tmp_path / "derived")
    text = Path(outs["vep_appm"]).read_text(encoding="utf-8")
    assert "NLRC5" in text
    assert "CIITA" in text
    assert "damaging_missense" in text
    assert "protein_truncating_variant" in text


def test_vep_appm_covers_required_gene_set(tmp_path):
    vep = tmp_path / "empty_appm_vep.tsv"
    vep.write_text(
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\tSYMBOL\n"
        "var1\tchr1:1\tA\tENSG0\tENST0\tTranscript\tdownstream_gene_variant\tMODIFIER\tNOTAPPM\n",
        encoding="utf-8",
    )
    outs = extract_appm_inputs_from_vcf(vep, tmp_path / "derived")
    rows = Path(outs["vep_appm"]).read_text(encoding="utf-8").splitlines()
    required = {"B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP", "JAK1", "JAK2", "STAT1", "NLRC5", "CIITA", "RFX5", "RFXANK", "RFXAP"}
    text = "\n".join(rows)
    for gene in required:
        assert gene in text
    assert "assessed_no_variant" in text
