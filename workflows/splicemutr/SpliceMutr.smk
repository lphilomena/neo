"""Production SpliceMutr cohort workflow for NeoAg.

Run with: splicemutr-neoag workflow workflows/splicemutr/SpliceMutr.smk
          --configfile configs/workflows/splicemutr.cohort.yaml --cores 16
"""

from pathlib import Path

ROOT = Path(workflow.basedir).resolve().parents[1]
OUT = Path(config["output_dir"]).resolve()
COHORT = OUT / "cohort"
LEAF = OUT / "leafcutter"
FORMED = OUT / "formed_transcripts"
COMBINED = OUT / "combined"
HOME = Path(config["splicemutr_home"]).resolve()
R = config.get("rscript", "Rscript")
PYTHON = config.get("python", "python")
SAMPLE_ID = config["target_sample_id"]


rule all:
    input:
        OUT / "splicemutr_candidates.tsv",
        OUT / "splicemutr_candidates.manifest.json",
        COHORT / "cohort_summary.json",
        COMBINED / "data_splicemutr_all_pep.txt"


rule prepare_cohort:
    input:
        samples=config["sample_sheet"]
    output:
        normalized=COHORT / "cohort_samples.normalized.tsv",
        groups=COHORT / "groups_file.txt",
        juncfiles=COHORT / "juncfiles.txt",
        summary=COHORT / "cohort_summary.json"
    params:
        min_normal=config.get("min_normal_samples", 2),
        min_tumor=config.get("min_tumor_samples", 1),
        min_reads=config.get("min_unique_reads", 10),
        low_power="--allow-low-power" if config.get("allow_low_power", False) else ""
    shell:
        """
        {PYTHON:q} {ROOT}/scripts/prepare_splicemutr_cohort.py \
          --samples {input.samples:q} --outdir {COHORT:q} \
          --min-normal-samples {params.min_normal} --min-tumor-samples {params.min_tumor} \
          --min-unique-reads {params.min_reads} {params.low_power}
        """


rule leafcutter_differential:
    input:
        groups=rules.prepare_cohort.output.groups,
        juncfiles=rules.prepare_cohort.output.juncfiles
    output:
        counts=LEAF / "cohort_perind_numers.counts.gz",
        significance=LEAF / "leafcutter_ds_cluster_significance.txt",
        effects=LEAF / "leafcutter_ds_effect_sizes.txt",
        rdata=LEAF / "data.Rdata",
        introns=LEAF / "splicemutr_introns.rds"
    params:
        cluster=config["leafcutter_cluster_script"],
        ds=config["leafcutter_ds_script"],
        prepare=config["leafviz_prepare_results_script"],
        exon=config["leafcutter_exon_file"],
        annotation=config["leafcutter_annotation_prefix"],
        max_intron=config.get("max_intron_length", 500000),
        min_samples=config.get("leafcutter_min_samples", 3),
        threads=config.get("threads", 8)
    shell:
        """
        mkdir -p {LEAF:q}
        {PYTHON:q} {params.cluster:q} -j {input.juncfiles:q} -r {LEAF:q} -o cohort -l {params.max_intron}
        {params.ds:q} -i {params.min_samples} --num_threads {params.threads} \
          --exon_file={params.exon:q} -o {LEAF}/leafcutter_ds {output.counts:q} {input.groups:q}
        {params.prepare:q} -o {output.rdata:q} -m {input.groups:q} {output.counts:q} \
          {output.significance:q} {output.effects:q} {params.annotation:q}
        {R:q} {HOME}/Rscripts/save_introns.R -i {output.rdata:q} -o splicemutr
        """


rule form_transcripts:
    input:
        introns=rules.leafcutter_differential.output.introns,
        txdb=config["txdb"]
    output:
        metadata=FORMED / "splicemutr_introns_data_splicemutr.rds",
        sequences=FORMED / "splicemutr_introns_sequences.fa"
    params:
        bsgenome=config["bsgenome_name"]
    shell:
        """
        mkdir -p {FORMED:q}
        {R:q} {HOME}/Rscripts/form_transcripts.R -o {FORMED}/splicemutr_introns \
          -t {input.txdb:q} -j {input.introns:q} -b {params.bsgenome:q} -f {HOME}/Rfunctions/functions.R
        """


rule coding_potential:
    input:
        metadata=rules.form_transcripts.output.metadata,
        sequences=rules.form_transcripts.output.sequences
    output:
        corrected=FORMED / "splicemutr_introns_data_splicemutr_cp_corrected.rds",
        filelist=FORMED / "filenames_cp.txt"
    shell:
        """
        {R:q} {HOME}/Rscripts/calc_coding_potential.R -s {input.metadata:q} -t {input.sequences:q} \
          -o {FORMED:q} -f {HOME}/Rfunctions/functions.R
        printf '%s\n' {output.corrected:q} > {output.filelist:q}
        """


rule combine_splicemutr:
    input:
        files=rules.coding_potential.output.filelist
    output:
        metadata=COMBINED / "data_splicemutr_all_pep.txt",
        proteins=COMBINED / "proteins.txt"
    shell:
        """
        mkdir -p {COMBINED:q}
        {R:q} {HOME}/Rscripts/combine_splicemutr.R -o {COMBINED:q} -s {input.files:q}
        """


rule export_candidates:
    input:
        metadata=rules.combine_splicemutr.output.metadata
    output:
        candidates=OUT / "splicemutr_candidates.tsv",
        manifest=OUT / "splicemutr_candidates.manifest.json"
    params:
        lengths=config.get("peptide_lengths", "8,9,10,11"),
        normal_arg=("--normal-junctions " + config["normal_junction_reference"])
                   if config.get("normal_junction_reference") else ""
    shell:
        """
        {PYTHON:q} {ROOT}/scripts/export_splicemutr_candidates.py --metadata {input.metadata:q} \
          --sample-id {SAMPLE_ID:q} --out {output.candidates:q} --peptide-lengths {params.lengths:q} \
          {params.normal_arg}
        """
