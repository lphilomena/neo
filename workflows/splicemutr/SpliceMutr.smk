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
LEAF_R = config.get("leafcutter_rscript", R)
PYTHON = config.get("python", "python")
SAMPLE_ID = config["target_sample_id"]
ANALYSIS_MODE = config.get("analysis_mode", "differential").lower()

if ANALYSIS_MODE not in {"differential", "outlier"}:
    raise ValueError("analysis_mode must be differential or outlier")


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


if ANALYSIS_MODE == "outlier":
    rule leafcutter_differential:
        input:
            groups=rules.prepare_cohort.output.groups,
            juncfiles=rules.prepare_cohort.output.juncfiles
        output:
            counts=LEAF / "cohort_perind_numers.counts.gz",
            pvalues=LEAF / "cohort_pVals.txt",
            introns=LEAF / "splicemutr_introns.rds"
        params:
            cluster=config["leafcutter_cluster_script"],
            outlier=config["leafcutter_outlier_script"],
            max_intron=config.get("max_intron_length", 500000),
            min_coverage=config.get("leafcutter_outlier_min_coverage", 20),
            threshold=config.get("leafcutter_outlier_pvalue", 0.05),
            adjust=config.get("leafcutter_outlier_adjust_method", "BH"),
            threads=config.get("threads", 8)
        shell:
            """
            mkdir -p {LEAF:q}
            {PYTHON:q} {params.cluster:q} -j {input.juncfiles:q} -r {LEAF:q} -o cohort -l {params.max_intron}
            {LEAF_R:q} {params.outlier:q} -p {params.threads} -c {params.min_coverage} \
              -o {LEAF}/cohort {output.counts:q}
            {R:q} {ROOT}/scripts/extract_leafcutter_outlier_introns.R \
              --pvalues {output.pvalues:q} --sample {SAMPLE_ID:q} \
              --threshold {params.threshold} --adjust-method {params.adjust:q} --out {output.introns:q}
            """
else:
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
            {LEAF_R:q} {params.ds:q} -i {params.min_samples} --num_threads {params.threads} \
              --exon_file={params.exon:q} -o {LEAF}/leafcutter_ds {output.counts:q} {input.groups:q}
            {LEAF_R:q} {params.prepare:q} -o {output.rdata:q} -m {input.groups:q} {output.counts:q} \
              {output.significance:q} {output.effects:q} {params.annotation:q}
            {R:q} {HOME}/Rscripts/save_introns.R -i {output.rdata:q} -o splicemutr
            """


rule form_transcripts:
    input:
        introns=rules.leafcutter_differential.output.introns,
        txdb=config["txdb"]
    output:
        metadata_files=FORMED / "formed_metadata_files.txt",
        done=FORMED / "form_transcripts.done"
    params:
        bsgenome=config["bsgenome_name"],
        bsgenome_package=config.get("bsgenome_package", "BSgenome.Hsapiens.UCSC.hg38"),
        chunks=config.get("form_transcript_chunks", min(config.get("threads", 8), 8)),
        max_parallel=config.get("form_transcript_max_parallel", min(config.get("threads", 8), 4))
    shell:
        """
        {ROOT}/scripts/run_splicemutr_form_parallel.sh \
          --rscript {R:q} --form-script {HOME}/Rscripts/form_transcripts.R \
          --split-script {ROOT}/scripts/split_splicemutr_introns.R \
          --introns {input.introns:q} --txdb {input.txdb:q} \
          --bsgenome-object {params.bsgenome:q} --bsgenome-package {params.bsgenome_package:q} \
          --functions {HOME}/Rfunctions/functions.R --outdir {FORMED:q} --chunks {params.chunks} \
          --max-parallel {params.max_parallel}
        """


rule coding_potential:
    input:
        metadata_files=rules.form_transcripts.output.metadata_files,
        done=rules.form_transcripts.output.done
    output:
        filelist=FORMED / "filenames_cp.txt"
    shell:
        """
        while IFS= read -r metadata; do
          prefix="${{metadata%_data_splicemutr.rds}}"
          {R:q} {HOME}/Rscripts/calc_coding_potential.R -s "$metadata" \
            -t "${{prefix}}_sequences.fa" -o {FORMED:q} -f {HOME}/Rfunctions/functions.R
        done < {input.metadata_files:q}
        find {FORMED:q} -maxdepth 1 -type f -name '*_cp_corrected.rds' -print | sort > {output.filelist:q}
        test -s {output.filelist:q}
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
