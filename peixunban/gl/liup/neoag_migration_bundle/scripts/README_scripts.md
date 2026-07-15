# Migration install scripts

Run scripts from inside the unpacked migration bundle. Set `NEOAG_TOOLS_ROOT` to the target install root.

Recommended order:

```bash
export NEOAG_TOOLS_ROOT=/path/to/test_env_tools
bash scripts/install_tier1_core.sh
bash scripts/install_tier2_tools.sh
bash scripts/run_doctor_bundle_test.sh
```

Optional heavier layers:

```bash
bash scripts/install_gatk_vep_facets_ascat.sh
bash scripts/install_fusion_tools.sh
bash scripts/install_optional_immuno_hla_tools.sh
```

Licensed tools require user-provided packages:

```bash
NETMHCPAN_TARBALL=/path/to/netMHCpan.tar.gz \
NETMHCSTABPAN_TARBALL=/path/to/netMHCstabpan.tar.gz \
bash scripts/install_licensed_tools.sh
```

Reference staging template:

```bash
NEOAG_REF_ROOT=/path/to/neoag_refs bash scripts/stage_references_template.sh
REFERENCE_MANIFEST=refs/reference_manifest.local_template.yaml bash scripts/run_full_doctor.sh
```

Notes:

- `install_all_nonlicensed.sh` can download/install many large packages. Use it only on a machine intended for full deployment.
- VEP cache is not downloaded unless `INSTALL_VEP_CACHE=1` is set.
- NetMHCpan, NetMHCstabpan, Polysolver/Novoalign and similar tools must follow local license rules.



MHCflurry compatibility:

- MHCflurry 2.0.6 requires legacy Keras session APIs when TensorFlow 2.x installs Keras 3. The tier2 install script installs tf-keras, and common.sh exports TF_USE_LEGACY_KERAS=1 for reproducible CPU execution.
- The tier2 install script runs a one-peptide MHCflurry smoke test and writes logs/mhcflurry_smoke.out.csv.

BigMHC compatibility:

- BigMHC_IM is launched by neoag-v03 with the neoag-core Python interpreter. install_tier1_core.sh now installs CPU Python dependencies needed by BigMHC: torch, numpy, pandas, scipy, scikit-learn, and psutil.
- On offline machines, pre-populate pip/conda caches for these packages before running the tier1 script, or install them into the neoag-core environment manually.

Root one-shot install:

- On a new machine where /root/neo/env_tool is the target, run scripts/install_to_root_neo_env_tool.sh as root. It clones neo-na0707_upload_release and na0707_upload_release, extracts the migration and patch bundles, installs Miniforge under /root/neo/env_tool/miniforge3, and runs tier1/tier2 installers.
- For offline machines, set MINIFORGE_INSTALLER=/path/to/Miniforge3-Linux-x86_64.sh and pre-populate conda/pip caches. For real production runs, copy or mount the companion neodata4git directory under /root/neo/neodata4git, or set NEODATA_ROOT to its path.

Real VCF release test:

- After tier1/tier2 installation and companion reference setup, run scripts/run_real_vcf_pipeline_test.sh to execute the validated real VCF pipeline test.
- By default it uses the internal validation VCF at /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/data/liver_0520_WGS_shortReads/somatic/M1ML150017383_L01_438.align.somatic.pass.vcf.gz with HLA-A*02:06,HLA-A*30:01,HLA-B*13:02,HLA-C*06:02,HLA-C*08:01.
- On a new machine, either mount/copy that VCF path or override REAL_VCF, TUMOR_SAMPLE_NAME, SAMPLE_ID, HLA_ALLELES, and REAL_VCF_OUTDIR before running the script.
- The script checks ranked peptides/events, MHCflurry evidence, PRIME evidence, BigMHC evidence, and the HTML evidence report.
