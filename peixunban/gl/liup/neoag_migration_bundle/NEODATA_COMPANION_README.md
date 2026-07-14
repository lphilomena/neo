# neodata4git companion references and licensed tools

This migration bundle is designed to install fresh conda/tool environments under:

/mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools

Real references and restricted/licensed tools are registered from the companion directory:

/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git

Use this companion directory only for authorized users/machines. Do not redistribute licensed tools publicly.

Key files added for production-style doctor checks:

- /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/activate_neoag_production_refs.sh
- /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/run_doctor_neodata4git.sh
- /mnt/zjl-bgi-zzb/peixunban/gl/liup/neoag_migration_bundle/refs/reference_manifest.neodata4git.yaml
- /mnt/zjl-bgi-zzb/peixunban/gl/liup/neoag_migration_bundle/configs/local/tools_manifest.neodata4git.yaml

Run validation:

source /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/activate_neoag_production_refs.sh
/mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/run_doctor_neodata4git.sh

Current expected doctor status: PARTIAL with no blocking issues. Remaining non-blocking items are optional/container tools, sample input placeholders, normal_ligandome, Arriba reference, and a VEP cache layout warning.

Update 2026-07-14:

- normal_ligandome is now staged as a real copied TSV:
  /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/refs/normal/normal_hla_ligands.hla_ligand_atlas.predicted_alleles.tsv
- arriba_reference is now staged as:
  /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/refs/arriba
  It links to real GRCh38 CTAT/STAR assets in neodata4git. Official optional Arriba databases were not found locally.
