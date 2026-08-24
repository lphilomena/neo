---
name: open-neo-install-check
description: Public macro Skill1 for Open-Neo installation, new-machine migration, environment/reference/tool validation, Doctor and smoke-test readiness. Use it before running analysis on a new or changed machine.
---

# Open-Neo Install Check

## Use when

- Deploying or migrating Open-Neo to another machine.
- Determining whether the machine is ready for `review`, `core`, `prediction`, or `full` use.
- Checking tools, models, caches, references, licensed assets, release boundaries, and minimal smoke tests.

## Do not use when

- The user only wants to interpret existing results; use `open-neo-review`.
- The environment is already verified and the user wants to run a case; use `open-neo-run`.

## Modes

- `plan`: collect inventory and prepare local manifest templates without executing smoke commands.
- `verify`: run read-only Doctor and selected smoke tests.
- `repair` / `install`: require explicit human approval; the macro Skill itself never bypasses license terms.
- `resume`: rerun an interrupted idempotent deployment, or reuse a matching PASS checkpoint; approval remains required.

## Required input

- `project_root`, or a verified release tarball supplied by the user.
- Optional `--conda-base DIR` / `NEOAG_CONDA_BASE`: an existing
  site-managed Conda or Miniforge root, including a NAS path. When supplied,
  the Skill uses that installation and does not search for or install another
  Miniconda/Miniforge.
- Optional `--install-claude-code`: install Claude Code with Anthropic's
  official native installer. The default release channel is `stable`; use
  `--claude-code-channel latest` or an exact `X.Y.Z` version when required.
  Actual installation requires `--approved --allow-download` and never performs
  authentication or stores an API key.

## Procedure

1. Require and verify the checksum for a release archive, reject traversal/link/device members, safely stage it under the output directory, and identify one project root.
2. Record Python, Java, Docker/Apptainer, Nextflow, disk, and platform information.
3. Generate comprehensive machine-local tools_manifest, reference_manifest, paths.env, and production_assets.local.tsv. The repository's configs/references/reference_manifest.yaml is authoritative; the local copy preserves its required, marker, sha256, version, and optional/required semantics while only remapping machine roots. Generic /srv targets are rewritten to the selected machine roots.
4. Automatically discover existing executables and references from the input manifests, documented environment variables, PATH, conda-style roots, and the standard portable data layout. Verify build-sensitive references conservatively; for example, a GRCh37 BAM-matcher SNP panel is rejected for GRCh38.
5. Generate configured manifests and validated command templates. A tool is `PARTIAL` when its executable exists but a required reference or confirmed sample-level invocation is missing; the Skill never invents cohort-specific SNAF/SpliceMutr workflows.
6. For approved repair/install, delegate to the portable neoag-remote-deploy new-machine installer; downloads remain opt-in and licensed assets remain external. Directory references are not accepted merely because the directory exists: their declared marker (for example versionInfo.json, Genome, ref_genome.fa.star.idx, or a SpecHLA database marker) must also exist. For `prediction`/`full`, the default installer profile is `--all-open`, sized so a successful Skill1 install can feed **`open-neo-run`** and multi-tool consensus without side-path installs. PRIME, MixMHCpred, BigMHC and DeepImmuno are installed and checked as advisory immunogenicity support when assets are available; missing advisory support must warn but not block READY. It includes the standard production path (core-env, VEP, GATK, immunogenicity, OptiType, FACETS, splice, LOHHLA, **fusion** with Nextflow / STAR / EasyFuse family, **BAM-matcher**, and CPU **torch**) plus open, redistributable alternative tools where available. `--standard` remains available as a lighter production-main-path install, and `--minimal` remains the review/core install.
7. Re-run discovery after installation, pin production runtime entrypoints, execute minimal smoke checks, and publish generated machine-local manifests under `configs/local/`. The runtime pinning records the local FACETS R environment and `snp-pileup`, Sequenza chr FASTA/GC wiggle, HMFTOOLS/PURPLE suite paths, SpecHLA `samtools`, NetMHCstabpan and NetChop 3.1d runtime paths, and a project-local `bedtools` wrapper for LOHHLA when the portable bedtools package is present. For fusion callers, Skill1 must ensure EasyFuse v2 installations with legacy module references have `environments/easyfuse_src.yml` and `environments/requantification_wo_easyfuse.yml`, resolve EasyFuse to the deployed `open-neo-deploy/env_tool/tools/EasyFuse` path when the project-local tool path is absent, install `easy-fuse` shims into STAR/Bowtie-only Nextflow conda envs, keep EasyFuse STAR_INDEX and STAR_CUSTOM STAR versions compatible, remove STAR-version-specific transient-index parameters such as `genomeType` and `genomeTransform*` both after STAR_INDEX creation and immediately before STAR_CUSTOM reads the transient index, pin `NXF_CONDA_CACHEDIR` to an absolute project work cache, avoid conflicting conda/mamba yes flags, pin EasyFuse STAR alignment to a readable `starfusion_index/ref_genome.fa.star.idx` containing `Genome` and `genomeParameters.txt`, patch EasyFuse FusionCatcher environments to use STAR 2.5.2b, tolerate the pinned FusionCatcher 1.33 reference build when only the cached 1.00 package is available, and build or repair a missing/partial Bowtie1 genome_index from the bundled Bowtie2 genome_index2, and provide the lincrnas.txt alias when only lncrnas.txt is present plus legacy optional FusionCatcher filter aliases/empty files, expose STAR-Fusion through a Perl runtime compatible with `Set::IntervalTree` and tolerate already-patched STAR-Fusion version-check scripts, include Perl `vendor_perl` paths for modules such as `common::sense`, resolve the deployed STAR-Fusion executable when it is not on PATH, repair required STAR-Fusion sidecar components `FusionFilter` and `FusionAnnotator`, expose the STAR executable used by STAR-Fusion, and configure Arriba to use fixed blacklist assets plus Arriba 2.5.1 argument names before production runs. Licensed tools are configured only when the user has supplied a legal local installation.
8. Run Doctor and evaluate required tools, alternative capability groups, critical references and sidecars against the requested tier. Missing required references can never produce READY.
9. Write deployment checkpoints, configuration fixes, and an installation status delta for safe resume/review. Full installs default to no wall-clock timeout; `--install-timeout SECONDS` can be supplied when an operator wants a bounded run. If the macro is interrupted or times out, it terminates the whole installer process group before recording an `INTERRUPTED` or `TIMEOUT` checkpoint so `resume` can continue from a controlled state.
10. Write `deployment_report.md`, machine-readable status files, and an audit log.

For PRIME, MixMHCpred and BigMHC, the approved installer treats the pinned
source tree and model files as deployment assets. It first synchronizes
`data/predictors/{prime,mixMHCpred_install,bigmhc}` from the selected asset
source and verifies tool-specific markers plus the pinned Git revision when
revision metadata is present. These directories remain outside GitHub. A
network snapshot download is only a fallback when the asset is unavailable;
downloaders detect old curl builds that lack `--retry-all-errors` and retain
portable retry/timeout behavior instead of failing on the unknown option.
This compatibility policy also applies to the Bioconductor data-cache helper
used by ASCAT/Sequenza, LOHHLA source retrieval, SpliceMutr assets, SNAF source snapshots, and normal
junction asset construction. Ubuntu 20.04's curl 7.68 is therefore supported;
an unavailable `--retry-all-errors` capability is not a machine-configuration
failure. The immunogenicity source downloader can additionally fall back to
`wget`.

For asset synchronization, approved `install`, `repair`, and `resume` runs must
use an explicit asset source. Set `--asset-source-root /mounted/assets` for a
local/mounted asset tree, or set both `--asset-source-host user@host` and
`--asset-source-root /path/on/source/host` for remote rsync. The repository
manifest uses the portable placeholder `/srv/neoag-assets/source`; Skill1 rewrites
that placeholder only after an explicit source root is supplied. It no longer
defaults to a lab-specific host or username. Approved execution records source
reachability and stops before installation when a required asset source cannot
be read. Optional assets are reported as missing/failed optional assets and do
not masquerade as installed references. Use `--asset-ssh-key ~/.ssh/id_ed25519`
or `OPEN_NEO_ASSET_SSH_KEY` when the asset host requires a specific key. The key
is passed only to the rsync asset sync step and is recorded in `manifests/paths.env`
without copying the key.

Licensed presentation assets should be staged under the selected licensed root
or the selected asset source before installation. For NetChop this means both
`netchop-3.1d.Linux.tar.gz` and the institutional license file are expected in
a readable `netchop/` directory, or `--netchop-archive` must point to the
authorized archive. SNAF source snapshot downloads are network fallback only;
Skill1 verifies the downloaded tarball and retries before failing. SpliceMutr hg38 BSgenome is treated as a synchronizable open asset; if it is absent from the selected asset source, approved runs with `--allow-download` download and stage the Bioconductor package before reference verification.

For a NAS-managed Conda installation, use:

```bash
open-neo install-check \
  --project-root . \
  --mode verify \
  --conda-base /nas/path/to/miniforge3 \
  --outdir work/install-check
```

The selected path is propagated to the portable installer as `--conda-base`
and recorded in `manifests/paths.env`.

To include Claude Code in an approved new-machine installation:

```bash
open-neo install-check \
  --project-root . \
  --mode install \
  --install-claude-code \
  --claude-code-channel stable \
  --approved \
  --allow-download \
  --outdir work/install-check
```

In `plan` mode these options are recorded in `deployment_command.json` without
downloading anything. Installation verifies `claude --version`; login remains a
separate interactive or enterprise-managed step.

## Reproducible derived assets

- Java is installed in the dedicated `neoag-runtime` environment; discover it through the configured tools manifest instead of assuming it is on the login-shell PATH.
- Normal expression background `data/normal/expression/normal_expression.gtex_v11_hpa_hspc.tsv` is a required prediction/full asset and must be synchronized to the machine-local fixed asset directory, exported as `NEOAG_NORMAL_EXPRESSION`, and checked during readiness; do not rely on the source NAS path at run time. Salmon index and tx2gene.tsv must come from the same GENCODE release. The portable v49 layout is data/rna/gencode_v49/{salmon_index,tx2gene.tsv}; salmon_index/versionInfo.json is mandatory and the local manifest uses the canonical key salmon_tx2gene.
- Install BAM-matcher with `scripts/install_bam_matcher.sh`. It pins the upstream revision and retains its isolated Python 2 environment; do not merge these legacy dependencies into the main project environment.
- Never use BAM-matcher's bundled GRCh37 loci with GRCh38 data. Build the portable panel with `scripts/build_bam_matcher_grch38_loci.sh`; the script requires exact dbSNP-ID mapping, biallelic SNVs, GRCh38 FASTA REF agreement, and emits a metadata manifest plus checksums.
- PRIME, MixMHCpred and BigMHC use pinned complete source assets, not model-only or wrapper-only copies. Required markers are `lib/run_PRIME.pl`, `MixMHCpred`, and `src/predict.py`; BigMHC must also contain `models/`. Resume reuses these synchronized directories instead of downloading them again.
- Full installs must leave FACETS, Sequenza, PURPLE/HMFTOOLS, SpecHLA, LOHHLA, and fusion callers callable without relying on interactive shell state. LOHHLA readiness also requires a readable licensed Polysolver tree, project-local BAM index/link handling, absolute LOHHLA inputs, fresh explicit LOHHLA work directories for reruns, and preservation/regeneration of root mpileup intermediates instead of reusing partial failed directories. `conf/tools.env.local.sh` receives a marked install-check block for these paths, and `bin/bedtools` is generated when the portable bedtools binary requires the bundled Conda C++ runtime. EasyFuse readiness requires its reference marker, Nextflow entrypoint, compatibility env YAML files when referenced by installed modules, deployed tool-directory fallback, FusionCatcher STAR 2.5.2b/build-version compatibility patching plus Bowtie1 genome_index complete-file validation/repair and lincrnas.txt aliasing plus legacy optional filter repair, an absolute Nextflow conda cache, and an explicit STAR index override to `starfusion_index/ref_genome.fa.star.idx` with both `Genome` and `genomeParameters.txt` readable. STAR-Fusion readiness requires the CTAT/EasyFuse starfusion index, deployed executable fallback, compatible Perl invocation for `Set::IntervalTree`, `vendor_perl` in `PERL5LIB` for modules such as `common::sense`, complete `FusionFilter`/`FusionAnnotator` sidecars with the `annot_filter.pass` filename compatibility patch, a STAR executable on PATH, and conservative default threads for reruns. Arriba readiness requires fixed blacklist assets, fixed reference lookup, and correct Arriba 2.5.1 argument mapping (`-O` for discarded output, `-p` for protein domains). NetMHCstabpan readiness must install/register a callable licensed/shim executable from the deployed licensed tools root or shared asset root, preserve the supplied license file when present, and provide a `gawk`-compatible runtime entry. NetMHCpan readiness must validate a real peptide-HLA prediction and prefer the container wrapper when native NetMHCpan 4.2c is incompatible with the host glibc. NetChop readiness must install from an authorized `netchop-3.1d.Linux.tar.gz` discovered in the licensed root or shared asset root, keep the license file alongside it when supplied, and pin `NETCHOP_HOME`/`NEOAG_NETCHOP_BIN` into `conf/tools.env.local.sh`. Re-running `install` or `resume` replaces only the marked block/wrapper and skips already synchronized assets when their checkpoints and markers are still valid. For prediction/full tiers, NetMHCpan, MHCflurry, NetMHCstabpan, and NetChop are all hard readiness requirements because raw heavy production manifests require all four presentation predictors.

## Outputs

- `environment_inventory.tsv`
- `doctor/doctor_status.json`
- `deployment_status.tsv`
- `claude_code_status.tsv` and `claude_code_status.json`, including requested
  state, readiness, version, binary path, and the nested install report path
- `tier_requirements.tsv`, `deployment_delta.tsv`
- `production_run_readiness.tsv`, covering full-run smoke compatibility for
  Sequenza, PURPLE/HMFTOOLS, FACETS, STAR/EasyFuse, BAM-matcher, HLA-LA, and
  optional NeoAg Gateway health
- `deployment_checkpoint.json` for mutating modes
- `deployment_report.md`
- `manifests/tools_manifest.local.yaml`
- `manifests/reference_manifest.local.yaml`
- `auto_configuration*/manifests/tools_manifest.configured.yaml`
- `auto_configuration*/manifests/reference_manifest.configured.yaml`
- `auto_configuration*/manifests/command_templates.yaml`
- `auto_configuration*/configuration_status.tsv`
- `auto_configuration*/smoke_tests.tsv`
- `auto_configuration*/recommended_fixes.md`
- `configs/local/*.generated.yaml` after an approved successful installation
- `skill_result.json`
- `run_state.json`

## Safety boundary

Do not install or redistribute licensed tools, download large references, overwrite production settings, delete files, or submit HPC jobs unless the user explicitly approves the operation.

## Contracts and failure handling

- Validate public inputs against `references/INPUT_SCHEMA.json` before execution.
- Emit the stable result contract described by `references/OUTPUT_SCHEMA.json`.
- Use the canonical failure codes and remediation in `references/FAILURE_CODES.md`.
- Every invocation writes `skill_result.json` and a sibling `run_state.json`; install and repair actions remain approval gated.
- The launcher must resolve Python 3.11+ with `tomllib` from `NEOAG_PYTHON`,
  configured Conda/deployment roots, the project virtual environment, or PATH;
  it must fail clearly instead of starting the macro with an older interpreter.
