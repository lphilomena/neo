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
- Optional `--conda-pkgs-source DIR` / `NEOAG_CONDA_PKGS_SOURCE`: a readable,
  pre-populated Conda package cache from the asset source or a previous machine.
  Skill1 copies missing cache entries into the machine-local writable
  `env_tool/conda_pkgs` before solving environments.
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
6. For approved repair/install, delegate to the portable neoag-remote-deploy new-machine installer; downloads remain opt-in and licensed assets remain external. Directory references are not accepted merely because the directory exists: their declared marker (for example versionInfo.json, Genome, ref_genome.fa.star.idx, or a SpecHLA database marker) must also exist. For `prediction`/`full`, the default installer profile is `--all-open`, sized so a successful Skill1 install can feed **`open-neo-run`** and multi-tool consensus without side-path installs. PRIME, MixMHCpred, BigMHC and DeepImmuno are installed and checked as advisory immunogenicity support when assets are available; missing advisory support must warn but not block READY. It includes the standard production path (core-env, VEP, GATK, immunogenicity, OptiType, FACETS, splice, LOHHLA, **fusion through one EasyFuse-native installation**, the complete **DNA-SV** group Manta / SvABA / GRIDSS, **BAM-matcher**, and CPU **torch**) plus open, redistributable alternative tools where available. EasyFuse owns its STAR, Arriba, STAR-Fusion and FusionCatcher module environments; Skill1 must not create independent duplicate caller installations by default. The splice group installs SNAF, the SNAF-author AltAnalyze 0.7.0.1 image pinned by immutable digest, and the pinned SpliceMutr workflow by default; each must pass its own readiness check. `--standard` remains available as a lighter production-main-path install, and `--minimal` remains the review/core install.

For full readiness, validate the resolved VEP cache layout rather than accepting a directory name: the configured root must contain `homo_sapiens/<version>_GRCh38`, and an empty or stale `~/.vep` must not override the deployed cache. BAM Matcher readiness requires a working Java runtime whose absolute path is written into the generated matcher configuration, plus matching chromosome naming across BAM/reference FASTA/loci VCF. PRIME/MixMHCpred readiness must probe the real installed entrypoints and import NumPy with the exact interpreter pinned by `NEOAG_PRIME_PYTHON`; auxiliary predictor failures remain explicit warnings rather than false READY results.
Treat the resolved Conda environment prefix as authoritative even when Skill1 is launched from an outer `.venv`: environment setup and smoke tests must invoke `${ENV_PREFIX}/bin/python`, `${ENV_PREFIX}/bin/python -m pip`, and target-environment executables directly rather than trusting `conda activate`, bare `python`, bare `pip`, or the inherited `PATH`. VEP follows the same isolation rule: invoke `${VEP_ENV_PREFIX}/bin/perl`, `${VEP_ENV_PREFIX}/bin/vep`, and `${VEP_ENV_PREFIX}/bin/vep_install` with conflicting `PERL5LIB`/local-lib variables removed, and require a real `perl -MDBI` plus `vep --help` smoke test before recording READY. A DBI module found only in another environment does not satisfy VEP readiness.
7. Re-run discovery after installation, pin production runtime entrypoints, execute minimal smoke checks, and publish generated machine-local manifests under `configs/local/`. The runtime pinning records the local FACETS R environment and `snp-pileup`, Sequenza chr FASTA/GC wiggle, HMFTOOLS/PURPLE suite paths, SpecHLA `samtools`, NetMHCstabpan and NetChop 3.1d runtime paths, and a project-local `bedtools` wrapper for LOHHLA when the portable bedtools package is present. For fusion callers, Skill1 must ensure EasyFuse v2 installations with legacy module references have `environments/easyfuse_src.yml` and `environments/requantification_wo_easyfuse.yml`, resolve EasyFuse to the deployed `open-neo-deploy/env_tool/tools/EasyFuse` path when the project-local tool path is absent, install `easy-fuse` shims into STAR/Bowtie-only Nextflow conda envs, keep EasyFuse STAR_INDEX and STAR_CUSTOM STAR versions compatible, remove STAR-version-specific transient-index parameters such as `genomeType` and `genomeTransform*` both after STAR_INDEX creation and immediately before STAR_CUSTOM reads the transient index, pin `NXF_CONDA_CACHEDIR` to an absolute project work cache, avoid conflicting conda/mamba yes flags, pin EasyFuse STAR alignment to a readable `starfusion_index/ref_genome.fa.star.idx` containing `Genome` and `genomeParameters.txt`, patch EasyFuse FusionCatcher environments to use STAR 2.5.2b, tolerate the pinned FusionCatcher 1.33 reference build when only the cached 1.00 package is available, and build or repair a missing/partial Bowtie1 genome_index from the bundled Bowtie2 genome_index2, and provide the lincrnas.txt alias when only lncrnas.txt is present plus legacy optional FusionCatcher filter aliases/empty files, expose STAR-Fusion through a Perl runtime compatible with `Set::IntervalTree` and tolerate already-patched STAR-Fusion version-check scripts, include Perl `vendor_perl` paths for modules such as `common::sense`, resolve the deployed STAR-Fusion executable when it is not on PATH, repair required STAR-Fusion sidecar components `FusionFilter` and `FusionAnnotator`, expose the STAR executable used by STAR-Fusion, and configure Arriba to use fixed blacklist assets plus Arriba 2.5.1 argument names before production runs. Licensed tools are configured only when the user has supplied a legal local installation.
8. Run Doctor and evaluate required tools, alternative capability groups, critical references and sidecars against the requested tier. Missing required references can never produce READY. If the selected Conda installation is readable but its named splice environments are not writable, create user-writable portable environments under the selected tools root and record `CONDA_ENVS_PATH` in generated wrappers; never require sudo or assume a username.
9. Write deployment checkpoints, configuration fixes, and an installation status delta for safe resume/review. Full installs default to no wall-clock timeout; `--install-timeout SECONDS` can be supplied when an operator wants a bounded run. If the macro is interrupted or times out, it terminates the whole installer process group before recording an `INTERRUPTED` or `TIMEOUT` checkpoint so `resume` can continue from a controlled state.
10. Write `deployment_report.md`, machine-readable status files, and an audit log.

The local `env_tool/conda_pkgs` directory is a persistent package cache, not an
assumed complete offline distribution. Skill1 first seeds it from
`--conda-pkgs-source` when supplied and automatically probes common
`conda_pkgs` locations below a mounted shared asset root. Packages already in
that cache are reused; only missing packages are fetched from the environment
file's configured channels, including Bioconda where required. Conda is
configured with extended connect/read timeouts, exponential backoff and 12 HTTP
retries. If a step still reports `IncompleteRead`, `Connection broken`, timeout,
HTTP 429 or transient 5xx errors, the idempotent installer step is retried while
retaining downloaded package-cache content. A non-network solver/package error
is not retried blindly.

For `conda env create/update -f ENVIRONMENT_YML`, the YAML `channels` section is
authoritative. Skill1 must not append `--override-channels` or `-c` flags to
these environment subcommands, because current Conda releases can reject that
combination and it can also change the environment's declared channel order.
Bioconductor data packages require a registered Conda package-cache entry, not
a manually extracted look-alike directory. Before creating an environment that
contains such a package, Skill1 prefetches the YAML transaction with
`conda create --download-only`, verifies `info/repodata_record.json`, downloads
and checksum-validates the data archive into the shared install cache, patches
the registered post-link helper to use that archive, and performs the final
link transaction with Conda offline. Readiness must load the corresponding R
data package (for ASCAT, both `ASCAT` and `GenomeInfoDbData`) so a failed or
network-dependent post-link cannot be reported as installed.

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

For public fixed-asset synchronization, Skill1 defaults to the public Hugging
Face Dataset `open-neo/open-neo-public-assets` when no explicit site asset source
is supplied. Approved `install`, `repair`, and `resume` runs with
`--allow-download` use the official `hf download` CLI and its local cache for
robust multi-hour downloads and cross-process resume. They never use `curl` for
the large Dataset archive. Skill1 installs `huggingface-hub` as a core runtime
dependency; if `hf` is unavailable, stop with an actionable prerequisite error
instead of restarting a partial large download through another client. Runs verify
every published SHA-256, stream-extract it under the selected reference root,
and synchronize the separately published RSEM and SpliceMutr BSgenome additions.
A repository-commit marker makes resume idempotent; complete files and an already
verified extraction are skipped. Select a Hugging Face endpoint with
`--hf-endpoint` or the `HF_ENDPOINT` environment variable (for example,
`export HF_ENDPOINT=https://hf-mirror.com`); the CLI option takes precedence
and the official `https://huggingface.co` endpoint remains the fallback.
Override the Dataset repository with `--public-asset-repo`, revision with `--public-asset-revision`, locations with
`--public-asset-root`/`--public-asset-cache`, or disable this fallback with
`--no-sync-public-assets`.

Do not delete the Hugging Face local cache after an interrupted run. Re-run the
same Skill mode with the same `--public-asset-cache`; `hf download` validates
cached chunks and continues rather than restarting from zero. The deployment
marker records `download_method=hf download` for auditability.

An explicit site source still takes precedence. Set `--asset-source-root
/mounted/assets` for a local/mounted tree, or set both `--asset-source-host
user@host` and `--asset-source-root /path/on/source/host` for remote rsync. The
repository manifest uses the portable placeholder `/srv/neoag-assets/source`;
Skill1 rewrites it only after an explicit source root is supplied. No workflow
may default to a lab-specific host or username. Use `--asset-ssh-key
~/.ssh/id_ed25519` or `OPEN_NEO_ASSET_SSH_KEY` when the selected site source
requires a key.

The public Dataset intentionally excludes `hla`, `lohhla`, `predictors`,
`presentation`, `tools`, licensed containers, NetMHCpan, NetMHCstabpan,
NetChop, POLYSOLVER and novoalign material. Skill1 never substitutes public
downloads for these assets; the operator must stage authorized copies under
`--licensed-root`, and readiness remains blocked when a hard licensed
requirement is absent.

Licensed presentation assets should be staged under the selected licensed root
or the selected asset source before installation. For NetChop this means both
`netchop-3.1d.Linux.tar.gz` and the institutional license file are expected in
a readable `netchop/` directory, or `--netchop-archive` must point to the
authorized archive. SNAF source snapshot downloads are network fallback only;
Skill1 verifies the downloaded tarball and retries before failing. SpliceMutr hg38 BSgenome is treated as a synchronizable open asset; if it is absent from the selected asset source, approved runs with `--allow-download` download and stage the Bioconductor package before reference verification.
For production SpliceMutr, readiness also checks the GTEx per-sample junction-count matrix used for public-proxy cohorts, the pinned LeafCutter/LeafViz scripts, a build-matched exon/LeafViz annotation set, and a TxDb SQLite generated from the same GENCODE GTF used by RNA alignment. The SpliceMutr R runtime and LeafCutter R runtime may be separate and must be recorded explicitly. Validate `BSgenome.Hsapiens.UCSC.hg38` as the package and `Hsapiens` as its exported genome object; never assume those two names are interchangeable. A compact aggregate normal-junction catalog alone is sufficient for safety lookup but not for constructing independent SpliceMutr control samples.

For a NAS-managed Conda installation, use:

```bash
open-neo install-check \
  --project-root . \
  --mode verify \
  --conda-base /nas/path/to/miniforge3 \
  --outdir work/install-check
```

For a new machine using the public fixed-asset mirror:

```bash
open-neo install-check \
  --project-root . \
  --mode install \
  --deployment-tier full \
  --approved \
  --allow-download \
  --public-asset-root /srv/open-neo/refs \
  --public-asset-cache /srv/open-neo/download-cache \
  --licensed-root /srv/open-neo/licensed_tools \
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
- Treat Manta, SvABA and GRIDSS as the required-all `dna_sv` capability group for `full`/`all-open` readiness. Install or register each caller independently, retain its native output directory, and probe the real executable or configured container entrypoint. Use `neoag-sv` for Manta/SvABA and an isolated `neoag-gridss` environment for GRIDSS's Java/R/Bioconductor dependency chain; expose stable wrappers so this physical split is transparent to Skill2. The group is READY only when all three members pass; one available caller must not satisfy the group. Record available and missing members in `tier_requirements.tsv` so resume repairs only the missing member rather than repeating completed installations.
- Default to `NEOAG_FUSION_INSTALL_MODE=easyfuse`: install one pinned EasyFuse source tree plus its `neoag-easyfuse` driver environment, and use EasyFuse's native module definitions for STAR, Arriba, STAR-Fusion and FusionCatcher. Readiness must check each internal caller definition and later preserve each caller's raw output, but must not clone or maintain standalone copies. `NEOAG_FUSION_INSTALL_MODE=standalone` (or `--standalone-fusion`) is an explicit compatibility fallback only; when selected, create/update `neoag-fusion` through the resumable Bioconductor data cache helper, defaulting to `genomeinfodbdata-1.2.11` and allowing `NEOAG_FUSION_BIOC_PACKAGE_KEY` to override the cache key.
- NetMHCstabpan readiness requires more than an executable path: verify the
  licensed data directory and runtime dependencies, then run one supported
  8-11mer HLA-I peptide through the real binary. Record executable, data root,
  version, allele, peptide length, exit status and parsed output in the smoke
  result. A help-only check is PARTIAL, not READY.
- Never use BAM-matcher's bundled GRCh37 loci with GRCh38 data. Build the portable panel with `scripts/build_bam_matcher_grch38_loci.sh`; the script requires exact dbSNP-ID mapping, biallelic SNVs, GRCh38 FASTA REF agreement, and emits a metadata manifest plus checksums.
- PRIME, MixMHCpred and BigMHC use pinned complete source assets, not model-only or wrapper-only copies. Required markers are `lib/run_PRIME.pl`, `MixMHCpred`, and `src/predict.py`; BigMHC must also contain `models/`. Resume reuses these synchronized directories instead of downloading them again.
- Full installs must leave FACETS, Sequenza, PURPLE/HMFTOOLS, SpecHLA, LOHHLA, and fusion callers callable without relying on interactive shell state. LOHHLA readiness also requires a readable licensed Polysolver tree, project-local BAM index/link handling, absolute LOHHLA inputs, fresh explicit LOHHLA work directories for reruns, and preservation/regeneration of root mpileup intermediates instead of reusing partial failed directories. Its selected R runtime must load `optparse`, `Rsamtools`, `Biostrings`, and `seqinr`, while the same non-interactive launcher PATH must resolve `bedtools` and `samtools`; testing these pieces in separate environments is insufficient. The LOHHLA smoke check must also parse a synthetic `CopyNumLoc` with R and verify that the leading tumor ID becomes a row name under the upstream `read.table(..., header=TRUE)` behavior. `conf/tools.env.local.sh` receives a marked install-check block for these paths, and `bin/bedtools` is generated when the portable bedtools binary requires the bundled Conda C++ runtime. EasyFuse readiness requires its pinned source, driver Nextflow entrypoint, reference marker, and internal Arriba/STAR-Fusion/FusionCatcher module environment definitions. Legacy EasyFuse layouts additionally require their referenced compatibility YAMLs and compatibility repairs; current module-native layouts must not be forced through legacy standalone installs. Runtime verification still checks each internal caller and records its raw output independently. NetMHCstabpan readiness must install/register a callable licensed/shim executable from the deployed licensed tools root or shared asset root, preserve the supplied license file when present, and provide a `gawk`-compatible runtime entry. NetMHCpan readiness must validate a real peptide-HLA prediction and prefer the container wrapper when native NetMHCpan 4.2c is incompatible with the host glibc. NetChop readiness must install from an authorized `netchop-3.1d.Linux.tar.gz` discovered in the licensed root or shared asset root, keep the license file alongside it when supplied, and pin `NETCHOP_HOME`/`NEOAG_NETCHOP_BIN` into `conf/tools.env.local.sh`. Re-running `install` or `resume` replaces only the marked block/wrapper and skips already synchronized assets when their checkpoints and markers are still valid. For prediction/full tiers, NetMHCpan, MHCflurry, NetMHCstabpan, and NetChop are all hard readiness requirements because raw heavy production manifests require all four presentation predictors.

## Outputs

- `environment_inventory.tsv`
- `doctor/doctor_status.json`
- `deployment_status.tsv`
- `claude_code_status.tsv` and `claude_code_status.json`, including requested
  state, readiness, version, binary path, and the nested install report path
- `tier_requirements.tsv`, `deployment_delta.tsv`
- `production_run_readiness.tsv`, covering full-run smoke compatibility for
  Sequenza, PURPLE/HMFTOOLS, FACETS, STAR/EasyFuse, the Manta/SvABA/GRIDSS
  DNA-SV group, BAM-matcher, HLA-LA, and
  optional NeoAg Gateway health
- `deployment_checkpoint.json` for mutating modes
- `deployment_report.md`
- `run_issue_log.json`, with one entry per observed installation/readiness problem: stable ID, affected Skill, symptom, root cause, attempted workaround, durable fix, validation evidence, and status. Append across resume runs; do not overwrite unresolved history.
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
- A file being present or executable is not sufficient readiness evidence. For binaries and wrappers, record a version/help or minimal functional smoke result in `run_issue_log.json`; malformed files and wrappers whose dependencies are missing remain unresolved.
- The launcher must resolve Python 3.11+ with `tomllib` from `NEOAG_PYTHON`,
  configured Conda/deployment roots, the project virtual environment, or PATH;
  it must fail clearly instead of starting the macro with an older interpreter.
