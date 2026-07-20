# run_pipeline — NeoAg Nextflow Pipeline Selection & Execution

## Trigger

Invoke this skill when the user wants to **run a neoantigen prediction pipeline** and
needs help choosing the correct Nextflow workflow based on their input data. Typical
user requests:

- "Run the pipeline on my tumor/normal BAM files"
- "Run neoantigen prediction from my VCF"
- "Run the full pipeline with QC"
- "I have BAM files and known HLA alleles, which pipeline should I use?"
- Any request involving `neoag-nextflow run` or `nextflow run`

## Workflow

### Step 1 — Collect Input Metadata

Ask the user (or inspect provided paths) to determine these key facts:

| Question | Why it matters |
|----------|---------------|
| What input files do you have? (BAM ×2, VCF, TOML config, FASTQ) | Determines starting point |
| Do you need HLA typing, or do you already know the alleles? | OptiType vs manual HLA |
| Do you need QC analysis? (LOH + tumor purity) | main_all vs main_all_qc |
| What's the sample ID? | Required by all pipelines |
| Where is the reference FASTA? | Required for variant calling |
| Do you have a dbSNP VCF? (optional, for FACETS in QC mode) | Only needed for QC |
| Do you want Docker or Conda runtime? | Profile selection |
| Need to resume a previous run? | `-resume` flag |

If the user provides a file path, verify it exists before proceeding.

### Step 2 — Match to Workflow

Use this decision table. Find the **first row that matches**:

```
 INPUT TYPE          │ HLA SOURCE      │ QC?  │ WORKFLOW
─────────────────────┼─────────────────┼──────┼──────────────────────────
 TOML config only    │ pre-configured  │  no  │ main_fromVCF.nf
 (VCF+HLA in TOML)   │ in TOML         │      │
─────────────────────┼─────────────────┼──────┼──────────────────────────
 TOML config         │ OptiType        │  no  │ main_fromVCF_nohla.nf
 + BAM or FASTQ      │ auto-typing     │      │
─────────────────────┼─────────────────┼──────┼──────────────────────────
 tumor BAM           │ OptiType        │  no  │ main_all.nf
 + normal BAM        │ auto-typing     │      │
─────────────────────┼─────────────────┼──────┼──────────────────────────
 tumor BAM           │ OptiType        │ yes  │ main_all_qc.nf
 + normal BAM        │ auto-typing     │      │
─────────────────────┼─────────────────┼──────┼──────────────────────────
 tumor BAM           │ user-provided   │  no  │ main_all_nohla.nf
 + normal BAM        │ (manual)        │      │
 + HLA alleles       │                 │      │
─────────────────────┴─────────────────┴──────┴──────────────────────────
```

### Step 3 — Map Parameters

Based on the chosen workflow, collect these parameters:

#### Workflow 1: `main_fromVCF.nf`
```
--run_config <path>       TOML config with [sample] and [inputs] sections
--sample_id <id>           Override TOML sample.id (optional)
--outdir <dir>             Output directory
```

#### Workflow 2: `main_fromVCF_nohla.nf`
```
--run_config <path>       TOML config (HLA will be auto-filled)
--input_bam <path>        BAM for HLA typing
--sample_id <id>           Sample identifier
--outdir <dir>             Output directory
 --skip_hla_typing          (optional) skip OptiType if HLA is in TOML
```

#### Workflow 3: `main_all.nf`
```
--normal_bam <path>       Normal/blood BAM
--tumor_bam <path>        Tumor BAM
--sample_id <id>           Sample identifier
--reference_fasta <path>  Reference genome FASTA
--outdir <dir>             Output directory
--tumor_sample_name <str> BAM read-group tumor sample name (default: TUMOR)
--normal_sample_name <str>BAM read-group normal sample name (default: NORMAL)
```

#### Workflow 4: `main_all_qc.nf`
```
Same as main_all.nf, plus:
--dbsnp_vcf <path>        dbSNP/common SNP VCF for FACETS (optional)
--skip_qc                  Skip QC subworkflows but keep main pipeline
```

#### Workflow 5: `main_all_nohla.nf`
```
--normal_bam <path>       Normal/blood BAM
--tumor_bam <path>        Tumor BAM
--hla_alleles <str>       Comma-separated, e.g. "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02"
--sample_id <id>           Sample identifier
--reference_fasta <path>  Reference genome FASTA
--outdir <dir>             Output directory
```

### Step 4 — Construct the Command

All workflows share these common flags:

```
-c conf/main_full.config         ← always include
-profile docker                  ← if NEOAG_RUNNER_MODE=docker
-profile conda                   ← if NEOAG_RUNNER_MODE=conda
-resume                          ← if resuming
-with-report <path>              ← optional execution report
```

**Command template:**

```bash
cd /mnt/disk_c/data_transfer/users/samba_wb/indev/neo

[ NEOAG_RUNNER_MODE=docker ] \
bin/neoag-nextflow run workflows/<chosen_workflow>.nf \
  --sample_id <SAMPLE_ID> \
  <...workflow-specific params...> \
  -c conf/main_full.config \
  [ -profile docker ] \
  [ -resume ]
```

Or use the unified launcher:

```bash
bash scripts/run_pipeline.sh \
  --workflow <chosen_workflow> \
  --sample_id <SAMPLE_ID> \
  <...workflow-specific params...>
```

### Step 5 — Validate Before Presenting

Before showing the command to the user, check:

1. **All required files exist** — verify each input BAM/VCF/TOML/FASTA path with `test -f`
2. **Reference genome has .fai and .dict** — `${ref}.fai` and `${ref%.fa*}.dict` must exist
3. **TOML config is valid** (if applicable) — check that `[sample]` has `id` and `[inputs]` has required keys
4. **bin/neoag-nextflow is executable** — verify the launcher exists
5. **conf/main_full.config exists** — verify the config file path

Report any missing files or issues. Do NOT proceed with missing requirements.

### Step 6 — Present to User

Show the user:

1. **The chosen workflow** and **why** it was selected (one sentence)
2. **The complete nextflow command** (copy-paste ready)
3. **Expected output locations**:
   - Scoring results: `<outdir>/`
   - HLA typing (if auto): `<outdir>/hla_typing/`
   - QC (if enabled): `<outdir>/qc/`
4. **Any warnings** (e.g., "dbsnp_vcf not provided — FACETS will produce stub output")

Wait for user confirmation. Do NOT execute unless the user explicitly asks.

---

## Reference: Quick Workflow Comparison

| Workflow | Input | HLA | Variant Calling | QC | Runtime (approx) |
|----------|-------|-----|-----------------|----|------------------|
| `main_fromVCF.nf` | TOML | in TOML | skipped | no | 2-4 h |
| `main_fromVCF_nohla.nf` | TOML + BAM | OptiType | skipped | no | 2-5 h |
| `main_all.nf` | BAM×2 | OptiType | Mutect2 | no | 8-24 h |
| `main_all_qc.nf` | BAM×2 | OptiType | Mutect2 | yes | 10-30 h |
| `main_all_nohla.nf` | BAM×2 + HLA | manual | Mutect2 | no | 6-20 h |

## Reference: Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `NEOAG_RUNNER_MODE` | `conda` or `docker` | `conda` |
| `NEOAG_REFERENCE_FASTA` | Reference genome path | (required) |
| `NEOAG_DBSNP_VCF` | dbSNP VCF for FACETS | (optional) |
| `NXF_HOME` | Nextflow metadata dir | `work/.nextflow_home` |

## Reference: TOML Config Format

A `main_fromVCF.nf` run config must contain:

```toml
[sample]
id = "SAMPLE001"
profile = "default"

[inputs]
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
...
```

See `conf/run.mycase.toml` and `conf/run.private.example.toml` for complete examples.
