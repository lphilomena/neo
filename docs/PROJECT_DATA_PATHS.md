# Project Data Paths

This table lists non-fixture data paths referenced by the project. Paths are relative to the repository root unless explicitly marked as external or private.

Current staged bundle root:

```bash
source /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/neodata.env.sh
```

Staging manifest:

`/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/manifests/data_staging_manifest.md`

## VEP cache

`NEOAG_VEP_CACHE` is the **cache root directory**, not the release subdirectory. VEP is called with `--dir_cache "$NEOAG_VEP_CACHE"` and `--cache_version 105`.

Expected layout:

```text
$NEOAG_VEP_CACHE/
└── homo_sapiens/
    └── 105_GRCh38/
        ├── info.txt
        └── 1/ 2/ ... 22/
```

| Field | Value |
| --- | --- |
| Species | `homo_sapiens` |
| Assembly | `GRCh38` |
| Release | `105` (`NEOAG_VEP_CACHE_VERSION`) |
| Package | `homo_sapiens_vep_105_GRCh38.tar.gz` |
| Typical size | about 12–16 GB |
| Source URL | `https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz` |

On this server:

| Role | Path |
| --- | --- |
| Canonical physical cache root | `/home/na/project/neoantigen/neoag_event_pipeline_v03_rc_artifact_quarantine_20260622_091158/data/vep` |
| Release directory | `/home/na/project/neoantigen/neoag_event_pipeline_v03_rc_artifact_quarantine_20260622_091158/data/vep/homo_sapiens/105_GRCh38` |
| Staged bundle cache root | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/vep` |
| Staged symlink | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/vep/homo_sapiens` → quarantine `data/vep/homo_sapiens` |
| Default in `conf/tools.env.sh` | `$NEOAG_TOOLS_ROOT/data/vep`, with automatic fallback to the quarantine cache root when `homo_sapiens/` is missing there |
| Fresh install target from `scripts/install_vep_cache.sh` | `~/.vep` (set `NEOAG_VEP_CACHE=$HOME/.vep` after install) |

Verify:

```bash
test -f "$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/info.txt"
```

| Category | Path used in project | Data / reference | Used for | Bundled in lightweight repo |
| --- | --- | --- | --- | --- |
| Current staged bundle | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata` | NeoAg artifact/reference bundle root | Source `neodata.env.sh` before real-data runs | No |
| Current staged env | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/neodata.env.sh` | Environment variables for staged bundle | Sets `NEOAG_TOOLS_ROOT`, `NEOAG_REFERENCE_FASTA`, `NEOAG_VEP_CACHE`, `CTAT_GENOME_LIB`, `NEOAG_EASYFUSE_REF` | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/Homo_sapiens_assembly38.fasta` | GRCh38 reference FASTA, linked from `/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | VEP, GATK, SV peptide building | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/Homo_sapiens_assembly38.fasta.fai` | GRCh38 FASTA index | VEP/GATK/SV workflows that need indexed FASTA | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/Homo_sapiens_assembly38.dict` | GRCh38 sequence dictionary | GATK Mutect2 and related WES workflows | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/GRCh38.fa` | GRCh38 chr-prefixed FASTA, linked from `/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa` | SV/fusion real-data scripts | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/GRCh38.fa.fai` | GRCh38 chr-prefixed FASTA index | SV/fusion real-data scripts | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/Homo_sapiens.GRCh38.pep.all.fa` | Ensembl protein FASTA | Peptide safety normal/reference proteome screen | No |
| Current staged reference | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/Homo_sapiens.GRCh38.110.gtf` | GENCODE/Ensembl GTF annotation, linked from EasyFuse reference | SV/fusion peptide generation and gene CNV/LOH helpers | No |
| Current staged VEP cache root | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/vep` | VEP cache root for staged bundle; `homo_sapiens` symlinked to quarantine cache | Sets `NEOAG_VEP_CACHE` in `neodata.env.sh` | No |
| Current staged VEP cache release | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/vep/homo_sapiens/105_GRCh38` | Ensembl indexed VEP cache, GRCh38 release 105 | Offline VEP annotation (`--dir_cache` + `--cache_version 105`) | No |
| Canonical VEP cache root | `/home/na/project/neoantigen/neoag_event_pipeline_v03_rc_artifact_quarantine_20260622_091158/data/vep` | Physical cache root used by local fallback and staged symlink | Offline VEP annotation | No |
| Canonical VEP cache release | `/home/na/project/neoantigen/neoag_event_pipeline_v03_rc_artifact_quarantine_20260622_091158/data/vep/homo_sapiens/105_GRCh38` | Ensembl indexed VEP cache files (`info.txt`, per-chromosome dirs) | Offline VEP annotation | No |
| Current staged CTAT | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/ctat/current` | CTAT genome library, linked from `/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play` | STAR-Fusion workflows | No |
| Current staged EasyFuse | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/easyfuse_ref_v4` | EasyFuse v4 reference, linked from `/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/easyfuse_ref_v4` | EasyFuse workflows | No |
| Current staged private data | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/external/chenxiaoliang_data` | Chenxiaoliang private data root, linked from `/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data` | Local patient scripts only | No |
| Current missing data | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/dbsnp_chr.vcf.gz` | dbSNP/common SNP VCF | Not staged yet; required for FACETS `snp-pileup` if using that path | No |
| Current missing data | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/af-only-gnomad.hg38.vcf.gz` | gnomAD AF VCF | Not staged yet; required for GATK Mutect2 filtering | No |
| Current missing data | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/pon.vcf.gz` | Panel-of-normals VCF | Not staged yet; required for GATK Mutect2 filtering when PoN is enabled | No |
| Current missing data | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/capture.bed` | Production capture BED | Not staged yet; WES/panel-specific | No |
| Current missing data | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata/data/ref/hg38/hla.txt` | Sample HLA allele file | Not staged because it is sample-specific | No |
| Short-path alignment reference | `/home/na/ref/hg38/GRCh38.fa` | Short path GRCh38 FASTA for `bwa-mem2` | `scripts/run_sv_wgs_chenxiaoliang.sh` default `REF_SHORT_DIR`; avoids bwa-mem2 path length issue | No |
| Short-path alignment reference | `/home/na/ref/hg38/GRCh38.fa.bwt.2bit.64` | `bwa-mem2` index for short path GRCh38 FASTA | WGS SV alignment with `ALIGNER=bwa-mem2` | No |
| Large reference | `data/ref/hg38/Homo_sapiens_assembly38.fasta` | GRCh38 reference FASTA | VEP, GATK, SV peptide building | No |
| Large reference | `data/ref/hg38/GRCh38.fa` | GRCh38 FASTA alternative path used by scripts | SV/fusion real-data scripts | No |
| Large reference | `data/ref/hg38/Homo_sapiens_assembly38.fasta.fai` | GRCh38 FASTA index | VEP/GATK/SV workflows that need indexed FASTA | No |
| Large reference | `data/ref/hg38/Homo_sapiens_assembly38.dict` | GRCh38 sequence dictionary | GATK Mutect2 and related WES workflows | No |
| Large reference | `data/ref/hg38/dbsnp_chr.vcf.gz` | dbSNP/common SNP VCF | FACETS `snp-pileup` and CNV/LOH workflows | No |
| Large reference | `data/ref/hg38/af-only-gnomad.hg38.vcf.gz` | gnomAD allele-frequency VCF | GATK Mutect2 filtering | No |
| Large reference | `data/ref/hg38/pon.vcf.gz` | Panel-of-normals VCF | GATK Mutect2 filtering | No |
| Large reference | `data/ref/hg38/Homo_sapiens.GRCh38.pep.all.fa` | Ensembl/GENCODE protein FASTA | Peptide safety normal/reference proteome screen | No |
| Large reference | `data/ref/hg38/Homo_sapiens.GRCh38.110.gtf` | GENCODE/Ensembl GTF annotation | SV/fusion peptide generation | No |
| Large reference | `data/ref/hg38/capture.bed` | WES/panel capture BED | WES SV Phase 1.5 and interval-restricted workflows | No |
| Large reference | `data/ref/hg38/hla.txt` | HLA allele file, one allele per line | Peptide prediction and SV workflows | No |
| Large reference | `data/vep/` | VEP cache root (`homo_sapiens/105_GRCh38/` underneath) | Offline VEP annotation | No |
| Large reference | `data/external/` | Site-managed external inputs | Private/local deployments | No |
| Site evidence | `/path/to/normal_expression.tsv` | Site-generated normal tissue expression table | Peptide safety evidence | No |
| Site evidence | `/path/to/normal_hla_ligands.tsv` | Site-generated normal HLA ligand table | Peptide safety evidence | No |
| Example/private data | `data/examples/HCC1395/HCC1395_inputs/annotated.expression.vcf.gz` | HCC1395 annotated VCF example path | HCC1395 example configs/tests when staged | No |
| Example/private data | `data/examples/HCC1395/HCC1395_inputs/HCC1395.splice_junctions.tsv` | HCC1395 splice junction example path | HCC1395 splice example configs/tests when staged | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.fasta` | Site reference FASTA configured by env | Production VEP/GATK/SV workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.fasta.fai` | Site reference FASTA index | Production VEP/GATK/SV workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.dict` | Site reference sequence dictionary | Production GATK workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/dbsnp_chr.vcf.gz` | Site dbSNP/common SNP VCF | FACETS `snp-pileup` and CNV/LOH workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/af-only-gnomad.hg38.vcf.gz` | Site gnomAD allele-frequency VCF | GATK Mutect2 filtering | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/pon.vcf.gz` | Site panel-of-normals VCF | GATK Mutect2 filtering | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/Homo_sapiens.GRCh38.pep.all.fa` | Site Ensembl/GENCODE protein FASTA | Peptide safety normal/reference proteome screen | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/Homo_sapiens.GRCh38.110.gtf` | Site GENCODE/Ensembl GTF annotation | SV/fusion peptide generation | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/capture.bed` | Site WES/panel capture BED | WES SV Phase 1.5 and interval-restricted workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/hg38/hla.txt` | Site HLA allele file, one allele per line | Peptide prediction and SV workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/vep` | Site VEP cache root configured by env; release dir is `${NEOAG_TOOLS_ROOT}/data/vep/homo_sapiens/105_GRCh38` | Production VEP workflows | No |
| External artifact root | `${NEOAG_TOOLS_ROOT}/data/ref/ctat` | CTAT genome library path | STAR-Fusion/fusion workflows | No |
| External artifact root | `${CTAT_GENOME_LIB}` | Active CTAT genome library directory | STAR-Fusion workflows | No |
| External path | `${NEOAG_DBSNP_VCF}` | Configured dbSNP/common SNP VCF | FACETS `snp-pileup` and CNV/LOH workflows | No |
| External path | `${NEOAG_NORMAL_PROTEOME_FASTA}` | Configured Ensembl/GENCODE protein FASTA | Peptide safety normal/reference proteome screen | No |
| External reference | `${NEOAG_EASYFUSE_REF}` | EasyFuse reference directory | EasyFuse workflows | No |
| EasyFuse reference | `${NEOAG_SHARED_REF_DIR}/easyfuse_ref_v4` | Shared EasyFuse v4 reference directory | Preferred source for `NEOAG_EASYFUSE_REF` in `conf/tools.env.sh` | No |
| EasyFuse reference archive | `${NEOAG_SHARED_REF_DIR}/easyfuse_ref_v4.tar.gz` | EasyFuse v4 reference archive | Staging/unpacking EasyFuse references | No |
| EasyFuse local default | `/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/easyfuse_ref_v4` | Local EasyFuse v4 reference directory | Chenxiaoliang EasyFuse scripts when `NEOAG_EASYFUSE_REF` is unset | No |
| EasyFuse reference file | `${NEOAG_EASYFUSE_REF}/BEFORE_EXECUTING_EASYFUSE` | EasyFuse reference readiness marker | Sanity check before EasyFuse runs | No |
| EasyFuse reference file | `${NEOAG_EASYFUSE_REF}/Homo_sapiens.GRCh38.110.gff3.db` | EasyFuse annotation database | `--annotation_db` for EasyFuse Nextflow run | No |
| EasyFuse reference file | `${NEOAG_EASYFUSE_REF}/Homo_sapiens.GRCh38.110.gtf.tsl` | Transcript support level annotation | `--reference_tsl` for EasyFuse Nextflow run | No |
| EasyFuse reference file | `${NEOAG_EASYFUSE_REF}/Homo_sapiens.GRCh38.110.gtf` | GENCODE/Ensembl gene annotation | SV/fusion peptide generation and gene CNV/LOH helpers | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/liver_0520_WGS_shortReads/somatic/*.vcf.gz` | Chenxiaoliang somatic VCFs | Local patient scripts only | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/liver_0520_WGS_shortReads/seq_liver_26052/*.bam` | Chenxiaoliang tumor WGS BAMs | Local FACETS/SpecHLA/SV scripts only | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/blood_0427/retransfer_6927_7362/*/*.bam` | Chenxiaoliang normal WGS BAMs | Local FACETS/SpecHLA/SV scripts only | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/liver_0520_WGS_shortReads/seq_liver_26052/*_1.fq.gz` | Chenxiaoliang tumor FASTQ R1 | Local SV/fusion scripts only | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/liver_0520_WGS_shortReads/seq_liver_26052/*_2.fq.gz` | Chenxiaoliang tumor FASTQ R2 | Local SV/fusion scripts only | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/blood_0427/retransfer_6927_7362/*_1.fq.gz` | Chenxiaoliang normal FASTQ R1 | Local SV scripts only | No |
| Private Chenxiaoliang data | `${CHENXIAOLIANG_DATA_ROOT}/data/blood_0427/retransfer_6927_7362/*_2.fq.gz` | Chenxiaoliang normal FASTQ R2 | Local SV scripts only | No |
| PolySolver reference | `/home/na/project/neoantigen/software/polysolver/data/abc_complete.fasta` | HLA allele FASTA for patient HLA FASTA generation | Local helper script default | No |
| PolySolver reference | `/home/na/project/neoantigen/software/polysolver/data/complete` | PolySolver complete data directory | Local helper script default | No |
| LOHHLA reference | `tools/lohhla/data/hla.dat` | LOHHLA HLA exon location file | `scripts/run_lohhla_example.sh` | No |

