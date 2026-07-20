from __future__ import annotations
import os
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable

ROOT = Path(__file__).resolve().parents[3]


class RunnerMode(Enum):
    CONDA = "conda"
    DOCKER = "docker"


def resolve_runner_mode() -> RunnerMode:
    """Resolve runner mode from NEOAG_RUNNER_MODE env var.

    Values: ``conda`` (default), ``docker``.
    In Docker mode, tools that lack a public image fall back to conda.
    """
    raw = os.environ.get("NEOAG_RUNNER_MODE", "conda").strip().lower()
    if raw in {"docker", "container"}:
        return RunnerMode.DOCKER
    return RunnerMode.CONDA


# Docker images for tools that have official or biocontainers images.
# Maps tool registry name → docker image:tag (explicit version, never :latest).
# Tools NOT listed here fall back to conda mode regardless of NEOAG_RUNNER_MODE.
DOCKER_IMAGES: dict[str, str] = {
    # ---- 官方 Docker 镜像 ----
    "vep": "ensemblorg/ensembl-vep:release_115.2",
    "pvacseq": "griffithlab/pvactools:6.1.1",
    "pvacfuse": "griffithlab/pvactools:6.1.1",
    "pvacsplice": "griffithlab/pvactools:6.1.1",
    "optitype": "fred2/optitype:release-v1.3.1",
    "gatk": "broadinstitute/gatk:4.6.2.0",
    "arriba": "uhrigs/arriba:2.5.1",
    "easyfuse": "tronbioinformatics/easyfuse:1.3.7",
    "gridss": "gridss/gridss:2.13.2",
    "delly": "dellytools/delly:v2.3.0",
    "sequenza": "sequenza/sequenza:3.99rc1",
    # ---- Biocontainers 镜像 ----
    "lohhla": "quay.io/biocontainers/lohhla:20171108--hdfd78af_3",
    "ascat": "quay.io/biocontainers/ascat:2.5.2--r40hdfd78af_3",
    "star_fusion": "quay.io/biocontainers/star-fusion:1.15.1--hdfd78af_1",
    "fusioncatcher": "quay.io/biocontainers/fusioncatcher:1.33b--hdfd78af_0",
    "manta": "quay.io/biocontainers/manta:1.6.0--h9ee0642_3",
    "svaba": "quay.io/biocontainers/svaba:1.2.0--h69ac913_1",
    # ---- 基础工具 (Biocontainers) ----
    "samtools": "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0",
    "tabix": "quay.io/biocontainers/tabix:1.11--hdfd78af_0",
    # ---- 自建镜像 (deploy/build_containers/) ----
    "netmhcpan": "neoag-netmhcpan:4.2c-ubuntu22.04",
    "netmhcstabpan": "neoag-netmhcstabpan:1.0-ubuntu22.04",
    "mhcflurry": "neoag-base-bioinfo:ubuntu22.04",
    "prime": "neoag-base-bioinfo:ubuntu22.04",
    "bigmhc_im": "neoag-base-bioinfo:ubuntu22.04",
    "deepimmuno": "neoag-base-bioinfo:ubuntu22.04",
    "pyclone": "neoag-base-bioinfo:ubuntu22.04",
    "facets": "neoag-purple-suite:ubuntu22.04",
    "spechla": "neoag-spechla:ubuntu22.04",
}


def tool_docker_image(tool_name: str) -> str | None:
    """Return the Docker image:tag for *tool_name*, or None if unavailable."""
    return DOCKER_IMAGES.get(tool_name)


@dataclass(frozen=True)
class ToolSpec:
    name: str
    category: str
    executable: str
    description: str
    version_args: tuple[str, ...] = ("--version",)
    fixture_outputs: dict[str, str] = field(default_factory=dict)
    docs_url: str = ""


def _fx(name: str) -> str:
    return str(ROOT / "data" / "fixtures" / name)


TOOL_REGISTRY: dict[str, ToolSpec] = {
    "pvacseq": ToolSpec(
        name="pvacseq",
        category="neoantigen",
        executable="pvacseq",
        description="pVACseq: SNV/InDel neoantigen prediction from VCF + HLA",
        version_args=(),  # --version is very slow to import
        fixture_outputs={"aggregated": _fx("pvacseq_aggregated.tsv")},
        docs_url="https://pvactools.readthedocs.io/en/latest/pvacseq.html",
    ),
    "pvacfuse": ToolSpec(
        name="pvacfuse",
        category="neoantigen",
        executable="pvacfuse",
        description="pVACfuse: fusion neoantigen prediction",
        version_args=(),
        fixture_outputs={"aggregated": _fx("pvacfuse_aggregated.tsv")},
        docs_url="https://pvactools.readthedocs.io/en/latest/pvacfuse.html",
    ),
    "pvacsplice": ToolSpec(
        name="pvacsplice",
        category="neoantigen",
        executable="pvacsplice",
        description="pVACsplice: splice junction neoantigen prediction from VCF + RegTools",
        version_args=(),
        fixture_outputs={"aggregated": _fx("pvacsplice_aggregated.tsv")},
        docs_url="https://pvactools.readthedocs.io/en/latest/pvacsplice.html",
    ),
    "netmhcpan": ToolSpec(
        name="netmhcpan",
        category="presentation",
        executable="netMHCpan",
        description="NetMHCpan binding and EL prediction",
        version_args=(),
        fixture_outputs={"xls": _fx("netmhcpan_example.xls")},
        docs_url="https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1",
    ),
    "mhcflurry": ToolSpec(
        name="mhcflurry",
        category="presentation",
        executable="mhcflurry-predict",
        description="MHCflurry presentation prediction",
        version_args=(),
        fixture_outputs={"csv": _fx("mhcflurry_predictions.csv")},
        docs_url="https://github.com/openvax/mhcflurry",
    ),
    "netmhcstabpan": ToolSpec(
        name="netmhcstabpan",
        category="presentation",
        executable="netMHCstabpan",
        description="NetMHCstabpan pMHC stability prediction",
        version_args=(),
        fixture_outputs={"tsv": _fx("netmhcstabpan_example.tsv")},
        docs_url="https://services.healthtech.dtu.dk/service.php?NetMHCstabpan-1.0",
    ),
    "prime": ToolSpec(
        name="prime",
        category="immunogenicity",
        executable="PRIME",
        description="PRIME class I immunogenicity predictor",
        version_args=(),
        fixture_outputs={"tsv": _fx("prime_example.tsv")},
        docs_url="https://github.com/GfellerLab/PRIME",
    ),
    "bigmhc_im": ToolSpec(
        name="bigmhc_im",
        category="immunogenicity",
        executable="bigmhc_predict",
        description="BigMHC_IM neoepitope immunogenicity predictor",
        version_args=(),
        fixture_outputs={"tsv": _fx("bigmhc_im_example.tsv")},
        docs_url="https://github.com/KarchinLab/bigmhc",
    ),
    "deepimmuno": ToolSpec(
        name="deepimmuno",
        category="immunogenicity",
        executable="deepimmuno-cnn.py",
        description="DeepImmuno-CNN class I immunogenicity predictor (9/10-mer pairs)",
        version_args=(),
        fixture_outputs={"tsv": _fx("deepimmuno_example.tsv")},
        docs_url="https://github.com/frankligy/DeepImmuno",
    ),
    "vep": ToolSpec(
        name="vep",
        category="annotation",
        executable="vep",
        description="Ensembl VEP variant consequence annotation",
        version_args=(),  # vep --version is slow; existence check only
        fixture_outputs={"appm": _fx("vep_appm.tsv")},
        docs_url="https://www.ensembl.org/vep",
    ),
    "lohhla": ToolSpec(
        name="lohhla",
        category="hla",
        executable="LOHHLA",
        description="HLA allele-specific LOH from WES/WGS",
        version_args=(),
        fixture_outputs={"hla_loh": _fx("hla_loh.tsv")},
        docs_url="https://bitbucket.org/mcferrine/lohhla/src/master/",
    ),
    "optitype": ToolSpec(
        name="optitype",
        category="hla",
        executable="optitype",
        description="OptiType: precision 4-digit HLA typing from NGS (DNA/RNA)",
        version_args=("--version",),
        fixture_outputs={"hla_csv": _fx("optitype_hla_typing.csv")},
        docs_url="https://github.com/FRED-2/OptiType",
    ),
    "facets": ToolSpec(
        name="facets",
        category="purity",
        executable="runFACETS.R",
        description="FACETS tumor purity and copy number",
        version_args=(),
        fixture_outputs={"purity": _fx("purity.tsv")},
        docs_url="https://github.com/mskcc/facets",
    ),
    "ascat": ToolSpec(
        name="ascat",
        category="purity",
        executable="ascat.R",
        description="ASCAT allele-specific copy number and purity",
        version_args=("--version",),
        fixture_outputs={"purity": _fx("purity.tsv")},
    ),
    "sequenza": ToolSpec(
        name="sequenza",
        category="purity",
        executable="sequenza-utils",
        description="Sequenza allele-specific copy number and purity/ploidy from WGS (bam2seqz + R sequenza)",
        version_args=(),
        fixture_outputs={"purity": _fx("purity.tsv")},
        docs_url="https://bitbucket.org/sequenza/sequenza/",
    ),
    "star_fusion": ToolSpec(
        name="star_fusion",
        category="fusion",
        executable="star-fusion-neoag",
        description="STAR-Fusion chimeric transcript detection",
        version_args=(),
        docs_url="https://github.com/STAR-Fusion/STAR-Fusion",
    ),
    "arriba": ToolSpec(
        name="arriba",
        category="fusion",
        executable="arriba",
        description="Arriba fusion detection from RNA-seq",
        version_args=(),
        docs_url="https://github.com/suhrig/arriba",
    ),
    "fusioncatcher": ToolSpec(
        name="fusioncatcher",
        category="fusion",
        executable="fusioncatcher-neoag",
        description="FusionCatcher somatic fusion detection from RNA-seq",
        version_args=(),
        docs_url="https://github.com/ndaniel/fusioncatcher",
    ),
    "easyfuse": ToolSpec(
        name="easyfuse",
        category="fusion",
        executable="easyfuse-neoag",
        description="EasyFuse fusion metacaller (STAR-Fusion + FusionCatcher + Arriba + ML)",
        version_args=(),
        docs_url="https://github.com/TRON-Bioinformatics/EasyFuse",
    ),
    "gatk": ToolSpec(
        name="gatk",
        category="variant",
        executable="gatk",
        description="GATK somatic variant calling / filtering",
        docs_url="https://gatk.broadinstitute.org/",
    ),
    "pyclone": ToolSpec(
        name="pyclone",
        category="clonality",
        executable="pyclone",
        description="PyClone-VI / PyClone clonal cluster CCF",
        version_args=("--version",),
        fixture_outputs={"ccf": ""},
    ),
}


CommandBuilder = Callable[["RunContext", Path], list[str]]


@dataclass
class RunContext:
    sample_id: str
    outdir: Path
    stub: bool = False
    executables: dict[str, str] = field(default_factory=dict)
    hla_alleles: list[str] = field(default_factory=list)
    raw_peptides: Path | None = None
    tumor_vcf: Path | None = None
    normal_vcf: Path | None = None
    tumor_sample_name: str | None = None
    normal_sample_name: str | None = None
    phased_vcf: Path | None = None
    prediction_algorithms: str = "NetMHCpan"
    fusion_tsv: Path | None = None
    splice_junction_tsv: Path | None = None
    reference_fasta: Path | None = None
    gencode_gtf: Path | None = None
    expression_tsv: Path | None = None
    variants_vcf: Path | None = None
    facets_rds: Path | None = None
    lohhla_prediction: Path | None = None
    pass_only: bool = True
    tool_provenance: dict = field(default_factory=dict)

    def exe(self, tool: str) -> str:
        import os
        if tool == "vep":
            return (
                self.executables.get(tool)
                or os.environ.get("NEOAG_VEP_BIN")
                or TOOL_REGISTRY[tool].executable
            )
        return self.executables.get(tool) or TOOL_REGISTRY[tool].executable
