"""Input adapters for the formal Splice Provenance Layer."""

from .easyquant import parse_easyquant
from .immunopepper import parse_immunopepper_kmers, parse_immunopepper_meta
from .irfinder import parse_irfinder
from .k4neo import parse_k4neo
from .mopepgen import parse_mopepgen, parse_mopepgen_gvf
from .pvacbind import parse_pvacbind
from .pvacsplice import parse_pvacsplice
from .regtools import parse_junction_source
from .splice2neo import parse_splice2neo
from .spladder import parse_spladder_gff3, parse_spladder_txt

__all__ = [
    "parse_easyquant",
    "parse_immunopepper_kmers",
    "parse_immunopepper_meta",
    "parse_irfinder",
    "parse_k4neo",
    "parse_mopepgen",
    "parse_mopepgen_gvf",
    "parse_pvacbind",
    "parse_pvacsplice",
    "parse_junction_source",
    "parse_splice2neo",
    "parse_spladder_gff3",
    "parse_spladder_txt",
]
