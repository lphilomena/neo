from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location("run_netchop_evidence", ROOT / "scripts/run_netchop_evidence.py")
assert SPEC and SPEC.loader
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_parse_outputs_reports_c_terminal_score(tmp_path: Path):
    output = tmp_path / "netchop.out"
    output.write_text(
        "   1   A  S   0.900000 SEQ1\n"
        "   2   C  .   0.200000 SEQ1\n"
        "   3   D  S   0.700000 SEQ1\n",
        encoding="utf-8",
    )
    maximum, mean, cterm, sites = MODULE.parse_outputs([output])["SEQ1"]
    assert maximum == "0.9"
    assert mean == "0.6"
    assert cterm == "0.7"
    assert sites == "2"
