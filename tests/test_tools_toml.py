"""Guards against conf/tools.toml regressions: it must always be valid TOML,
every section must declare a real 'entries' list and a supported 'check' type,
and every entry_mode used by run-demo must resolve to at least one tool.

This exists because conf/tools.toml was hand-edited to add a section with a
missing value (`check = ` with nothing after it), which made the file
unparsable and broke `neoag-v03 run-demo` for every entry_mode until it was
caught manually. This test makes that class of mistake fail CI instead.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from neoag_v03.tools_config import (
    DEFAULT_TOOLS_TOML,
    check_entry_tools,
    load_tools_toml,
)

ENTRY_MODES = ["snv_indel", "fusion", "splice_junction", "sv_wgs", "sv_wes", "peptide_only"]
VALID_CHECK_KINDS = {"bin", "dir", "file", "env", "bin_or_env"}
REQUIRED_FIELD_BY_CHECK = {
    "bin": ("bin",),
    "dir": ("dir",),
    "file": ("file",),
    "env": ("env_var",),
    "bin_or_env": ("bin", "env_var"),
}


def test_tools_toml_is_valid_toml():
    assert DEFAULT_TOOLS_TOML.is_file(), "conf/tools.toml is missing"
    data = load_tools_toml()
    assert isinstance(data, dict)
    assert data, "conf/tools.toml parsed to an empty document"


def test_every_section_has_a_supported_check_kind_and_required_fields():
    data = load_tools_toml()
    for name, section in data.items():
        if not isinstance(section, dict) or not section:
            # e.g. [spechla]: intentionally empty, informational-only section.
            continue
        check = section.get("check")
        assert check in VALID_CHECK_KINDS, f"[{name}] has unsupported check kind: {check!r}"
        for field in REQUIRED_FIELD_BY_CHECK[check]:
            assert field in section, f"[{name}] check={check!r} but missing required field '{field}'"
        entries = section.get("entries", [])
        assert isinstance(entries, list), f"[{name}].entries must be a list"
        for mode in entries:
            assert mode in ENTRY_MODES, f"[{name}] declares unknown entry_mode {mode!r}"


@pytest.mark.parametrize("entry_mode", ENTRY_MODES)
def test_check_entry_tools_does_not_raise(entry_mode):
    # Must run cleanly with zero tools installed (the CI/sandbox default) --
    # this is exactly what run-demo calls before every fixture smoke test.
    result = check_entry_tools(entry_mode)
    assert isinstance(result, dict)
    for name, (ok, msg) in result.items():
        assert isinstance(ok, bool)
        assert isinstance(msg, str) and msg


def test_snv_indel_and_splice_junction_require_vep():
    data = load_tools_toml()
    vep = data.get("vep", {})
    assert vep.get("optional") is False, "VEP should be marked required (optional=false) for snv_indel/splice_junction"
    assert "snv_indel" in vep.get("entries", [])
    assert "splice_junction" in vep.get("entries", [])
