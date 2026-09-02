from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_lohhla_copy_number_table_has_r_read_table_rowname_shape():
    script = (ROOT / "scripts/run_hla_loh_multi_tool.sh").read_text(encoding="utf-8")
    assert "Ploidy\\ttumorPurity\\ttumorPloidy\\n" in script
    assert "'$TUMOR_ID'" not in script
    assert '"$TUMOR_ID" "$PLOIDY" "$PURITY" "$PLOIDY"' in script
    assert "printf '\\ttumorPurity\\ttumorPloidy\\n'" not in script


def test_lohhla_runtime_is_resolved_as_one_noninteractive_toolchain():
    script = (ROOT / "scripts/run_lohhla_sample.sh").read_text(encoding="utf-8")
    for package in ("optparse", "Rsamtools", "Biostrings", "seqinr"):
        assert f'requireNamespace("{package}",quietly=TRUE)' in script
    assert "resolve_bedtools_dir" in script
    assert "LOHHLA_BEDTOOLS_DIR" in script


def test_lohhla_discovers_polysolver_from_portable_asset_roots():
    runner = (ROOT / "scripts/run_lohhla_sample.sh").read_text(encoding="utf-8")
    tools_env = (ROOT / "conf/tools.env.sh").read_text(encoding="utf-8")
    for text in (runner, tools_env):
        assert "NEOAG_ASSET_ROOT" in text
        assert "data/lohhla/polysolver" in text
        assert "scripts/shell_call_hla_type" in text
        assert "data/abc_complete.fasta" in text


def test_lohhla_uses_complete_picard_command_jar_runtime():
    script = (ROOT / "scripts/run_lohhla_sample.sh").read_text(encoding="utf-8")
    for jar in ("SamToFastq.jar", "SortSam.jar", "FilterSamReads.jar"):
        assert f'${{PSHOME}}/binaries/{jar}' in script
    assert 'LOHHLA_GATK_RUNTIME_DIR="${PSHOME}/binaries"' in script


def test_easyfuse_patches_fusioncatcher_for_every_layout():
    script = (ROOT / "scripts/run_easyfuse_sample.sh").read_text(encoding="utf-8")
    module_native = script.index("EasyFuse module-native layout")
    run_workflow = script.index("run_nextflow()")
    preflight = script[module_native:run_workflow]
    assert 'bash "${ROOT}/scripts/patch_easyfuse_fusioncatcher_compat.sh"' in preflight
    retry = script[script.index("run_nextflow || {"):]
    assert 'bash "${ROOT}/scripts/patch_easyfuse_fusioncatcher_compat.sh"' in retry
