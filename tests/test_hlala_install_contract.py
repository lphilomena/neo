from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_hlala_installer_uses_real_pinned_bioconda_program():
    script = (
        ROOT
        / ".agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh"
    ).read_text()
    assert 'HLALA_VERSION="1.0.4"' in script
    assert '"hla-la=$HLALA_VERSION"' in script
    assert 'local env_prefix="$home/.conda"' in script
    assert 'ln -sfn "$env_prefix/bin/HLA-LA.pl"' in script
    assert "run_hla_la_container.sh" in script
    assert 'export HLALA_ENV_PREFIX="$TOOLS_ROOT/tools/HLA-LA/.conda"' in script
    assert 'export NEOAG_HLALA_BACKEND="auto"' in script


def test_hlala_production_activation_declares_native_environment():
    script = (
        ROOT
        / ".agents/skills/neoag-remote-deploy/scripts/10_rewrite_production_activation.sh"
    ).read_text()
    for name in (
        "HLALA_ENV_PREFIX",
        "HLA_LA_ENV_PREFIX",
        "HLALA_BIN",
        "HLALA_GRAPH",
        "NEOAG_HLALA_BACKEND",
    ):
        assert name in script
    assert "${HLALA_ENV_PREFIX}/bin/HLA-LA.pl" in script


def test_hlala_verification_requires_executable_runtime_and_prepared_graph():
    script = (ROOT / "scripts/verify_hla_la_container.sh").read_text()
    for value in (
        "bin/HLA-LA.pl",
        "bin/perl",
        "bin/samtools",
        "serializedGRAPH",
        "PRG/graph.txt",
        "--testing 1",
    ):
        assert value in script


def test_hlala_sample_runner_validates_real_bestguess_output():
    script = (ROOT / "scripts/run_hla_la_sample.sh").read_text()
    for value in (
        "--customGraphDir",
        "--maxThreads",
        "R1_bestguess_G.txt",
        "run_metadata.tsv",
        ".complete",
    ):
        assert value in script


def test_hla_typing_skill_documents_real_hlala_execution_contract():
    skill = (ROOT / ".agents/skills/neoag-hla-typing-loh/SKILL.md").read_text()
    for value in (
        "scripts/verify_hla_la_container.sh",
        "scripts/run_hla_la_sample.sh",
        "R1_bestguess_G.txt",
        "run_metadata.tsv",
        ".complete",
    ):
        assert value in skill
