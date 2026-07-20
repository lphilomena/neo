"""Environment acceptance: all registered tools resolve on a configured host."""

from neoag.tools.runner import TOOL_REGISTRY, check_all_tools


def test_check_all_tools_ok_when_sourced():
    statuses = check_all_tools()
    missing = [s for s in statuses if not s.available]
    assert not missing, "\n".join(
        f"{s.name}: {s.message} ({s.executable})" for s in missing
    )
    assert len(statuses) == len(TOOL_REGISTRY)
