from neoag.controlled_execution.gateway import ROUTE_SPECS, _risk


def test_open_neo_gateway_routes_registered():
    assert "/open/install-check" in ROUTE_SPECS
    assert "/open/run" in ROUTE_SPECS
    assert "/open/review" in ROUTE_SPECS


def test_open_neo_gateway_mode_risk():
    assert _risk("/open/run", {"mode": "plan"}) == "LOW"
    assert _risk("/open/run", {"mode": "execute"}) == "HIGH"
    assert _risk("/open/install-check", {"mode": "verify"}) == "LOW"
    assert _risk("/open/install-check", {"mode": "install"}) == "HIGH"
