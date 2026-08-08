"""sim.validate.validation_export -- the payload behind the Phase 5
Validation view. Computed once per test session (it takes real wall time:
two Baruah cases, a density-multiplier bisection, and Swarm C's four-latitude
sensitivity sweep), then checked against the numbers established at Gate 2/3.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim.atmosphere import SpaceWeather
from sim.validate import analytic_validation_summary, validation_export

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


def test_analytic_summary_matches_the_unit_tests():
    """Tests 1-3, recomputed independently of tests/test_validation.py."""
    s = analytic_validation_summary()
    assert s["test_1_energy_conservation"]["passed"]
    assert s["test_1_energy_conservation"]["relative_change_in_a"] < 1e-12
    assert s["test_2_thrust_spiral"]["passed"]
    assert s["test_2_thrust_spiral"]["relative_error"] < 1e-4
    assert s["test_3_critical_density_fixed_point"]["passed"]
    assert s["test_3_critical_density_fixed_point"]["relative_rate"] < 1e-9


@pytest.fixture(scope="module")
def payload(sw) -> dict:
    return validation_export(sw, dt=20.0)


def test_baruah_cases_match_gate_2_numbers(payload):
    cases = {c["ram_area_m2"]: c for c in payload["test_4_baruah_reproduction"]["cases"]}
    big, small = cases[4.48], cases[1.00]

    assert big["outcome"] == "REENTERED"
    assert 15.0 < big["reentry_error_pct"] < 21.0, "outside the 20% acceptance band"

    assert small["outcome"] in ("INDETERMINATE", "REACHED_SHELL", "REENTERED")
    assert abs(small["decay_error_pct"]) < 20.0, "outside the 20% acceptance band"


def test_density_multipliers_agree_to_within_a_few_percent(payload):
    """README.md's central finding: 4.48x apart in area, ~2% apart in correction."""
    k = payload["test_4_baruah_reproduction"]["implied_density_multiplier"]
    k_big, k_small = k["4.48"], k["1.00"]
    assert 1.1 < k_big < 1.3
    assert 1.1 < k_small < 1.3
    spread = abs(k_big - k_small) / ((k_big + k_small) / 2)
    assert spread < 0.05, f"multipliers disagree by {spread * 100:.1f}%, expected ~2%"


def test_swarm_c_is_flagged_weakest(payload):
    swarm = payload["swarm_c_secondary"]
    assert swarm["flagged_weakest"] is True
    assert len(swarm["flagged_reason"]) > 0
    assert set(swarm["latitude_sensitivity"]) == {"0.0", "45.0", "70.0", "87.4"}


def test_payload_is_plain_json(payload):
    text = json.dumps(payload, allow_nan=False)
    assert json.loads(text) == payload


def test_validation_key_lands_in_sweeps_json_shape():
    """Structural check against sim/sweeps.py's payload -- see test_export.py
    for the full read-back gate; this only pins that the key is named right."""
    from sim.export import load_json

    p = Path(__file__).resolve().parent.parent / "out" / "sweeps.json"
    if not p.exists():
        pytest.skip("out/sweeps.json not generated in this environment")
    d = load_json(p)
    assert "validation" in d
    assert set(d["validation"]) == {
        "analytic_tests", "test_4_baruah_reproduction", "swarm_c_secondary",
    }
