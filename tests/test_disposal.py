"""Disposal delta-v, propellant, and the compliance verdict. V2_BRIEF.md §3.

The delta-v maths is closed-form, so it is tested against hand calculations
rather than against itself. The compliance verdict is tested for the
properties that must hold rather than for the specific altitudes it currently
returns, which depend on the atmosphere model and are expected to move.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from sim.atmosphere import ClimatologyDensity, SpaceWeather
from sim.constants import G0, MU, R_E
from sim.disposal import (
    SECONDS_PER_YEAR,
    Compliance,
    assess_compliance,
    circular_velocity,
    decay_time_s,
    default_rule,
    highest_complying_altitude_m,
    hohmann_lower_circular,
    load_rules,
    perigee_lowering_dv,
    propellant_mass,
    vis_viva,
)

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"
LAT = 53.0

# Bc = m/(Cd*A) = 100 kg/m^2, the ballistic coefficient V2_BRIEF.md §3's
# altitude table assumes, so results here are comparable to it.
MASS, CD, AREA = 220.0, 2.2, 1.0
ISP = 1666.0


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


@pytest.fixture(scope="module")
def atmos(sw):
    return {
        lvl: ClimatologyDensity.for_level(lvl, sw, lat_deg=LAT)
        for lvl in ("low", "mean", "high")
    }


# --------------------------------------------------------------------------
# Closed-form orbital mechanics
# --------------------------------------------------------------------------

def test_vis_viva_reduces_to_circular_velocity():
    """At r = a the ellipse is a circle. The two formulas must agree."""
    r = R_E + 550e3
    assert vis_viva(r, r) == pytest.approx(circular_velocity(r), rel=1e-12)


def test_hohmann_matches_hand_calculation():
    """**Phase 7 gate:** disposal delta-v checked against a hand calculation.

    Worked by hand for 705 km -> 400 km circular, MU = 3.986004418e14,
    R_E = 6378.137 km:

        r1  = 7083.137 km,  r2 = 6778.137 km,  a_t = (r1+r2)/2 = 6930.637 km
        v1  = sqrt(MU/r1)              = 7501.6374 m/s
        vt1 = sqrt(MU*(2/r1 - 1/a_t))  = 7418.6462 m/s
        dv1 = v1 - vt1                 =   82.9911 m/s
        v2  = sqrt(MU/r2)              = 7668.5582 m/s
        vt2 = sqrt(MU*(2/r2 - 1/a_t))  = 7752.4676 m/s
        dv2 = vt2 - v2                 =   83.9094 m/s
        total                          =  166.9006 m/s
    """
    dv1, dv2 = hohmann_lower_circular(R_E + 705e3, R_E + 400e3)
    assert dv1 == pytest.approx(82.9911, abs=1e-3)
    assert dv2 == pytest.approx(83.9094, abs=1e-3)
    assert dv1 + dv2 == pytest.approx(166.9006, abs=1e-3)


def test_hohmann_agrees_with_the_small_step_approximation():
    """Independent magnitude check: for |dr| << r, total dv -> (v/2)*|dr|/r.

    A different derivation reaching the same number, which is what makes it
    worth having alongside the hand calculation -- an algebra slip in
    `hohmann_lower_circular` would have to be reproduced in the first-order
    expansion too in order to hide here.
    """
    r1, r2 = R_E + 705e3, R_E + 400e3
    exact = sum(hohmann_lower_circular(r1, r2))
    approx = 0.5 * circular_velocity(r1) * abs(r1 - r2) / r1
    # First order low by ~3% at |dr|/r ~ 0.04, and low is the expected sign.
    assert approx < exact
    assert abs(exact - approx) / exact < 0.05


def test_hohmann_is_symmetric_in_direction():
    """Raising costs exactly what lowering costs. Magnitudes, so no sign flip."""
    r1, r2 = R_E + 400e3, R_E + 705e3
    assert sum(hohmann_lower_circular(r1, r2)) == pytest.approx(
        sum(hohmann_lower_circular(r2, r1)), rel=1e-12
    )


def test_hohmann_to_the_same_orbit_is_free():
    dv1, dv2 = hohmann_lower_circular(R_E + 500e3, R_E + 500e3)
    assert dv1 == pytest.approx(0.0, abs=1e-9)
    assert dv2 == pytest.approx(0.0, abs=1e-9)


def test_perigee_lowering_is_cheaper_than_the_circular_transfer():
    """The module's central caveat, asserted rather than just documented.

    If this ever inverted, the "conservative bound" claim in the module
    docstring and in `margin_note` would be false.
    """
    r_op = R_E + 705e3
    circular = sum(hohmann_lower_circular(r_op, R_E + 400e3))
    elliptical = perigee_lowering_dv(r_op, R_E + 200e3)
    assert elliptical < circular


def test_perigee_lowering_rejects_a_raise():
    with pytest.raises(ValueError, match="that is a raise"):
        perigee_lowering_dv(R_E + 400e3, R_E + 700e3)


def test_rocket_equation_against_hand_calculation():
    """m_prop = m0*(1 - exp(-dv/(Isp*g0))), checked term by term."""
    dv, m0, isp = 166.9006, 220.0, 1666.0
    expected = m0 * (1.0 - math.exp(-dv / (isp * G0)))
    assert propellant_mass(dv, m0, isp) == pytest.approx(expected, rel=1e-12)
    assert propellant_mass(dv, m0, isp) == pytest.approx(2.2360, abs=1e-3)


def test_rocket_equation_limits():
    assert propellant_mass(0.0, 220.0, ISP) == pytest.approx(0.0)

    # Monotone in delta-v, and asymptotic to the vehicle mass from below.
    m0 = 220.0
    a = propellant_mass(1e4, m0, ISP)
    b = propellant_mass(3e4, m0, ISP)
    assert 0.0 < a < b < m0

    # At absurd delta-v the exponential underflows and the result is exactly
    # m0. That is the correct limit -- propellant can never exceed the vehicle
    # -- so the invariant is "never more than m0", not "strictly less".
    assert propellant_mass(1e9, m0, ISP) == m0


def test_rocket_equation_rejects_bad_inputs():
    with pytest.raises(ValueError, match="Isp must be positive"):
        propellant_mass(100.0, 220.0, 0.0)
    with pytest.raises(ValueError, match="mass must be positive"):
        propellant_mass(100.0, -1.0, ISP)
    with pytest.raises(ValueError, match="delta-v must be non-negative"):
        propellant_mass(-1.0, 220.0, ISP)


# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

def test_rules_load_with_their_citations():
    rules = load_rules()
    fcc = rules["fcc_5yr"]
    assert fcc.window_years == 5.0
    assert "25.283" in fcc.citation
    assert fcc.clock_starts == "end_of_mission"
    assert fcc.applies_below_altitude_km == 2000.0
    assert fcc.quote, "the rule text itself should travel with the rule"
    assert rules["iadc_25yr"].window_years == 25.0


def test_default_rule_is_the_fcc_rule():
    assert default_rule().id == "fcc_5yr"


def test_rule_scope_is_altitude_bounded():
    fcc = default_rule()
    assert fcc.applies_at(700.0)
    assert not fcc.applies_at(2500.0)


# --------------------------------------------------------------------------
# Decay time and the complying altitude
# --------------------------------------------------------------------------

def test_decay_time_increases_with_altitude(atmos):
    c = atmos["mean"]
    times = [
        decay_time_s(c, h, MASS, CD, AREA, t_max_s=200 * SECONDS_PER_YEAR)
        for h in (300e3, 400e3, 500e3)
    ]
    assert all(t is not None for t in times)
    assert times[0] < times[1] < times[2]


def test_decay_time_decreases_with_solar_activity(atmos):
    """More solar activity, denser thermosphere, faster decay."""
    t = {
        lvl: decay_time_s(atmos[lvl], 450e3, MASS, CD, AREA,
                          t_max_s=200 * SECONDS_PER_YEAR)
        for lvl in ("low", "mean", "high")
    }
    assert t["high"] < t["mean"] < t["low"]


def test_highest_complying_altitude_brackets_the_window(atmos):
    """The returned altitude complies; a little above it does not.

    This is the property that makes the bisection result meaningful, and it
    does not depend on the specific altitude, so it survives an atmosphere
    model change.
    """
    c = atmos["mean"]
    window = default_rule().window_s
    h = highest_complying_altitude_m(c, MASS, CD, AREA, window)
    assert h is not None

    below = decay_time_s(c, h - 5e3, MASS, CD, AREA, t_max_s=window)
    assert below is not None and below <= window

    above = decay_time_s(c, h + 20e3, MASS, CD, AREA, t_max_s=window)
    assert above is None or above > window


def test_the_briefs_550km_figure_is_the_25_year_threshold(atmos):
    """V2_BRIEF.md §3 says "below ~550 km" natural decay is the question.

    Computing it shows 550 km is the *25-year* threshold, not the FCC
    five-year one, which lands near 450 km at mean solar activity. The brief
    told us to verify its regulatory paragraph rather than trust it; this is
    the check, kept as a test so the distinction cannot quietly rot.

    Bands are wide because the point is which rule each figure belongs to,
    not the precise altitude.
    """
    c = atmos["mean"]
    rules = load_rules()
    h_5yr = highest_complying_altitude_m(c, MASS, CD, AREA, rules["fcc_5yr"].window_s)
    h_25yr = highest_complying_altitude_m(c, MASS, CD, AREA, rules["iadc_25yr"].window_s)

    assert h_5yr < h_25yr, "a longer window must permit a higher orbit"
    assert 400e3 < h_5yr < 500e3, f"5-year threshold at {h_5yr / 1e3:.0f} km"
    assert 500e3 < h_25yr < 600e3, f"25-year threshold at {h_25yr / 1e3:.0f} km"


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

def test_low_orbit_complies_naturally(atmos):
    res = assess_compliance(
        atmos["mean"], 350e3, MASS, CD, AREA, ISP,
        propellant_available_kg=5.0, solar_activity="mean",
    )
    assert res.verdict is Compliance.COMPLIANT_NATURAL
    assert res.delta_v_ms == 0.0
    assert res.propellant_required_kg == 0.0
    assert res.natural_decay_years is not None
    assert res.natural_decay_years < res.window_years


def test_high_orbit_needs_a_burn_it_can_afford(atmos):
    res = assess_compliance(
        atmos["mean"], 705e3, MASS, CD, AREA, ISP,
        propellant_available_kg=5.0, solar_activity="mean",
    )
    assert res.verdict is Compliance.COMPLIANT_WITH_DISPOSAL
    assert res.natural_decay_years is None
    assert res.disposal_altitude_km < 705.0
    assert res.delta_v_ms > 0.0
    assert 0.0 < res.propellant_required_kg <= res.propellant_available_kg
    assert res.margin_note, "the circular-orbit caveat must travel with the result"


def test_insufficient_propellant_is_caught(atmos):
    """The failure V2_BRIEF.md §3 says this tool exists to catch."""
    res = assess_compliance(
        atmos["mean"], 705e3, MASS, CD, AREA, ISP,
        propellant_available_kg=0.1, solar_activity="mean",
    )
    assert res.verdict is Compliance.NON_COMPLIANT_INSUFFICIENT_PROPELLANT
    assert res.propellant_required_kg > res.propellant_available_kg
    assert any("short by" in n for n in res.notes)


def test_above_the_leo_region_is_out_of_scope(atmos):
    """47 CFR 25.283(e) stops at 2000 km; the verdict must say so, not guess."""
    res = assess_compliance(
        atmos["mean"], 2500e3, MASS, CD, AREA, ISP,
        propellant_available_kg=50.0, solar_activity="mean",
    )
    assert res.verdict is Compliance.OUT_OF_SCOPE
    assert res.delta_v_ms is None


def test_solar_activity_changes_the_disposal_cost(atmos):
    """The band is the answer, so the levels must actually separate.

    A denser thermosphere means disposal can stop higher up and costs less.
    """
    results = {
        lvl: assess_compliance(
            atmos[lvl], 705e3, MASS, CD, AREA, ISP,
            propellant_available_kg=5.0, solar_activity=lvl,
        )
        for lvl in ("low", "mean", "high")
    }
    assert (
        results["high"].disposal_altitude_km
        > results["mean"].disposal_altitude_km
        > results["low"].disposal_altitude_km
    )
    assert (
        results["high"].delta_v_ms
        < results["mean"].delta_v_ms
        < results["low"].delta_v_ms
    )


def test_result_serialises_for_the_json_export(atmos):
    """Phase 11 will need this; the contract is that nothing is lost."""
    res = assess_compliance(
        atmos["mean"], 705e3, MASS, CD, AREA, ISP,
        propellant_available_kg=5.0, solar_activity="mean",
    )
    d = res.as_dict()
    assert d["verdict"] == "COMPLIANT_WITH_DISPOSAL"
    assert d["rule"]["citation"] == "47 CFR § 25.283(e)"
    assert d["perigee_lowering_delta_v_ms"] > 0.0
    assert d["margin_note"]
