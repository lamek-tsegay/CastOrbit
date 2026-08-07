"""Regression tests for validation Test 4 (PHYSICS.md §8).

These pin *what the model produces*, not what the paper says. The published
targets are compared against and reported, but the assertions here are
deliberately loose bands around the observed output, so that a future change to
the atmosphere, integrator or grid that moves the answer shows up as a failure
rather than passing silently.

Marked slow: each case propagates ~2 days at 10 s.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sim.atmosphere import SpaceWeather
from sim.validate import (
    CD_VALIDATION,
    MASS_KG,
    TARGETS,
    THRUST_N,
    T_REFERENCE_S,
    run_case,
)

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


def test_paper_parameters_are_not_drifting():
    """The Cd trap: validation must run at 1.0, not 2.2 (PHYSICS.md §8)."""
    assert CD_VALIDATION == 1.0
    assert MASS_KG == 227.0
    assert THRUST_N == 0.0, "the paper models drag alone; safe mode means F = 0"
    assert T_REFERENCE_S == 139500.0  # 2022-02-03 18:13 -> 2022-02-05 08:58 UT


def test_4_48_case_reenters_late_by_about_seven_hours(sw, capsys):
    """A = 4.48 m2: reenters, ~7 h later than the published 08:58 UT."""
    r = run_case(4.48, storm_time=True, sw=sw, t_max_s=3 * 86400.0)
    assert r.outcome == "REENTERED"
    t_reentry = r.reentry_time_s
    err = (t_reentry - T_REFERENCE_S) / T_REFERENCE_S

    with capsys.disabled():
        print(f"\n  Test 4 (A=4.48 m2): reentry at {t_reentry / 3600:.2f} h, "
              f"target {T_REFERENCE_S / 3600:.2f} h, decay-timing error {err * 100:+.1f}%")

    assert 44.0 < t_reentry / 3600 < 48.0, "reentry time moved outside the observed band"
    assert err < 0.20, "decay timing now exceeds the 20% acceptance in PHYSICS.md §8"


def test_1_00_case_decays_about_five_and_a_half_km(sw, capsys):
    """A = 1.00 m2: survives, decaying ~5.5 km against the published 6.76 km."""
    r = run_case(1.00, storm_time=True, sw=sw, t_max_s=2 * 86400.0)
    alt = r.altitude_at_reference_km
    assert alt is not None
    decay = 210.0 - alt
    target_decay = 210.0 - TARGETS[1.00]
    err = (decay - target_decay) / target_decay

    with capsys.disabled():
        print(f"  Test 4 (A=1.00 m2): altitude {alt:.2f} km, target "
              f"{TARGETS[1.00]:.2f} km; decay {decay:.2f} km vs {target_decay:.2f} km "
              f"({err * 100:+.1f}%)")

    assert 203.5 < alt < 205.5, "altitude at reference moved outside the observed band"
    assert abs(err) < 0.20, "decay now exceeds the 20% acceptance in PHYSICS.md §8"


def test_storm_time_changes_timing_by_under_an_hour(sw, capsys):
    """3-hourly ap vs daily Ap barely moves the integrated result.

    Worth pinning: the §6.2 history machinery is correct but, for this event,
    nearly irrelevant, because AP_AVG is by construction the mean of AP1-AP8
    and the decay integrates over the spikes.
    """
    on = run_case(4.48, storm_time=True, sw=sw, t_max_s=3 * 86400.0)
    off = run_case(4.48, storm_time=False, sw=sw, t_max_s=3 * 86400.0)
    shift_h = (on.reentry_time_s - off.reentry_time_s) / 3600.0

    with capsys.disabled():
        print(f"  Test 4 (storm_time on vs off, A=4.48 m2): "
              f"{on.reentry_time_s / 3600:.2f} h vs {off.reentry_time_s / 3600:.2f} h "
              f"-> {shift_h:+.2f} h")

    assert abs(shift_h) < 1.0, "3-hourly ap now shifts reentry by more than an hour"
    assert shift_h < 0.0, "3-hourly ap should give slightly faster decay here"
