"""Generalised mission entry point. V2_BRIEF.md §7, Phase 7.

The riskiest thing this module does is choose an atmosphere model on the
caller's behalf, so most of these tests are about that choice being correct,
visible, and reported.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from sim.atmosphere import SpaceWeather
from sim.disposal import Compliance
from sim.mission import Mission, choose_atmosphere, fly, fly_solar_band
from sim.satellite import Outcome, SatelliteConfig, ThrusterMode

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"
DAY_YEARS = 1.0 / 365.25


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


def a_satellite(**kw) -> SatelliteConfig:
    base = dict(mass_kg=260.0, area_m2=1.0, cd=2.2, thrust_n=0.0,
                isp_s=1666.0, dry_mass_kg=220.0)
    base.update(kw)
    return SatelliteConfig(**base)


def a_mission(sw=None, **kw) -> Mission:
    base = dict(
        epoch=datetime(2027, 3, 1, tzinfo=timezone.utc),
        insertion_altitude_m=705e3,
        inclination_deg=53.0,
        satellite=a_satellite(),
        mission_duration_years=7.0,
        thruster_mode=ThrusterMode.SAFE_MODE,
        label="test",
    )
    base.update(kw)
    return Mission(**base)


# --------------------------------------------------------------------------
# Mission validation
# --------------------------------------------------------------------------

def test_naive_epoch_is_coerced_to_utc():
    m = a_mission(epoch=datetime(2027, 3, 1))
    assert m.epoch.tzinfo is timezone.utc


def test_retrograde_inclination_is_rejected_with_guidance():
    with pytest.raises(ValueError, match="180 - i"):
        a_mission(inclination_deg=97.4)


def test_target_below_insertion_is_rejected():
    with pytest.raises(ValueError, match="lowering is a disposal manoeuvre"):
        a_mission(insertion_altitude_m=500e3, target_altitude_m=400e3)


def test_safe_mode_forces_zero_thrust():
    """PHYSICS.md §5: safe mode is F = 0 regardless of the rated thrust."""
    sat = a_satellite(thrust_n=0.071)
    assert a_mission(satellite=sat, thruster_mode=ThrusterMode.SAFE_MODE).thrust_n == 0.0
    assert a_mission(satellite=sat, thruster_mode=ThrusterMode.NOMINAL).thrust_n == 0.071


# --------------------------------------------------------------------------
# Atmosphere selection -- the judgement call
# --------------------------------------------------------------------------

def test_short_historical_window_uses_real_space_weather(sw):
    m = a_mission(
        epoch=datetime(2022, 2, 3, 18, 13, tzinfo=timezone.utc),
        mission_duration_years=3 * DAY_YEARS,
    )
    _, choice = choose_atmosphere(m, sw)
    assert choice.model == "DensityGrid"
    assert choice.uses_real_space_weather
    assert choice.solar_activity is None


def test_window_past_the_end_of_the_file_falls_back_to_climatology(sw):
    m = a_mission(epoch=datetime(2027, 3, 1, tzinfo=timezone.utc))
    _, choice = choose_atmosphere(m, sw)
    assert choice.model == "ClimatologyDensity"
    assert not choice.uses_real_space_weather
    assert choice.solar_activity == "mean"
    assert str(sw.last_day) in choice.reason, (
        "the reason should name the date the data runs out"
    )


def test_long_window_inside_the_file_is_still_too_expensive(sw):
    """Coverage is necessary but not sufficient -- cost rules out long grids.

    A two-year window ending before the file does would be affordable to
    build only as ~5800 three-hour blocks.
    """
    start = datetime(2015, 1, 1, tzinfo=timezone.utc)
    m = a_mission(epoch=start, mission_duration_years=2.0)
    _, choice = choose_atmosphere(m, sw)
    assert choice.model == "ClimatologyDensity"
    assert "exceeds" in choice.reason


def test_atmosphere_choice_is_reported_in_the_result(sw):
    """A silently swapped atmosphere model would be the worst failure here."""
    r = fly(a_mission(), sw)
    d = r.as_dict()
    assert d["atmosphere"]["model"] == "ClimatologyDensity"
    assert d["atmosphere"]["reason"]
    assert d["atmosphere"]["uses_real_space_weather"] is False


# --------------------------------------------------------------------------
# The V1 case, through the V2 entry point
# --------------------------------------------------------------------------

def test_v1_baruah_case_reproduces_through_the_mission_api(sw, capsys):
    """End-to-end: the validated V1 result, reached via the generalised path.

    Nothing about V1 is special-cased in `mission.py`, so this checks that
    generalising the engine did not quietly change what it computes for the
    one case where the answer is known.
    """
    m = Mission(
        epoch=datetime(2022, 2, 3, 18, 13, tzinfo=timezone.utc),
        insertion_altitude_m=210e3,
        inclination_deg=53.22,
        satellite=SatelliteConfig(mass_kg=227.0, area_m2=4.48, cd=1.0,
                                  thrust_n=0.0, isp_s=None),
        mission_duration_years=3 * DAY_YEARS,
        thruster_mode=ThrusterMode.SAFE_MODE,
        label="baruah-4.48",
    )
    r = fly(m, sw, storm_time=True, dt_max=1800.0)

    assert r.atmosphere.uses_real_space_weather
    assert r.outcome is Outcome.REENTERED

    v1_fixed_step_s = 164900.0
    err = abs(r.outcome_time_s - v1_fixed_step_s) / v1_fixed_step_s
    with capsys.disabled():
        print(f"\n  Mission API vs V1 fixed step: {r.outcome_time_s:.1f} s vs "
              f"{v1_fixed_step_s:.1f} s -> {err * 100:.4f}%")
    assert err < 1e-3, f"V2 entry point moved the V1 answer by {err * 100:.4f}%"


# --------------------------------------------------------------------------
# Lifecycle
# --------------------------------------------------------------------------

def test_thrusting_satellite_reaches_its_shell(sw):
    """Arbitrary target: raise from 210 km to a 320 km shell."""
    m = Mission(
        epoch=datetime(2027, 3, 1, tzinfo=timezone.utc),
        insertion_altitude_m=210e3,
        inclination_deg=53.0,
        satellite=a_satellite(thrust_n=0.071, area_m2=1.0),
        mission_duration_years=0.5,
        target_altitude_m=320e3,
        thruster_mode=ThrusterMode.NOMINAL,
        label="raise",
    )
    r = fly(m, sw, assess_disposal=False)
    assert r.outcome is Outcome.REACHED_SHELL
    assert r.outcome_time_s is not None
    assert r.final_altitude_km >= 320.0


def test_high_orbit_gets_a_compliance_verdict(sw):
    r = fly(a_mission(), sw)
    assert r.outcome is not Outcome.REENTERED
    assert r.compliance is not None
    assert r.compliance.verdict in {
        Compliance.COMPLIANT_WITH_DISPOSAL,
        Compliance.NON_COMPLIANT_INSUFFICIENT_PROPELLANT,
    }
    assert r.compliance.delta_v_ms > 0.0


def test_disposal_is_priced_from_the_end_of_life_state(sw):
    """Not from the nominal shell -- the satellite has drifted and burned mass.

    A satellite that decayed during its mission is cheaper to dispose of than
    one that held station, and the verdict must reflect where it actually is.
    """
    r = fly(a_mission(), sw, solar_activity="high")
    assert r.compliance is not None
    # It sank during the mission, so disposal starts below the insertion
    # altitude, not at 705 km.
    assert r.compliance.operational_altitude_km < 705.0
    assert r.compliance.operational_altitude_km == pytest.approx(
        r.final_altitude_km, rel=1e-9
    )


def test_solar_band_sweeps_rather_than_choosing(sw):
    """V2_BRIEF.md §6: uncertain values get swept. The spread is the answer."""
    band = fly_solar_band(a_mission(), sw)
    assert set(band) == {"low", "mean", "high"}

    targets = {k: v.compliance.disposal_altitude_km for k, v in band.items()}
    assert targets["low"] < targets["mean"] < targets["high"], (
        "a denser thermosphere should let disposal stop higher up"
    )
    costs = {k: v.compliance.delta_v_ms for k, v in band.items()}
    assert costs["high"] < costs["mean"] < costs["low"]


def test_result_serialises_whole(sw):
    d = fly(a_mission(), sw).as_dict()
    for key in ("mission", "outcome", "final_altitude_km", "atmosphere",
                "step_stats", "compliance"):
        assert key in d
    assert d["step_stats"]["tolerance_respected"] is True
    assert d["compliance"]["rule"]["citation"] == "47 CFR § 25.283(e)"
