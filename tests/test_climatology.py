"""Long-horizon climatology atmosphere. V2_BRIEF.md §3-4.

`ClimatologyDensity` is the model that makes 25-year compliance runs possible.
It is also the model most able to launder an assumption into a result, because
it replaces real space weather with a stated level, so these tests concentrate
on the places where that could go wrong quietly.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

from sim.atmosphere import (
    SOLAR_ACTIVITY_PERCENTILES,
    ClimatologyDensity,
    SpaceWeather,
    density_from_indices,
)

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"
LAT = 53.22


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


@pytest.fixture(scope="module")
def levels(sw) -> dict[str, ClimatologyDensity]:
    return {
        name: ClimatologyDensity.for_level(name, sw, lat_deg=LAT)
        for name in SOLAR_ACTIVITY_PERCENTILES
    }


def test_density_falls_monotonically_with_altitude(levels):
    for name, c in levels.items():
        assert np.all(np.diff(c.rho_profile) < 0.0), (
            f"{name} profile is not monotonically decreasing in altitude"
        )


def test_activity_levels_are_ordered_everywhere(levels):
    """high > mean > low at every altitude, or the levels are mislabelled."""
    lo, mid, hi = levels["low"], levels["mean"], levels["high"]
    assert np.all(hi.rho_profile > mid.rho_profile)
    assert np.all(mid.rho_profile > lo.rho_profile)


def test_solar_cycle_spread_is_large_at_altitude(levels):
    """The spread is the headline uncertainty, so pin its size.

    An order of magnitude between solar min and max at 700 km is the
    well-known behaviour of the thermosphere. If this collapsed, the swept
    band would be reporting false precision.
    """
    lo, hi = levels["low"], levels["high"]
    ratio_700 = hi(0.0, 700e3) / lo(0.0, 700e3)
    ratio_400 = hi(0.0, 400e3) / lo(0.0, 400e3)
    assert ratio_700 > 10.0, f"only {ratio_700:.1f}x spread at 700 km"
    assert ratio_400 > 5.0, f"only {ratio_400:.1f}x spread at 400 km"
    assert ratio_700 > ratio_400, "spread should widen with altitude"


def test_interpolation_error_is_under_one_percent(levels):
    """PHYSICS.md §6.3's bar, applied to the climatology table."""
    for name, c in levels.items():
        err = c.max_relative_error()
        assert err < 0.01, f"{name}: interpolation error {err * 100:.3f}%"


def test_levels_are_read_from_the_data_not_hardcoded(sw, levels):
    """The activity levels must trace to SW-All.csv, like every other number."""
    for name, c in levels.items():
        p = SOLAR_ACTIVITY_PERCENTILES[name]
        assert c.f107 == pytest.approx(sw.f107_percentile(p))
        assert c.ap == pytest.approx(sw.ap_percentile(p))
    # And they must actually differ, or the percentiles are not being applied.
    assert levels["low"].f107 < levels["mean"].f107 < levels["high"].f107


def test_reference_year_does_not_change_the_profile(sw):
    """Twelve-month averaging should make the sampled year immaterial.

    If this fails, the average is not covering the seasonal cycle and the
    profile carries a hidden date dependence.
    """
    a = ClimatologyDensity.for_level("mean", sw, lat_deg=LAT, reference_year=1996)
    b = ClimatologyDensity.for_level("mean", sw, lat_deg=LAT, reference_year=2014)
    rel = np.max(np.abs(a.rho_profile - b.rho_profile) / b.rho_profile)
    assert rel < 0.02, f"profile moved {rel * 100:.2f}% between reference years"


def test_matches_direct_pymsis_at_the_same_indices(sw, levels):
    """The averaging must not drift from what pymsis returns for one sample.

    A single point-in-time call sits inside the seasonal and local-time spread
    the average removes, so this is a band check, not an equality check --
    but a factor-of-two band still catches a units error, a wrong percentile
    or a botched reshape, which is what it is for.
    """
    c = levels["mean"]
    alts_km = np.array([300.0, 500.0, 700.0])
    direct = density_from_indices(
        [datetime(2000, 6, 15, 12, 0)],
        alts_km,
        LAT,
        0.0,
        [c.f107],
        [c.f107a],
        [[c.ap] * 7],
    )[0]
    averaged = c.lookup(alts_km * 1e3)
    ratio = averaged / direct
    assert np.all((ratio > 0.5) & (ratio < 2.0)), (
        f"climatology/pymsis ratio {ratio} outside the 0.5-2.0 band"
    )


def test_is_interchangeable_with_density_grid(levels):
    """Must satisfy the same (t_s, h_m) -> float contract as DensityGrid.

    The propagator's `deriv` closure calls one or the other; a signature
    mismatch would only show up at runtime deep inside an integration.
    """
    c = levels["mean"]
    value = c(0.0, 400e3)
    assert isinstance(value, float)
    assert value > 0.0
    # Time is ignored by construction -- assert that, so it cannot start
    # mattering silently.
    assert c(0.0, 400e3) == c(1e9, 400e3)


def test_altitude_is_clamped_outside_the_table(levels):
    """RK4 interior stages can leave the table; extrapolation must not diverge."""
    c = levels["mean"]
    assert c(0.0, -50e3) == pytest.approx(c(0.0, c.alts_km[0] * 1e3))
    assert c(0.0, 5000e3) == pytest.approx(c(0.0, c.alts_km[-1] * 1e3))


def test_table_covers_the_leo_region_the_disposal_rule_defines(levels):
    """47 CFR 25.283(e) scopes LEO as below 2000 km; the table must reach it."""
    for c in levels.values():
        assert c.alts_km[-1] >= 2000.0


def test_unknown_level_is_rejected(sw):
    with pytest.raises(ValueError, match="unknown solar activity level"):
        ClimatologyDensity.for_level("catastrophic", sw, lat_deg=LAT)


def test_space_weather_file_really_does_run_out(sw):
    """The premise of this whole module, asserted rather than assumed.

    If the data file ever gains decades of forecast, the climatology's
    justification changes and this test should be the thing that notices.
    """
    horizon_years = (
        datetime(sw.last_day.year, sw.last_day.month, sw.last_day.day,
                 tzinfo=timezone.utc)
        - datetime.now(timezone.utc)
    ).days / 365.25
    assert horizon_years < 25.0, (
        f"SW-All.csv now reaches {horizon_years:.1f} years out; a 25-year run "
        "may no longer need the climatology model"
    )
