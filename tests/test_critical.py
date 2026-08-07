"""Critical density and critical altitude tests.

PHYSICS.md §4.1 supplies a critical density table that must be reproduced, and
a headline result -- the ~10 km margin between the insertion altitude and the
storm-time critical altitude -- that must be *computed*, not hardcoded.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

from sim.atmosphere import SpaceWeather, density_from_indices
from sim.constants import R_E
from sim.critical import critical_altitude, critical_density

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"

# PHYSICS.md §4.1: F = 0.071 N (DERIVED), a = R_E + 210 km.
THRUST = 0.071
A_REF = R_E + 210e3
REF_DATE = datetime(2022, 2, 4, tzinfo=timezone.utc)
REF_LAT, REF_LON = 53.22, 0.0
REF_F107, REF_F107A = 127.0, 110.0

# PHYSICS.md §4.1, "Critical densities". (Cd, A) -> rho_crit
CRITICAL_TABLE = {
    (1.0, 1.00): 2.347e-09,
    (1.0, 4.48): 5.239e-10,
    (2.2, 1.00): 1.067e-09,
    (2.2, 4.48): 2.381e-10,
}
TABLE_RTOL = 5e-4  # the table is quoted to four significant figures


def test_critical_density_table_4_1(capsys):
    """PHYSICS.md §4.1 -- the critical density table, reproduced exactly."""
    with capsys.disabled():
        print("\n  §4.1 critical density table:")
        print("    Cd     A (m2)   computed      published")
    for (cd, area), published in CRITICAL_TABLE.items():
        computed = critical_density(A_REF, THRUST, cd, area)
        with capsys.disabled():
            print(f"    {cd:<6.1f} {area:<8.2f} {computed:.4e}    {published:.3e}")
        assert computed == pytest.approx(published, rel=TABLE_RTOL), (cd, area)


def test_critical_density_scaling():
    """rho_crit is inversely proportional to Cd*A and linear in F and a."""
    base = critical_density(A_REF, THRUST, 1.0, 1.0)
    assert critical_density(A_REF, THRUST, 2.0, 1.0) == pytest.approx(base / 2.0)
    assert critical_density(A_REF, THRUST, 1.0, 2.0) == pytest.approx(base / 2.0)
    assert critical_density(A_REF, 2 * THRUST, 1.0, 1.0) == pytest.approx(2 * base)


@pytest.fixture(scope="module")
def storm_rho():
    """rho(h) for the §4.1 storm conditions, as a callable for bisection."""
    def rho_of_h(h_m: float) -> float:
        return float(
            density_from_indices(
                [REF_DATE], [h_m / 1e3], REF_LAT, REF_LON,
                [REF_F107], [REF_F107A], [[56] * 7],
            )[0, 0]
        )
    return rho_of_h


def test_storm_critical_altitude_is_near_200_km(storm_rho, capsys):
    """PHYSICS.md §4.1 headline -- the margin was about 10 km.

    For the worst-case configuration (Cd = 2.2, A = 4.48 m2) the storm density
    reaches rho_crit at approximately 200 km. The satellites were inserted at
    210 km. The spec is explicit: compute this, do not hardcode it.
    """
    h_crit = critical_altitude(storm_rho, THRUST, cd=2.2, area=4.48)
    margin_km = (210e3 - h_crit) / 1e3

    with capsys.disabled():
        print(
            f"\n  §4.1 headline (computed, not hardcoded):\n"
            f"    worst case Cd=2.2, A=4.48 m2 -> rho_crit = "
            f"{critical_density(R_E + h_crit, THRUST, 2.2, 4.48):.4e} kg/m3\n"
            f"    storm critical altitude      = {h_crit / 1e3:.2f} km\n"
            f"    insertion altitude           = 210.00 km\n"
            f"    margin                       = {margin_km:.2f} km"
        )

    assert 195e3 < h_crit < 205e3, f"h_crit = {h_crit / 1e3:.2f} km, expected ~200 km"
    assert 5.0 < margin_km < 15.0, f"margin {margin_km:.2f} km, expected ~10 km"


def test_critical_altitude_is_a_root(storm_rho):
    """The returned altitude really does satisfy rho(h) == rho_crit(h)."""
    h_crit = critical_altitude(storm_rho, THRUST, cd=2.2, area=4.48, tol_m=1e-4)
    assert storm_rho(h_crit) == pytest.approx(
        critical_density(R_E + h_crit, THRUST, 2.2, 4.48), rel=1e-6
    )


def test_larger_area_raises_critical_altitude(storm_rho):
    """More drag area means the balance point sits higher up -- monotonicity."""
    small = critical_altitude(storm_rho, THRUST, cd=2.2, area=1.00)
    large = critical_altitude(storm_rho, THRUST, cd=2.2, area=4.48)
    assert large > small


def test_critical_altitude_requires_bracketed_root(storm_rho):
    """A configuration with no balance point in [100, 600] km must raise."""
    # Huge thrust, tiny area: thrust wins everywhere above 100 km.
    with pytest.raises(ValueError, match="no critical altitude above"):
        critical_altitude(storm_rho, thrust=50.0, cd=1.0, area=0.01)
    # Negligible thrust, huge area: drag wins everywhere below 600 km.
    with pytest.raises(ValueError, match="no critical altitude below"):
        critical_altitude(storm_rho, thrust=1e-6, cd=2.2, area=100.0)
