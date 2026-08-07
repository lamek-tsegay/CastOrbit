"""Space weather loader and NRLMSIS wrapper tests.

PHYSICS.md §6.2 requires a loader test that reproduces the Feb 3-5 2022 table:
"An off-by-one in the history window is silent and will corrupt every result."
PHYSICS.md §4.1 supplies reference densities that the implementation must
reproduce.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

from sim.atmosphere import SpaceWeather, density_from_indices

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"

# PHYSICS.md §4.1: lat 53.22, 2022-02-04, F10.7 = 127, F10.7A = 110.
REF_DATE = datetime(2022, 2, 4, tzinfo=timezone.utc)
REF_LAT, REF_LON = 53.22, 0.0
REF_F107, REF_F107A = 127.0, 110.0

# PHYSICS.md §6.2, "Verified values for the validation window".
SW_TABLE = {
    "2022-02-03": dict(
        ap_avg=26, ap3=(7, 32, 48, 56, 27, 15, 15, 12), f107=126.5, f107a=109.1
    ),
    "2022-02-04": dict(
        ap_avg=32, ap3=(27, 27, 22, 22, 27, 56, 48, 27), f107=129.6, f107a=108.8
    ),
    "2022-02-05": dict(
        ap_avg=11, ap3=(12, 15, 7, 9, 15, 12, 12, 6), f107=125.9, f107a=108.5
    ),
}

# PHYSICS.md §4.1, "Actual densities from NRLMSIS 2.1".
DENSITY_TABLE = {
    180.0: (4.118e-10, 4.685e-10),
    200.0: (2.030e-10, 2.360e-10),
    210.0: (1.468e-10, 1.726e-10),
    230.0: (8.010e-11, 9.632e-11),
    260.0: (3.495e-11, 4.347e-11),
}

# Tolerance is set by the table's own precision: the values are quoted to four
# significant figures, so 5e-4 relative is "reproduced exactly as printed".
TABLE_RTOL = 5e-4


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


def test_loader_reproduces_physics_md_table(sw, capsys):
    """PHYSICS.md §6.2 -- the verified Feb 3-5 2022 values, read from the file."""
    with capsys.disabled():
        print("\n  §6.2 loader table:")
    for iso, expected in SW_TABLE.items():
        day = sw.day(datetime.fromisoformat(iso).date())
        assert day.ap_avg == expected["ap_avg"], iso
        assert day.ap3 == tuple(float(v) for v in expected["ap3"]), iso
        assert day.f107_obs == expected["f107"], iso
        assert day.f107_obs_center81 == expected["f107a"], iso
        with capsys.disabled():
            print(
                f"    {iso}  AP_AVG={day.ap_avg:>4.0f}  "
                f"AP1-8={[int(v) for v in day.ap3]}  "
                f"F10.7={day.f107_obs}  CENTER81={day.f107_obs_center81}"
            )


def test_f107_uses_previous_day(sw):
    """PHYSICS.md §6.2 -- f107s takes the previous day's observed value."""
    # 4 Feb's F10.7 input is 3 Feb's observation, 126.5, not 4 Feb's 129.6.
    assert sw.f107(datetime(2022, 2, 4, tzinfo=timezone.utc)) == 126.5
    assert sw.f107(datetime(2022, 2, 5, tzinfo=timezone.utc)) == 129.6
    # f107a is the centred 81-day average of the day itself.
    assert sw.f107a(datetime(2022, 2, 4, tzinfo=timezone.utc)) == 108.8


def test_f107a_is_centred_not_trailing(sw):
    """PHYSICS.md §6.2 -- CENTER81, not LAST81.

    On 2022-02-04 these differ by 7.6 sfu (108.8 vs 101.2), so picking the
    wrong column is a real error, not a rounding difference.
    """
    assert sw.f107a(datetime(2022, 2, 4, tzinfo=timezone.utc)) == 108.8
    assert sw.f107a(datetime(2022, 2, 4, tzinfo=timezone.utc)) != 101.2


def test_ap_array_slot_mapping(sw):
    """PHYSICS.md §6.2 -- the 7-element ap history window, checked by hand.

    This is the off-by-one guard the spec asks for. At 2022-02-04 09:00 UT the
    current slot is AP4 of 4 Feb (index 3, value 22) and the preceding slots
    walk backwards through 4 Feb and into 3 Feb.
    """
    when = datetime(2022, 2, 4, 9, 0, tzinfo=timezone.utc)
    aps = sw.ap_array(when)

    feb3 = SW_TABLE["2022-02-03"]["ap3"]  # 7, 32, 48, 56, 27, 15, 15, 12
    feb4 = SW_TABLE["2022-02-04"]["ap3"]  # 27, 27, 22, 22, 27, 56, 48, 27

    assert aps[0] == 32                   # daily Ap for 4 Feb
    assert aps[1] == feb4[3] == 22        # 09-12 UT slot, current
    assert aps[2] == feb4[2] == 22        # 06-09 UT, 3 h before
    assert aps[3] == feb4[1] == 27        # 03-06 UT, 6 h before
    assert aps[4] == feb4[0] == 27        # 00-03 UT, 9 h before

    # Index 5: slots 4-11 back == all eight slots of 3 Feb.
    assert aps[5] == pytest.approx(np.mean(feb3))
    # Index 6: slots 12-19 back == all eight slots of 2 Feb.
    feb2 = sw.day(datetime(2022, 2, 2).date()).ap3
    assert aps[6] == pytest.approx(np.mean(feb2))


def test_ap_array_slot_boundaries(sw):
    """Every 3 h UT slot maps to the right AP column."""
    feb4 = SW_TABLE["2022-02-04"]["ap3"]
    for slot, expected in enumerate(feb4):
        when = datetime(2022, 2, 4, slot * 3 + 1, tzinfo=timezone.utc)
        assert sw.ap_array(when)[1] == expected, f"slot {slot}"


def test_density_table_4_1(capsys):
    """PHYSICS.md §4.1 -- "Actual densities from NRLMSIS 2.1".

    Quiet is ap = 5, storm is ap = 56. Both use pymsis's *default* options,
    under which NRLMSIS reads only the daily Ap; see the `atmosphere` module
    docstring. The "ap3 = 80" annotation in the §4.1 caption has no effect on
    these numbers and is not required to reproduce them.
    """
    alts = np.array(sorted(DENSITY_TABLE))
    n = len(alts)

    def run(ap):
        return density_from_indices(
            [REF_DATE] * 1, alts, REF_LAT, REF_LON,
            [REF_F107], [REF_F107A], [[ap] * 7],
        )[0]

    quiet, storm = run(5), run(56)

    with capsys.disabled():
        print("\n  §4.1 density table (computed / published):")
        print("    alt      quiet computed    published     storm computed    published")
    for i, alt in enumerate(alts):
        q_pub, s_pub = DENSITY_TABLE[float(alt)]
        with capsys.disabled():
            print(
                f"    {alt:5.0f} km  {quiet[i]:.4e}      {q_pub:.3e}    "
                f"{storm[i]:.4e}      {s_pub:.3e}"
            )
        assert quiet[i] == pytest.approx(q_pub, rel=TABLE_RTOL), f"quiet {alt} km"
        assert storm[i] == pytest.approx(s_pub, rel=TABLE_RTOL), f"storm {alt} km"


def test_three_hourly_ap_needs_storm_time_option(capsys):
    """The 3-hourly ap array is inert unless `storm_time=True`.

    Documents the trap rather than asserting a published number: with default
    options, changing aps[1:] from 5 to 80 changes the density by nothing at
    all. This is why PHYSICS.md §6.2's history assembly must be paired with the
    storm-time option in Phase 2.
    """
    alt = np.array([210.0])
    common = dict(lat_deg=REF_LAT, lon_deg=REF_LON)

    def run(aps, storm_time):
        return density_from_indices(
            [REF_DATE], alt, REF_LAT, REF_LON,
            [REF_F107], [REF_F107A], [aps], storm_time=storm_time,
        )[0, 0]

    flat = [56] * 7
    spiky = [56, 80, 80, 80, 80, 56, 56]

    assert run(flat, False) == run(spiky, False), "default options should ignore aps[1:]"
    assert run(flat, True) != run(spiky, True), "storm-time options should read aps[1:]"

    with capsys.disabled():
        print(
            f"\n  3-hourly ap sensitivity at 210 km:\n"
            f"    default options : flat {run(flat, False):.4e}  "
            f"spiky {run(spiky, False):.4e}  (identical)\n"
            f"    storm_time=True : flat {run(flat, True):.4e}  "
            f"spiky {run(spiky, True):.4e}  "
            f"({(run(spiky, True) / run(flat, True) - 1) * 100:+.2f}%)"
        )


def test_loader_rejects_noncontiguous_data():
    """A gap in the daily series would silently corrupt every ap history."""
    sw_full = SpaceWeather.load(DATA)
    days = dict(sw_full._days)
    del days[datetime(2000, 1, 1).date()]
    with pytest.raises(ValueError, match="not contiguous"):
        SpaceWeather(days)
