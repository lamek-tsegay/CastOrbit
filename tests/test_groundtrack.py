"""Ground-track geometry for the globe view.

sim/groundtrack.py is explicitly NOT part of the validated physics model (see
its module docstring) -- these tests check that it is self-consistent
geometry, not that it matches any PHYSICS.md target.
"""

from __future__ import annotations

import numpy as np
import pytest

from sim.groundtrack import INCLINATION_DEG, cause_of_loss, ground_track


def test_latitude_stays_within_the_inclination_bound():
    """A circular orbit at inclination i never exceeds latitude i."""
    t_s = np.linspace(0, 5 * 5400.0, 200)  # several orbital periods
    h_km = np.full_like(t_s, 210.0)
    lat, _ = ground_track(t_s, h_km, satellite_id=3, n_satellites=49,
                          deploy_time_s=100.0)
    assert max(lat) <= INCLINATION_DEG + 1e-6
    assert min(lat) >= -INCLINATION_DEG - 1e-6


def test_longitude_stays_wrapped():
    t_s = np.linspace(0, 10 * 86400.0, 500)  # many days, many wraps
    h_km = np.full_like(t_s, 210.0)
    _, lon = ground_track(t_s, h_km, satellite_id=0, n_satellites=49,
                          deploy_time_s=0.0)
    assert all(-180.0 <= v <= 180.0 for v in lon)


def test_satellites_are_visually_separated():
    """Different id/deploy_time must not collapse onto the same ground track."""
    t_s = np.array([0.0, 600.0, 1200.0])
    h_km = np.array([210.0, 210.0, 210.0])
    lat_a, lon_a = ground_track(t_s, h_km, 0, 49, 0.0)
    lat_b, lon_b = ground_track(t_s, h_km, 24, 49, 1000.0)
    assert lat_a != lat_b or lon_a != lon_b


def test_decaying_altitude_still_produces_valid_geometry():
    """A satellite whose radius shrinks over time must not produce NaN/inf."""
    t_s = np.linspace(0, 2 * 86400.0, 300)
    h_km = np.linspace(210.0, 100.0, 300)  # monotonic decay to reentry
    lat, lon = ground_track(t_s, h_km, 10, 49, 300.0)
    assert all(np.isfinite(lat)) and all(np.isfinite(lon))
    assert max(lat) <= INCLINATION_DEG + 1e-6


def test_same_number_of_points_as_input():
    t_s = np.array([0.0, 100.0, 200.0, 300.0])
    h_km = np.array([210.0, 209.0, 208.0, 207.0])
    lat, lon = ground_track(t_s, h_km, 5, 49, 0.0)
    assert len(lat) == len(t_s)
    assert len(lon) == len(t_s)


@pytest.mark.parametrize("outcome", ["REENTERED", "INDETERMINATE", "REACHED_SHELL",
                                     "PROPELLANT_EXHAUSTED"])
def test_cause_of_loss_covers_every_outcome(outcome):
    text = cause_of_loss(outcome, "2022-02-05T00:00:00+00:00", 3.5)
    assert isinstance(text, str) and len(text) > 0


def test_cause_of_loss_mentions_the_actual_cd_times_area():
    text = cause_of_loss("REENTERED", "2022-02-05T00:00:00+00:00", 6.19)
    assert "6.19" in text
