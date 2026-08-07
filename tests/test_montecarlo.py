"""Monte Carlo sampling and the ensemble propagator. PHYSICS.md §9.

The load-bearing test here is `test_ensemble_matches_scalar_propagation`: the
vectorised batch integrator must reproduce the scalar `propagate` path used for
the Baruah validation. If those two ever diverge, every Phase 3 number is
suspect while the Phase 2 numbers still look fine.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pytest

from sim.atmosphere import SpaceWeather
from sim.constants import R_E
from sim.dynamics import da_dt, derivatives
from sim.montecarlo import (
    CD_RANGE,
    MASS_MEAN_KG,
    N_SATELLITES,
    RAM_AREA_RANGE,
    Batch,
    build_grid,
    run_batch,
    sample_batch,
)
from sim.satellite import Outcome

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"
EPOCH = datetime(2022, 2, 3, 18, 13, tzinfo=timezone.utc) - timedelta(seconds=1800)


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


@pytest.fixture(scope="module")
def storm_grid(sw):
    return build_grid(sw, EPOCH, storm=True, duration_s=6 * 86400.0)


def test_sample_batch_respects_section_9_distributions():
    """PHYSICS.md §9 distribution table."""
    rng = np.random.default_rng(0)
    b = sample_batch(rng, 210.0, n=20000)

    assert len(b) == 20000
    assert RAM_AREA_RANGE[0] <= b.area_m2.min()
    assert b.area_m2.max() <= RAM_AREA_RANGE[1]
    assert CD_RANGE[0] <= b.cd.min() and b.cd.max() <= CD_RANGE[1]
    assert 0.0 <= b.deploy_time_s.min() and b.deploy_time_s.max() <= 3600.0

    assert b.mass_kg.mean() == pytest.approx(MASS_MEAN_KG, rel=1e-2)
    assert b.mass_kg.std() == pytest.approx(0.03 * MASS_MEAN_KG, rel=5e-2)
    assert b.insertion_altitude_m.mean() == pytest.approx(210e3, rel=1e-3)
    assert b.insertion_altitude_m.std() == pytest.approx(2e3, rel=5e-2)
    assert b.thrust_n.mean() == pytest.approx(0.071, rel=1e-2)


def test_default_batch_is_49_satellites():
    rng = np.random.default_rng(0)
    assert len(sample_batch(rng, 210.0)) == N_SATELLITES == 49


def test_dynamics_are_array_safe():
    """The ensemble and the scalar runs must share one equation, not two."""
    a = np.full(4, R_E + 210e3)
    m = np.full(4, 227.0)
    rho = np.full(4, 1.5e-10)
    vec = da_dt(a, m, rho, np.full(4, 0.071), np.full(4, 2.2), np.full(4, 2.0))
    scalar = da_dt(R_E + 210e3, 227.0, 1.5e-10, 0.071, 2.2, 2.0)
    assert vec.shape == (4,)
    assert np.allclose(vec, scalar, rtol=0, atol=0)


def test_ensemble_matches_scalar_propagation(storm_grid, capsys):
    """The vectorised batch loop must agree with the scalar `propagate`.

    One satellite, identical parameters, no deployment offset, safe mode off.
    Any divergence means the masking or the staging in `run_batch` is wrong.
    """
    from sim.integrator import propagate
    from sim.montecarlo import ISP_S

    area, mass, cd, thrust = 3.0, 227.0, 2.2, 0.071
    h0, dt, t_max = 210e3, 60.0, 2 * 86400.0

    batch = Batch(
        area_m2=np.array([area]),
        mass_kg=np.array([mass]),
        insertion_altitude_m=np.array([h0]),
        deploy_time_s=np.array([0.0]),
        thrust_n=np.array([thrust]),
        cd=np.array([cd]),
    )
    ens = run_batch(batch, storm_grid, dt=dt, t_max_s=t_max, sample_every=10**9)

    def deriv(t, y):
        rho = storm_grid(t, y[0] - R_E)
        return derivatives(y, rho, thrust, cd, area, isp=ISP_S)

    scalar = propagate(
        deriv, np.array([R_E + h0, mass]), dt=dt, t_max=t_max,
        shell_altitude_m=550e3,
    )

    rel = abs(ens.final_altitude_m[0] - (scalar.a_m[-1] - R_E)) / (
        scalar.a_m[-1] - R_E
    )
    with capsys.disabled():
        print(f"\n  ensemble vs scalar after 2 days: "
              f"{ens.final_altitude_m[0] / 1e3:.6f} km vs "
              f"{(scalar.a_m[-1] - R_E) / 1e3:.6f} km, rel diff {rel:.2e}")
    assert rel < 1e-12


def test_undeployed_satellites_do_not_move(storm_grid):
    """A satellite must hold its insertion altitude until its deployment time."""
    batch = Batch(
        area_m2=np.array([4.0, 4.0]),
        mass_kg=np.array([227.0, 227.0]),
        insertion_altitude_m=np.array([210e3, 210e3]),
        deploy_time_s=np.array([0.0, 3600.0]),
        thrust_n=np.array([0.0, 0.0]),
        cd=np.array([2.2, 2.2]),
    )
    r = run_batch(batch, storm_grid, dt=60.0, t_max_s=3600.0, sample_every=10**9)
    # The second satellite is deployed exactly at the end, so it must be intact.
    assert r.final_altitude_m[1] == pytest.approx(210e3, abs=1e-6)
    assert r.final_altitude_m[0] < 210e3 - 100.0   # the first has been decaying


def test_safe_mode_decays_and_nominal_climbs(storm_grid):
    """PHYSICS.md §5: with F = 0 the fleet sinks; with thrusters on it climbs."""
    rng = np.random.default_rng(1)
    batch = sample_batch(rng, 210.0)

    safe = run_batch(batch, storm_grid, dt=60.0, t_max_s=2 * 86400.0,
                     safe_mode_exit_s=np.inf, sample_every=10**9)
    nominal = run_batch(batch, storm_grid, dt=60.0, t_max_s=2 * 86400.0,
                        safe_mode_exit_s=None, sample_every=10**9)

    assert np.all(safe.final_altitude_m < 210e3)
    assert np.mean(nominal.final_altitude_m) > 210e3


def test_propellant_exhausted_is_never_reported(storm_grid):
    """Documented limitation: no propellant mass is published, so it is not modelled."""
    rng = np.random.default_rng(2)
    r = run_batch(sample_batch(rng, 210.0), storm_grid, dt=120.0,
                  t_max_s=86400.0, sample_every=10**9)
    assert r.counts()[Outcome.PROPELLANT_EXHAUSTED.value] == 0


def test_quiet_grid_is_thinner_than_storm(sw):
    """The quiet counterfactual must actually be less dense at 210 km."""
    storm = build_grid(sw, EPOCH, storm=True, duration_s=2 * 86400.0)
    quiet = build_grid(sw, EPOCH, storm=False, duration_s=2 * 86400.0)
    for t in (0.0, 43200.0, 86400.0):
        assert quiet(t, 210e3) < storm(t, 210e3)


def test_batch_is_reproducible_from_seed(storm_grid):
    a = run_batch(sample_batch(np.random.default_rng(7), 195.0), storm_grid,
                  dt=120.0, t_max_s=3 * 86400.0, sample_every=10**9)
    b = run_batch(sample_batch(np.random.default_rng(7), 195.0), storm_grid,
                  dt=120.0, t_max_s=3 * 86400.0, sample_every=10**9)
    assert a.survival_fraction == b.survival_fraction
    assert np.array_equal(a.final_altitude_m, b.final_altitude_m)
