"""Mission definition and the lifecycle run. V2_BRIEF.md §7, Phase 7.

V1 could answer one question about one event. This module is the generalised
entry point: an arbitrary altitude, inclination, epoch and target shell, flown
on the same engine, ending in a compliance verdict.

Nothing here is new physics. It selects an atmosphere model, calls
`propagate_adaptive`, and hands the result to `disposal.assess_compliance`.
The equations are still `PHYSICS.md` §3.2 and the disposal maths documented in
`docs/PHYSICS_V2.md`.

**The one judgement call is which atmosphere to use**, and it is made
explicitly rather than by default. A five-day run in February 2022 should use
the real 3-hourly indices; a 25-year disposal run cannot, because the space
weather file ends a few weeks out. `choose_atmosphere` decides, records why,
and the reason travels in the result -- a run that quietly swapped its
atmosphere model would be the single most misleading thing this module could
do.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone

import numpy as np

from .atmosphere import (
    SOLAR_ACTIVITY_PERCENTILES,
    ClimatologyDensity,
    DensityGrid,
    SpaceWeather,
)
from .constants import R_E
from .disposal import (
    SECONDS_PER_YEAR,
    ComplianceResult,
    DisposalRule,
    assess_compliance,
)
from .dynamics import derivatives
from .integrator import propagate_adaptive
from .satellite import Outcome, SatelliteConfig, ThrusterMode

# Beyond this the per-3h-block DensityGrid stops being affordable, regardless
# of whether the file happens to cover the window. Sixty days is ~480 blocks.
MAX_REAL_WEATHER_DAYS = 60.0

# DensityGrid needs 57 h of prior ap history for the §6.2 slot arithmetic.
_AP_HISTORY_MARGIN = timedelta(hours=60)


@dataclass
class Mission:
    """What to fly. Arbitrary altitude, inclination, epoch and target.

    `inclination_deg` doubles as the density sampling latitude, which is the
    convention V1 used for the 53.22 deg Starlink case. A satellite at
    inclination i spends most of its time between +/-i, and single-point
    sampling is the PHYSICS.md §10.2 limitation -- carried forward here
    unchanged rather than quietly improved, so V2 results stay comparable to
    the validated V1 ones.
    """

    epoch: datetime
    insertion_altitude_m: float
    inclination_deg: float
    satellite: SatelliteConfig
    mission_duration_years: float = 5.0
    target_altitude_m: float | None = None       # operational shell, if raising
    thruster_mode: ThrusterMode = ThrusterMode.NOMINAL
    label: str = "unnamed"

    def __post_init__(self):
        if self.epoch.tzinfo is None:
            self.epoch = self.epoch.replace(tzinfo=timezone.utc)
        if self.insertion_altitude_m <= 0:
            raise ValueError("insertion altitude must be positive")
        if not -90.0 <= self.inclination_deg <= 90.0:
            # NRLMSIS takes a geodetic latitude. An inclination above 90 deg is
            # a retrograde orbit whose sampling latitude is 180 - i.
            raise ValueError(
                f"inclination {self.inclination_deg} deg is outside the "
                "latitude range NRLMSIS is sampled at; pass 180 - i for a "
                "retrograde orbit"
            )
        if self.target_altitude_m is not None:
            if self.target_altitude_m <= self.insertion_altitude_m:
                raise ValueError(
                    "target shell must be above the insertion altitude; "
                    "lowering is a disposal manoeuvre, not a mission"
                )

    @property
    def thrust_n(self) -> float:
        """F = 0 in safe mode, rated thrust otherwise. PHYSICS.md §5."""
        if self.thruster_mode is ThrusterMode.SAFE_MODE:
            return 0.0
        return self.satellite.thrust_n

    @property
    def duration_s(self) -> float:
        return self.mission_duration_years * SECONDS_PER_YEAR


@dataclass
class AtmosphereChoice:
    """Which density model was used, and why. Never inferred by the caller."""

    model: str                  # "DensityGrid" | "ClimatologyDensity"
    reason: str
    solar_activity: str | None  # None when real indices were used
    uses_real_space_weather: bool

    def as_dict(self) -> dict:
        return {
            "model": self.model,
            "reason": self.reason,
            "solar_activity": self.solar_activity,
            "uses_real_space_weather": self.uses_real_space_weather,
        }


def choose_atmosphere(
    mission: Mission,
    sw: SpaceWeather,
    solar_activity: str = "mean",
    storm_time: bool = False,
    alt_max_km: float = 2000.0,
) -> tuple[object, AtmosphereChoice]:
    """Pick the density model this mission can actually be flown on.

    Real 3-hourly indices win whenever they exist and are affordable, because
    they are what V1 was validated against. The climatology takes over when
    either condition fails, and says which one.
    """
    window_days = mission.duration_s / 86400.0
    end = mission.epoch + timedelta(seconds=mission.duration_s)
    last_covered = datetime(
        sw.last_day.year, sw.last_day.month, sw.last_day.day, tzinfo=timezone.utc
    )
    first_covered = mission.epoch - _AP_HISTORY_MARGIN

    covered = end <= last_covered and first_covered.date() >= sw._first
    affordable = window_days <= MAX_REAL_WEATHER_DAYS

    if covered and affordable:
        grid = DensityGrid(
            mission.epoch,
            sw,
            lat_deg=mission.inclination_deg,
            duration_s=mission.duration_s,
            storm_time=storm_time,
            alt_max_km=min(alt_max_km, 1000.0),
        )
        return grid, AtmosphereChoice(
            model="DensityGrid",
            reason=(
                f"{window_days:.1f}-day window lies inside SW-All.csv "
                f"(through {sw.last_day}), so real 3-hourly indices are used."
            ),
            solar_activity=None,
            uses_real_space_weather=True,
        )

    if not covered:
        reason = (
            f"Window ends {end.date()}, past the end of SW-All.csv "
            f"({sw.last_day}). Future space weather is not predictable, so "
            f"indices are held at the '{solar_activity}' activity level "
            f"(p{SOLAR_ACTIVITY_PERCENTILES[solar_activity]:.0f} of the "
            "observed record) and the answer should be read as a band across "
            "levels, not a single value."
        )
    else:
        reason = (
            f"{window_days:.0f}-day window exceeds the "
            f"{MAX_REAL_WEATHER_DAYS:.0f}-day limit for per-3h-block density "
            f"grids. Indices held at the '{solar_activity}' level."
        )

    clim = ClimatologyDensity.for_level(
        solar_activity, sw, lat_deg=mission.inclination_deg, alt_max_km=alt_max_km
    )
    return clim, AtmosphereChoice(
        model="ClimatologyDensity",
        reason=reason,
        solar_activity=solar_activity,
        uses_real_space_weather=False,
    )


@dataclass
class LifecycleResult:
    """The full V2 answer: what the orbit did, then whether it complies."""

    mission_label: str
    outcome: Outcome
    outcome_time_s: float | None
    final_altitude_km: float
    atmosphere: AtmosphereChoice
    step_stats: dict
    compliance: ComplianceResult | None = None
    t_s: np.ndarray = field(default=None, repr=False)
    h_km: np.ndarray = field(default=None, repr=False)

    def as_dict(self) -> dict:
        return {
            "mission": self.mission_label,
            "outcome": self.outcome.value,
            "outcome_time_s": self.outcome_time_s,
            "final_altitude_km": self.final_altitude_km,
            "atmosphere": self.atmosphere.as_dict(),
            "step_stats": self.step_stats,
            "compliance": None if self.compliance is None
            else self.compliance.as_dict(),
        }


def fly(
    mission: Mission,
    sw: SpaceWeather,
    solar_activity: str = "mean",
    storm_time: bool = False,
    tol: float = 1e-4,
    dt_max: float | None = None,
    assess_disposal: bool = True,
    rule: DisposalRule | None = None,
) -> LifecycleResult:
    """Propagate the mission, then price its disposal.

    `dt_max` defaults to 1800 s on real space weather -- the density grid's own
    node spacing, so the integrator cannot step past the atmosphere it was
    given -- and to 30 days on the climatology, where the profile is constant
    in time and only the altitude change needs resolving.
    """
    atmosphere, choice = choose_atmosphere(
        mission, sw, solar_activity=solar_activity, storm_time=storm_time
    )
    if dt_max is None:
        dt_max = 1800.0 if choice.uses_real_space_weather else 30 * 86400.0

    sat = mission.satellite
    thrust = mission.thrust_n

    def deriv(t: float, y: np.ndarray) -> np.ndarray:
        return derivatives(
            y,
            atmosphere(t, y[0] - R_E),
            thrust=thrust,
            cd=sat.cd,
            area=sat.area_m2,
            isp=sat.isp_s,
        )

    traj = propagate_adaptive(
        deriv,
        np.array([R_E + mission.insertion_altitude_m, sat.mass_kg]),
        t_max=mission.duration_s,
        tol=tol,
        dt_max=dt_max,
        shell_altitude_m=mission.target_altitude_m,
        dry_mass_kg=sat.dry_mass_kg,
    )

    compliance = None
    if assess_disposal:
        # Disposal starts from wherever the satellite actually ended up, at
        # whatever mass it has left -- not from the nominal shell.
        h_eol_m = float(traj.a_m[-1] - R_E)
        mass_eol = float(traj.mass_kg[-1])
        available = (
            0.0 if sat.dry_mass_kg is None
            else max(0.0, mass_eol - sat.dry_mass_kg)
        )
        if traj.outcome is not Outcome.REENTERED:
            compliance = assess_compliance(
                atmosphere,
                h_eol_m,
                mass_at_eol_kg=mass_eol,
                cd=sat.cd,
                area_m2=sat.area_m2,
                isp_s=sat.isp_s or 0.0,
                propellant_available_kg=available,
                solar_activity=choice.solar_activity or "observed",
                rule=rule,
            ) if sat.isp_s else None

    return LifecycleResult(
        mission_label=mission.label,
        outcome=traj.outcome,
        outcome_time_s=traj.outcome_time_s,
        final_altitude_km=float(traj.h_km[-1]),
        atmosphere=choice,
        step_stats=traj.stats.as_dict(),
        compliance=compliance,
        t_s=traj.t_s,
        h_km=traj.h_km,
    )


SOLAR_BAND_CAVEAT = (
    "Bounds, not scenarios. Each level holds solar activity fixed at a "
    "percentile of the 1957-onward observed record for the whole run. No real "
    "25-year period sits at one level -- the solar cycle is ~11 years, so a "
    "sustained solar maximum is physically impossible, not merely unlikely. "
    "Read low/high as bounds the true answer lies between, and never as "
    "best-case and worst-case futures."
)


@dataclass
class SolarBandResult:
    """Results at every solar activity level, with the caveat attached.

    This exists as a type rather than a plain dict so that
    `SOLAR_BAND_CAVEAT` cannot be separated from the numbers on the way to a
    UI. Anything rendering the band gets `caveat` in the same payload; there
    is no accessor that returns the spread without it.
    """

    levels: dict[str, LifecycleResult]
    caveat: str = SOLAR_BAND_CAVEAT

    def __getitem__(self, level: str) -> LifecycleResult:
        return self.levels[level]

    def __iter__(self):
        return iter(self.levels)

    def items(self):
        return self.levels.items()

    def as_dict(self) -> dict:
        return {
            "caveat": self.caveat,
            "levels": {k: v.as_dict() for k, v in self.levels.items()},
        }


def fly_solar_band(
    mission: Mission,
    sw: SpaceWeather,
    **kwargs,
) -> SolarBandResult:
    """Fly the same mission at every solar activity level.

    "Uncertain values get swept, not chosen" (V2_BRIEF.md §6). Future solar
    activity is the largest uncertainty in any multi-year run, so the intended
    way to use `fly` for a long mission is through this, and to read the
    spread as the answer.

    Returns a `SolarBandResult`, which carries `SOLAR_BAND_CAVEAT` alongside
    the numbers -- see that class for why it is a type and not a dict.
    """
    return SolarBandResult(
        levels={
            level: fly(mission, sw, solar_activity=level, **kwargs)
            for level in SOLAR_ACTIVITY_PERCENTILES
        }
    )
