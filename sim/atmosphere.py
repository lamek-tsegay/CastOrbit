"""NRLMSIS 2.1 atmosphere, CelesTrak space weather loader, and density cache.

PHYSICS.md §6. Three separable pieces:

  * `SpaceWeather`  -- parses `data/SW-All.csv` and assembles the F10.7 and
    7-element `ap` inputs that NRLMSIS expects (§6.2).
  * `density_from_indices` / `density` -- thin wrappers over `pymsis` (§6.1).
  * `DensityGrid` -- precomputed (time, altitude) cache with interpolation,
    because a per-step pymsis call is far too slow (§6.3).

TRAP -- the 3-hourly `ap` history is ignored by default.
    NRLMSIS 2.1 uses only the daily Ap (`aps[0]`) unless the model option
    `geomagnetic_activity` is set to -1. All the careful slot bookkeeping in
    §6.2 is therefore inert under pymsis's defaults. Pass `storm_time=True` to
    switch the model into the storm-time formulation that actually reads
    `aps[1:]`. The reference densities in PHYSICS.md §4.1 were computed with
    the *default* options, so `storm_time=False` is what reproduces that table.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from pymsis import msis
from scipy.interpolate import RegularGridInterpolator

# Column indices in SW-All.csv, verified against the header in PHYSICS.md §6.2.
_COL_DATE = 0
_COL_AP1 = 12          # AP1..AP8 occupy 12..19
_COL_AP_AVG = 20
_COL_F107_OBS = 24
_COL_F107_DATA_TYPE = 26
_COL_F107_OBS_CENTER81 = 27

_SLOTS_PER_DAY = 8     # AP1..AP8 map to UT slots 00-03, 03-06, ... 21-00


@dataclass(frozen=True)
class SpaceWeatherDay:
    """One row of SW-All.csv, reduced to the fields NRLMSIS needs."""

    day: date
    ap_avg: float
    ap3: tuple[float, ...]        # AP1..AP8, in UT slot order
    f107_obs: float | None        # None when the row is a prediction with no OBS
    f107_obs_center81: float | None
    f107_data_type: str           # OBS / INT / PRD / PRM -- provenance for the writeup

    @property
    def is_observed(self) -> bool:
        return self.f107_data_type == "OBS"


class SpaceWeather:
    """CelesTrak SW-All.csv reader.

    PHYSICS.md §6.2. Three traps handled here:
      1. `F10.7_OBS_CENTER81` (centred), not `LAST81` (trailing).
      2. `f107s` takes the *previous day's* `F10.7_OBS` (NRLMSIS convention).
      3. The file has CRLF line endings.
    """

    def __init__(self, days: dict[date, SpaceWeatherDay]):
        self._days = days
        self._first = min(days)
        self._last = max(days)
        n_days = (self._last - self._first).days + 1
        if n_days != len(days):
            raise ValueError(
                f"SW-All.csv is not contiguous: {len(days)} rows span {n_days} days"
            )
        # Flat 3-hourly ap series, so history windows are plain slicing and an
        # off-by-one is testable rather than silent.
        flat = np.empty(n_days * _SLOTS_PER_DAY, dtype=float)
        for d, row in days.items():
            i = (d - self._first).days * _SLOTS_PER_DAY
            flat[i : i + _SLOTS_PER_DAY] = row.ap3
        self._ap_series = flat

    @classmethod
    def load(cls, path: str | Path) -> "SpaceWeather":
        """Parse SW-All.csv. `newline=''` lets csv handle the CRLF endings.

        The file ends with long-range *monthly* F10.7 predictions (data type
        PRM) that carry no ap columns at all and are not daily-spaced. Those
        rows are skipped: a row without AP_AVG cannot drive NRLMSIS, and
        keeping them would break both the parse and the contiguity check.
        """
        days: dict[date, SpaceWeatherDay] = {}
        skipped = 0
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            if header[_COL_DATE] != "DATE" or header[_COL_AP_AVG] != "AP_AVG":
                raise ValueError(f"unexpected SW-All.csv header layout: {header[:25]}")
            for row in reader:
                if not row or not row[_COL_DATE]:
                    continue
                if not row[_COL_AP_AVG].strip():
                    skipped += 1
                    continue
                d = date.fromisoformat(row[_COL_DATE].strip())
                days[d] = SpaceWeatherDay(
                    day=d,
                    ap_avg=float(row[_COL_AP_AVG]),
                    ap3=tuple(float(row[_COL_AP1 + i]) for i in range(_SLOTS_PER_DAY)),
                    f107_obs=_opt_float(row[_COL_F107_OBS]),
                    f107_obs_center81=_opt_float(row[_COL_F107_OBS_CENTER81]),
                    f107_data_type=row[_COL_F107_DATA_TYPE].strip(),
                )
        if not days:
            raise ValueError(f"no data rows parsed from {path}")
        sw = cls(days)
        sw.n_skipped_rows = skipped
        return sw

    def day(self, when: datetime | date) -> SpaceWeatherDay:
        d = when.date() if isinstance(when, datetime) else when
        try:
            return self._days[d]
        except KeyError:
            raise KeyError(f"{d} outside SW-All.csv range {self._first}..{self._last}")

    def f107(self, when: datetime | date) -> float:
        """F10.7 for NRLMSIS: the *previous day's* observed value (§6.2)."""
        d = when.date() if isinstance(when, datetime) else when
        value = self.day(d - timedelta(days=1)).f107_obs
        if value is None:
            raise ValueError(f"no observed F10.7 for {d - timedelta(days=1)}")
        return value

    def f107a(self, when: datetime | date) -> float:
        """81-day *centred* average of observed F10.7 (§6.2)."""
        value = self.day(when).f107_obs_center81
        if value is None:
            raise ValueError(f"no centred 81-day F10.7 for {when}")
        return value

    def _slot_index(self, when: datetime) -> int:
        """Global index into the flat 3-hourly ap series."""
        d = when.date()
        if d < self._first or d > self._last:
            raise KeyError(f"{d} outside SW-All.csv range")
        return (d - self._first).days * _SLOTS_PER_DAY + when.hour // 3

    def ap_array(self, when: datetime) -> list[float]:
        """The 7-element `aps` entry NRLMSIS expects, per the §6.2 table.

        Index 0 is the daily Ap; 1-4 are the current and three preceding
        3-hourly ap values; 5 and 6 are means of the eight 3-hourly values
        12-33 h and 36-57 h before. Indices 1-6 are only read by the model when
        `storm_time=True` (see the module docstring).
        """
        i = self._slot_index(when)
        if i < 19:
            raise ValueError(
                f"{when} is too close to the start of the file for a 57 h ap history"
            )
        s = self._ap_series
        return [
            self.day(when).ap_avg,          # 0: daily Ap
            s[i],                           # 1: current 3-hourly ap
            s[i - 1],                       # 2: 3 h before
            s[i - 2],                       # 3: 6 h before
            s[i - 3],                       # 4: 9 h before
            float(s[i - 11 : i - 3].mean()),   # 5: mean, slots 4-11 back (12-33 h)
            float(s[i - 19 : i - 11].mean()),  # 6: mean, slots 12-19 back (36-57 h)
        ]


def _opt_float(text: str) -> float | None:
    text = text.strip()
    return float(text) if text else None


def _naive_utc(dates: list[datetime]) -> list[datetime]:
    """Strip tzinfo after converting to UTC.

    numpy cannot represent timezones in datetime64 and pymsis only warns about
    it. Everything in this project is UTC, so the conversion is explicit here
    rather than left to a silent coercion.
    """
    out = []
    for d in dates:
        if d.tzinfo is not None:
            d = d.astimezone(timezone.utc).replace(tzinfo=None)
        out.append(d)
    return out


def _options(storm_time: bool) -> list[float] | None:
    """NRLMSIS switch array. See the module-level trap note."""
    if not storm_time:
        return None
    return msis.create_options(geomagnetic_activity=-1)


def density_from_indices(
    dates: list[datetime],
    alts_km: np.ndarray | list[float],
    lat_deg: float,
    lon_deg: float,
    f107s: list[float],
    f107as: list[float],
    aps: list[list[float]],
    storm_time: bool = False,
) -> np.ndarray:
    """Total mass density, kg/m^3, from explicit space weather indices.

    PHYSICS.md §6.1. `alts_km` is in kilometres, not metres.

    Returns shape (len(dates), len(alts_km)). One latitude/longitude only --
    the single-point density sampling limitation of PHYSICS.md §10.2.
    """
    alts_km = np.atleast_1d(np.asarray(alts_km, dtype=float))
    out = msis.calculate(
        _naive_utc(dates), [lon_deg], [lat_deg], alts_km, f107s, f107as, aps,
        options=_options(storm_time),
    )
    # pymsis returns (ndates, nlons, nlats, nalts, nspecies); index 0 is total
    # mass density. Squeeze only the single lon/lat axes, never the others.
    rho = np.asarray(out)[..., 0].reshape(len(dates), alts_km.size)
    return rho


def density(
    dates: list[datetime],
    alts_m: np.ndarray | list[float],
    lat_deg: float,
    sw: SpaceWeather,
    lon_deg: float = 0.0,
    storm_time: bool = False,
) -> np.ndarray:
    """Density in kg/m^3 with indices pulled from SW-All.csv.

    PHYSICS.md §6.1-6.2. `alts_m` is in **metres** here (converted internally),
    matching the rest of the simulator's SI convention.
    """
    alts_km = np.atleast_1d(np.asarray(alts_m, dtype=float)) / 1e3
    return density_from_indices(
        dates,
        alts_km,
        lat_deg,
        lon_deg,
        [sw.f107(d) for d in dates],
        [sw.f107a(d) for d in dates],
        [sw.ap_array(d) for d in dates],
        storm_time=storm_time,
    )


class DensityGrid:
    """Precomputed (time, altitude) density cache with log-linear interpolation.

    PHYSICS.md §6.3. Density is queried every RK4 stage; a direct pymsis call
    per stage is far too slow. The grid is built once per run and interpolated.

    Interpolation is linear in log(rho) rather than in rho, because density is
    close to exponential in altitude. That is a deviation from the letter of
    §6.3 ("interpolate") in the direction of lower error, and it is what makes
    the <1% accuracy requirement comfortable rather than marginal. Use
    `max_relative_error()` to verify against direct calls.

    The time axis is split into *blocks*. NRLMSIS inputs are piecewise
    constant in time -- the daily Ap and F10.7 switch at midnight UT, and the
    3-hourly ap switches every 3 h -- so density has genuine step
    discontinuities in time. Linear interpolation across a step has a peak
    error equal to the size of the step no matter how fine the spacing, which
    is why a single uniform 30 min grid misses the §6.3 <1% target by a factor
    of five (worst case ~5.6% near 550 km). Instead one interpolator is built
    per 3 h block, using that block's constant indices for every node in it.
    No interpolation cell ever spans a discontinuity, and the step is
    reproduced exactly rather than smeared.
    """

    _BLOCK_S = 3 * 3600.0   # NRLMSIS index blocks are 3 h of UT

    def __init__(
        self,
        epoch: datetime,
        sw: SpaceWeather,
        lat_deg: float,
        lon_deg: float = 0.0,
        duration_s: float = 5 * 86400.0,
        alt_min_km: float = 100.0,
        alt_max_km: float = 600.0,
        alt_step_km: float = 1.0,
        time_step_s: float = 1800.0,
        storm_time: bool = False,
        ap_override: float | None = None,
    ):
        """
        Args:
            ap_override: if set, replaces every element of the `aps` array with
                this value, leaving F10.7 alone. This is how "quiet" conditions
                are constructed for the §9 sweeps: same epoch, same solar flux,
                geomagnetic activity forced flat. PHYSICS.md §4.1 uses ap = 5
                for quiet and ap = 56 for storm, so those are the natural
                values. Set to None to use the real SW-All.csv history.
        """
        if epoch.tzinfo is None:
            epoch = epoch.replace(tzinfo=timezone.utc)
        self.epoch = epoch
        self.lat_deg = lat_deg
        self.lon_deg = lon_deg
        self.storm_time = storm_time
        self.ap_override = ap_override
        self.duration_s = float(duration_s)

        self.alts_km = np.arange(
            alt_min_km, alt_max_km + 0.5 * alt_step_km, alt_step_km, dtype=float
        )
        self._alt0 = float(self.alts_km[0])
        self._alt_step = float(alt_step_km)
        self._n_alt = int(self.alts_km.size)

        self._starts: list[float] = []
        self._interps: list[RegularGridInterpolator] = []
        self._nodes: list[np.ndarray] = []
        self._logrho: list[np.ndarray] = []
        self.n_nodes = 0

        for t0, t1 in self._blocks(epoch, self.duration_s):
            nodes = np.arange(t0, t1, time_step_s, dtype=float)
            if nodes.size == 0 or nodes[-1] < t1:
                nodes = np.append(nodes, t1)
            if nodes.size < 2:
                nodes = np.array([t0, t1], dtype=float)

            # Indices are constant across the block by construction, so they are
            # read once at the block start and reused for every node. Evaluating
            # the closing node this way yields the block's left-limit density,
            # which is exactly what interpolation up to the boundary needs.
            ref = self._to_datetime(t0)
            f107 = sw.f107(ref)
            f107a = sw.f107a(ref)
            aps = (
                [float(ap_override)] * 7 if ap_override is not None
                else sw.ap_array(ref)
            )

            n = nodes.size
            rho = density_from_indices(
                [self._to_datetime(t) for t in nodes],
                self.alts_km,
                lat_deg,
                lon_deg,
                [f107] * n,
                [f107a] * n,
                [aps] * n,
                storm_time=storm_time,
            )
            if not np.all(np.isfinite(rho)) or np.any(rho <= 0.0):
                raise ValueError("pymsis returned non-positive or non-finite density")

            log_rho = np.log(rho)
            self._starts.append(t0)
            self._nodes.append(nodes)
            self._logrho.append(log_rho)
            self._interps.append(
                RegularGridInterpolator(
                    (nodes, self.alts_km),
                    log_rho,
                    method="linear",
                    bounds_error=False,
                    fill_value=None,  # linear extrapolation, never NaN
                )
            )
            self.n_nodes += n

        self._starts_arr = np.array(self._starts, dtype=float)

    @classmethod
    def _blocks(cls, epoch: datetime, duration_s: float) -> list[tuple[float, float]]:
        """[start, end) offsets in seconds for each constant-index interval.

        The first block is partial when the epoch is not on a 3 h UT boundary
        (the Feb 2022 launch epoch, 18:13 UT, is not).
        """
        secs_into_day = (
            epoch.hour * 3600 + epoch.minute * 60 + epoch.second
            + epoch.microsecond * 1e-6
        )
        first_boundary = (-secs_into_day) % cls._BLOCK_S
        edges = [0.0]
        b = first_boundary
        if b == 0.0:
            b = cls._BLOCK_S
        while b < duration_s:
            edges.append(b)
            b += cls._BLOCK_S
        edges.append(max(duration_s, edges[-1] + cls._BLOCK_S))
        return list(zip(edges[:-1], edges[1:]))

    def _to_datetime(self, t_s: float) -> datetime:
        return self.epoch + timedelta(seconds=float(t_s))

    def _block_index(self, t_s: float) -> int:
        k = int(np.searchsorted(self._starts_arr, t_s, side="right")) - 1
        return min(max(k, 0), len(self._interps) - 1)

    def __call__(self, t_s: float, h_m: float) -> float:
        """Density at `t_s` seconds after epoch and altitude `h_m` metres.

        Altitude is clamped to the tabulated range. Below the 100 km floor the
        grid has no data, and log-linear extrapolation there runs away badly
        enough to drive the state negative inside a single RK4 step. 100 km is
        the termination altitude (PHYSICS.md §3.3), so the clamp only ever
        affects the interior stages of the one step that crosses it, and the
        reported reentry time is unaffected to well under a step.
        """
        h_km = min(max(h_m / 1e3, self.alts_km[0]), self.alts_km[-1])
        interp = self._interps[self._block_index(t_s)]
        return float(np.exp(interp((t_s, h_km))))

    def max_relative_error(self, sw: SpaceWeather, n_alt: int = 40) -> float:
        """Largest relative interpolation error against direct pymsis calls.

        Probes cell midpoints in both axes -- where linear interpolation is
        worst -- across every block. PHYSICS.md §6.3 requires this to be under
        1%. Probes use `sw` directly, so this checks the block bookkeeping as
        well as the interpolation.
        """
        h_probe_km = np.linspace(self.alts_km[0] + 0.5, self.alts_km[-1] - 0.5, n_alt)
        worst = 0.0
        for k, interp in enumerate(self._interps):
            nodes = interp.grid[0]
            t_probe = (nodes[:-1] + nodes[1:]) / 2.0
            if t_probe.size == 0:
                continue
            dates = [self._to_datetime(t) for t in t_probe]
            exact = density_from_indices(
                dates, h_probe_km, self.lat_deg, self.lon_deg,
                [sw.f107(d) for d in dates],
                [sw.f107a(d) for d in dates],
                [sw.ap_array(d) for d in dates],
                storm_time=self.storm_time,
            )
            tt, hh = np.meshgrid(t_probe, h_probe_km, indexing="ij")
            approx = np.exp(interp(np.stack([tt, hh], axis=-1)))
            worst = max(worst, float(np.max(np.abs(approx - exact) / exact)))
        return worst
