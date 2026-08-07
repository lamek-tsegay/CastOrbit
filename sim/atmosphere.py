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


