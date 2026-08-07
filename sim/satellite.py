"""Satellite configuration and run outcome taxonomy.

PHYSICS.md §5 (thruster mode) and §3.3 (termination conditions).

Phase 1 needs only the enums and a plain parameter container; the Monte Carlo
sampling that populates them lives in `montecarlo.py` (Phase 3).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ThrusterMode(Enum):
    """PHYSICS.md §5.

    Thruster state is an explicit, time-dependent input, never a constant.
    The February 2022 loss is only reproducible because the fleet was commanded
    into SAFE_MODE: F = 0 with a knife-edge ram area, so drag acted alone.
    """

    NOMINAL = "nominal"      # F = rated thrust, A = nominal ram area
    SAFE_MODE = "safe_mode"  # F = 0,            A = knife-edge ram area


class Outcome(Enum):
    """PHYSICS.md §3.3 and §9."""

    REACHED_SHELL = "REACHED_SHELL"
    REENTERED = "REENTERED"
    PROPELLANT_EXHAUSTED = "PROPELLANT_EXHAUSTED"
    INDETERMINATE = "INDETERMINATE"


@dataclass
class SatelliteConfig:
    """Per-satellite physical parameters.

    `area_m2` is the ram area for the *current* thruster mode; a run that
    switches mode switches both `thrust_n` and `area_m2` together (PHYSICS.md §5).
    """

    mass_kg: float
    area_m2: float
    cd: float
    thrust_n: float = 0.0
    isp_s: float | None = None       # None => no mass depletion (dm/dt = 0)
    dry_mass_kg: float | None = None  # propellant exhausted when mass hits this
