"""Physical constants.

PHYSICS.md §2. These values are fixed by the specification; do not substitute
alternatives (e.g. a mean Earth radius) without recording the change there.
"""

MU = 3.986004418e14   # Earth gravitational parameter, m^3/s^2
R_E = 6378.137e3      # Earth equatorial radius, m
G0 = 9.80665          # standard gravity, m/s^2

# PHYSICS.md §3.3: the "unrecoverable" altitude used by Baruah et al. (2024).
# Kept here so it is impossible to disagree between the integrator and the
# outcome taxonomy.
REENTRY_ALTITUDE = 100e3  # m
