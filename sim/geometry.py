"""Parametric geometry and the projected-area solver. V2_BRIEF.md §7, Phase 8.

**Why this module is the point of V2.** The README's central finding is that
`rho`, `Cd` and `A` are inseparable from a decay curve -- only the product is
observable. V1 could only respond by sweeping `A` across a disputed range. If
`A` is instead *derived* from geometry the user specified, one of the three is
constrained by construction rather than assumed, and the degeneracy is
attacked rather than documented (V2_BRIEF.md §2).

So the job here is narrow and specific: given a shape and an attitude, return
the area presented to the oncoming flow.

---

**The physics this is valid for.** At 200-700 km the flow is free-molecular:
the mean free path is kilometres, so molecules hit the satellite and are
re-emitted without colliding with each other on the way in. There is no
boundary layer and no wake pressure recovery. Drag is therefore proportional
to the area geometrically projected onto the plane perpendicular to the
velocity vector, which is exactly what this module computes, and the
accommodation physics is bundled into `Cd` (PHYSICS.md §2).

**What is deliberately not modelled:**

  * *Aerodynamic torque and attitude dynamics.* Attitude is an input, per
    PHYSICS.md §1's exclusion. This module says what area an attitude
    presents, never what attitude the satellite will settle into.
  * *Concave self-shadowing at grazing angles*, beyond what the union solver
    below catches -- the union is exact for the projection, but a rear surface
    sitting in the aerodynamic shadow of a front one still receives some
    re-emitted flux in reality.
  * *Panel flexure.* Solar arrays are treated as rigid flat plates.

---

**Two independent solvers, deliberately.** `projected_area` sums each part's
analytic projection; `projected_area_union` projects every vertex, takes each
part's convex hull in the projection plane, and rasterises the union. They
must agree exactly when no two parts overlap in projection, and the union is
the correct one when they do. Their difference *is* the self-shadowing
correction, which is a quantity worth having rather than an error to
eliminate. This mirrors the "two independent code paths agree" check the
project already relies on between `montecarlo.py` and `critical.py`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

Vec3 = tuple[float, float, float]


def unit(v: Sequence[float]) -> np.ndarray:
    """Normalise a vector, rejecting the zero vector rather than returning NaN."""
    arr = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(arr))
    if n == 0.0 or not math.isfinite(n):
        raise ValueError(f"cannot normalise {v!r}")
    return arr / n


def direction_from_angles(azimuth_deg: float, elevation_deg: float) -> np.ndarray:
    """Unit vector from spherical angles, body frame.

    Azimuth is measured in the x-y plane from +x; elevation from that plane
    towards +z. Provided so attitudes can be written readably in specs and
    tests instead of as raw triples.
    """
    az, el = math.radians(azimuth_deg), math.radians(elevation_deg)
    return np.array([
        math.cos(el) * math.cos(az),
        math.cos(el) * math.sin(az),
        math.sin(el),
    ])


# --------------------------------------------------------------------------
# Parts
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Box:
    """Rectangular bus, axis-aligned in the body frame.

    `length_m` runs along +x, `width_m` along +y, `thickness_m` along +z.
    """

    length_m: float
    width_m: float
    thickness_m: float
    center_m: Vec3 = (0.0, 0.0, 0.0)
    name: str = "bus"

    def __post_init__(self):
        for field_name in ("length_m", "width_m", "thickness_m"):
            if getattr(self, field_name) <= 0.0:
                raise ValueError(f"{field_name} must be positive")

    def projected_area(self, v_hat: np.ndarray) -> float:
        """Area presented to `v_hat`, analytically.

        For a convex body the projected area is `(1/2) * sum_i A_i |n_i . v|`
        over all faces. For a box the six faces pair up into three terms:

            A = |vx|*(W*T) + |vy|*(L*T) + |vz|*(L*W)
        """
        return float(
            abs(v_hat[0]) * self.width_m * self.thickness_m
            + abs(v_hat[1]) * self.length_m * self.thickness_m
            + abs(v_hat[2]) * self.length_m * self.width_m
        )

    def vertices(self) -> np.ndarray:
        cx, cy, cz = self.center_m
        hx, hy, hz = self.length_m / 2, self.width_m / 2, self.thickness_m / 2
        return np.array([
            [cx + sx * hx, cy + sy * hy, cz + sz * hz]
            for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)
        ])


@dataclass(frozen=True)
class Panel:
    """Flat rectangular plate -- a solar array wing or a radiator.

    `normal` is the plate's face normal in the body frame. `span_m` and
    `chord_m` are its in-plane dimensions. A real array has a small but
    non-zero thickness, so `thickness_m` defaults to 1 cm rather than 0: a
    perfectly zero-thickness plate presents exactly zero area edge-on, which
    would make a feathered array look drag-free and flatter the knife-edge
    result.
    """

    span_m: float
    chord_m: float
    normal: Vec3 = (0.0, 0.0, 1.0)
    center_m: Vec3 = (0.0, 0.0, 0.0)
    thickness_m: float = 0.01
    name: str = "panel"

    def __post_init__(self):
        for field_name in ("span_m", "chord_m"):
            if getattr(self, field_name) <= 0.0:
                raise ValueError(f"{field_name} must be positive")
        if self.thickness_m < 0.0:
            raise ValueError("thickness_m must be non-negative")

    @property
    def area_m2(self) -> float:
        return self.span_m * self.chord_m

    def _frame(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Orthonormal (span, chord, normal) triad for this plate."""
        n = unit(self.normal)
        seed = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(seed, n))) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        e_span = unit(np.cross(seed, n))
        e_chord = np.cross(n, e_span)
        return e_span, e_chord, n

    def projected_area(self, v_hat: np.ndarray) -> float:
        """`A*|cos(theta)|` for the face, plus the two edge contributions.

        Treated as a thin box, so it reduces to the flat-plate cosine law when
        `thickness_m` is small -- which is the Phase 8 gate's analytic check.
        """
        e_span, e_chord, n = self._frame()
        return float(
            abs(float(np.dot(v_hat, n))) * self.span_m * self.chord_m
            + abs(float(np.dot(v_hat, e_span))) * self.chord_m * self.thickness_m
            + abs(float(np.dot(v_hat, e_chord))) * self.span_m * self.thickness_m
        )

    def vertices(self) -> np.ndarray:
        e_span, e_chord, n = self._frame()
        c = np.asarray(self.center_m, dtype=float)
        hs, hc, ht = self.span_m / 2, self.chord_m / 2, self.thickness_m / 2
        return np.array([
            c + ss * hs * e_span + sc * hc * e_chord + sn * ht * n
            for ss in (-1, 1) for sc in (-1, 1) for sn in (-1, 1)
        ])


Part = Box | Panel


# --------------------------------------------------------------------------
# Projection plane geometry
# --------------------------------------------------------------------------

def _projection_basis(v_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two orthonormal vectors spanning the plane perpendicular to `v_hat`."""
    seed = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(seed, v_hat))) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    e1 = unit(np.cross(seed, v_hat))
    e2 = np.cross(v_hat, e1)
    return e1, e2


def _convex_hull_2d(points: np.ndarray) -> np.ndarray:
    """Andrew's monotone chain. Returns hull vertices counter-clockwise.

    Written out rather than pulled from scipy: ARCHITECTURE.md §4 restricts
    scipy to interpolation, and a 2D hull is twenty lines.
    """
    pts = np.unique(np.round(points, 12), axis=0)
    pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]
    if pts.shape[0] <= 2:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[np.ndarray] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: list[np.ndarray] = []
    for p in pts[::-1]:
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def _polygon_area(hull: np.ndarray) -> float:
    """Shoelace formula."""
    if hull.shape[0] < 3:
        return 0.0
    x, y = hull[:, 0], hull[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _points_in_hull(hull: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Boolean mask of which `pts` lie inside a counter-clockwise convex hull."""
    if hull.shape[0] < 3:
        return np.zeros(pts.shape[0], dtype=bool)
    inside = np.ones(pts.shape[0], dtype=bool)
    for i in range(hull.shape[0]):
        a, b = hull[i], hull[(i + 1) % hull.shape[0]]
        edge = b - a
        rel = pts - a
        inside &= (edge[0] * rel[:, 1] - edge[1] * rel[:, 0]) >= -1e-12
    return inside


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

@dataclass
class Geometry:
    """An assembly of parts in a shared body frame."""

    parts: list[Part]
    label: str = "unnamed"

    def __post_init__(self):
        if not self.parts:
            raise ValueError("geometry needs at least one part")

    def projected_area(self, v_hat: Sequence[float]) -> float:
        """Sum of the parts' analytic projections. Ignores mutual shadowing.

        This is an upper bound on the true projected area: two parts that
        overlap in projection are counted twice. Use `projected_area_union`
        for the shadowed value, and `shadowing_fraction` for the difference.
        """
        v = unit(v_hat)
        return float(sum(p.projected_area(v) for p in self.parts))

    def hulls(self, v_hat: Sequence[float]) -> list[np.ndarray]:
        """Each part's outline in the projection plane perpendicular to `v_hat`."""
        v = unit(v_hat)
        e1, e2 = _projection_basis(v)
        out = []
        for part in self.parts:
            w = part.vertices()
            out.append(_convex_hull_2d(np.column_stack([w @ e1, w @ e2])))
        return out

    def projected_area_union(
        self, v_hat: Sequence[float], grid: int = 600
    ) -> float:
        """Union of the projected outlines, by rasterising the projection plane.

        Independent of `projected_area`: it goes through vertices, hulls and a
        pixel count rather than through per-part area formulas, so agreement
        between the two is meaningful rather than circular.

        Single parts short-circuit to the exact hull area, so the flat-plate
        and single-box cases carry no discretisation error at all. For
        multi-part assemblies the rasterisation error scales as ~1/grid.
        """
        hulls = [h for h in self.hulls(v_hat) if h.shape[0] >= 3]
        if not hulls:
            return 0.0
        if len(hulls) == 1:
            return _polygon_area(hulls[0])

        allpts = np.vstack(hulls)
        lo, hi = allpts.min(axis=0), allpts.max(axis=0)
        pad = 0.02 * np.maximum(hi - lo, 1e-9)
        lo, hi = lo - pad, hi + pad

        xs = np.linspace(lo[0], hi[0], grid)
        ys = np.linspace(lo[1], hi[1], grid)
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        pts = np.column_stack([gx.ravel(), gy.ravel()])

        covered = np.zeros(pts.shape[0], dtype=bool)
        for h in hulls:
            covered |= _points_in_hull(h, pts)

        cell = (xs[1] - xs[0]) * (ys[1] - ys[0])
        return float(covered.sum() * cell)

    def shadowing_fraction(self, v_hat: Sequence[float], grid: int = 600) -> float:
        """How much of the summed area is double-counted overlap, 0 to 1.

        Zero when no two parts overlap in projection. This is the quantity the
        two solvers exist to expose; a design whose parts shadow each other
        heavily has an effective area the naive sum badly overstates.
        """
        summed = self.projected_area(v_hat)
        if summed <= 0.0:
            return 0.0
        return max(0.0, 1.0 - self.projected_area_union(v_hat, grid=grid) / summed)

    # ---- extremes -------------------------------------------------------

    def extreme_area(
        self,
        maximise: bool,
        use_union: bool = False,
        n_coarse: int = 4000,
        n_refine: int = 200,
        refine_rounds: int = 10,
    ) -> tuple[float, np.ndarray]:
        """Search the sphere of attitudes for the min or max projected area.

        Fibonacci-lattice coarse scan, then repeated local resampling in a
        shrinking spherical cap around the best direction. The function is
        piecewise-smooth in the direction and has no narrow spikes -- it is a
        sum of `|v.n|` terms -- so a coarse global scan followed by local
        refinement finds the extremum without needing a gradient.

        **This is a numerical search, not a closed form**, so the result
        carries a small error. With the defaults the cap shrinks to ~5e-6 rad,
        giving agreement with the analytic box extremes to better than 1e-4
        relative (`tests/test_geometry.py::test_box_extremes_match_the_closed_form`).
        Extrema at a sharp corner of the `|v.n|` surface -- which is where a
        box's minimum sits -- converge more slowly than smooth ones, so that
        test is the binding check on these defaults.

        Returns `(area, direction)`.
        """
        measure = (
            (lambda v: self.projected_area_union(v)) if use_union
            else (lambda v: self.projected_area(v))
        )
        sign = 1.0 if maximise else -1.0

        best_dir = None
        best_val = -math.inf
        for v in _fibonacci_hemisphere(n_coarse):
            val = sign * measure(v)
            if val > best_val:
                best_val, best_dir = val, v

        spread = math.pi / math.sqrt(n_coarse)
        for _ in range(refine_rounds):
            for v in _cap_samples(best_dir, spread, n_refine):
                val = sign * measure(v)
                if val > best_val:
                    best_val, best_dir = val, v
            spread *= 0.4

        return sign * best_val, best_dir

    def knife_edge_area(self, **kw) -> tuple[float, np.ndarray]:
        """Minimum area over all attitudes -- the low-drag orientation.

        For the February 2022 fleet this is the physically relevant one:
        SpaceX commanded knife-edge precisely to minimise this number
        (PHYSICS.md §5).
        """
        return self.extreme_area(maximise=False, **kw)

    def broadside_area(self, **kw) -> tuple[float, np.ndarray]:
        """Maximum area over all attitudes -- the worst-case drag orientation."""
        return self.extreme_area(maximise=True, **kw)

    def area_range_m2(self, **kw) -> tuple[float, float]:
        """`(knife_edge, broadside)`, the bounds any attitude must lie between."""
        return self.knife_edge_area(**kw)[0], self.broadside_area(**kw)[0]


def _fibonacci_hemisphere(n: int) -> Iterable[np.ndarray]:
    """Near-uniform directions over a hemisphere.

    A hemisphere suffices: projected area is invariant under `v -> -v`,
    because every term depends on `|v.n|`.
    """
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = (i + 0.5) / n            # 0..1, the upper hemisphere
        r = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * i
        yield np.array([r * math.cos(theta), r * math.sin(theta), z])


def _cap_samples(center: np.ndarray, spread: float, n: int) -> Iterable[np.ndarray]:
    """`n` directions within a spherical cap of angular radius `spread`."""
    e1, e2 = _projection_basis(center)
    rng = np.random.default_rng(0xC0FFEE)
    for _ in range(n):
        ang = spread * math.sqrt(rng.random())
        phi = 2.0 * math.pi * rng.random()
        yield unit(
            math.cos(ang) * center
            + math.sin(ang) * (math.cos(phi) * e1 + math.sin(phi) * e2)
        )


# --------------------------------------------------------------------------
# Building geometry from data/satellite_specs.json
# --------------------------------------------------------------------------

def _spec_value(block: dict, key: str, which: str) -> float:
    """Pull `value`, `range[0]` or `range[1]` for a spec field.

    `which` is "value", "min" or "max". A field with no `range` returns its
    `value` for all three, so sweeping a spec never silently invents a spread
    that the source did not claim.
    """
    field = block[key]
    if which == "value" or "range" not in field:
        return float(field["value"])
    lo, hi = field["range"]
    return float(lo if which == "min" else hi)


def geometry_from_spec(
    specs: dict,
    variant: str = "v1_5",
    which: str = "value",
    include_array: bool = True,
) -> Geometry:
    """Build the parametric geometry for a spec variant.

    Every dimension comes from `data/satellite_specs.json` with its own source
    and confidence tag; nothing is written into this function. `which`
    selects the low, nominal or high end of each swept dimension --
    V2_BRIEF.md §6's "uncertain values get swept, not chosen", applied to
    geometry.

    The bus lies in the x-y plane with its thickness along z, and the array is
    coplanar with it, extending along +x from the chassis edge. That is the
    Starlink flat-panel arrangement: a single wing unfurling in the plane of
    the chassis, not a pair of wings on a boom.

    `include_array=False` gives the chassis alone, which is what the
    `_finding` note in the spec file compares against Baruah's published
    maximum.
    """
    g = specs[variant]["geometry"]
    bl = _spec_value(g, "bus_length_m", which)
    bw = _spec_value(g, "bus_width_m", which)
    bt = _spec_value(g, "bus_thickness_m", which)

    parts: list[Part] = [Box(length_m=bl, width_m=bw, thickness_m=bt, name="bus")]

    if include_array:
        span = _spec_value(g, "solar_array_span_m", which)
        chord = _spec_value(g, "solar_array_chord_m", which)
        thick = _spec_value(g, "solar_array_thickness_m", which)
        # Coplanar with the bus, hinged at the +x edge, so the two never
        # overlap in a face-on projection.
        parts.append(
            Panel(
                span_m=span,
                chord_m=chord,
                normal=(0.0, 0.0, 1.0),
                center_m=(bl / 2 + span / 2, 0.0, 0.0),
                thickness_m=thick,
                name="solar_array",
            )
        )

    return Geometry(parts, label=f"{variant}:{which}{'' if include_array else ':bus-only'}")
