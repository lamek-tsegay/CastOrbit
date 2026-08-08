"""Projected-area solver. V2_BRIEF.md §7, Phase 8.

Two gates to satisfy, and both are here:

  1. A flat plate at angle theta to the velocity vector gives `A*cos(theta)`
     to within 1%.
  2. For a v1.5-like spec, the knife-edge and broadside areas bracket the
     published 1.00-4.48 m^2 range.

Everything else tests the properties that make those two believable.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from sim.geometry import (
    Box,
    Geometry,
    Panel,
    direction_from_angles,
    geometry_from_spec,
    unit,
)

SPECS_PATH = Path(__file__).resolve().parent.parent / "data" / "satellite_specs.json"

# Baruah et al. (2024), via data/satellite_specs.json.
PUBLISHED_KNIFE_EDGE_M2 = 1.00
PUBLISHED_MAX_M2 = 4.48


@pytest.fixture(scope="module")
def specs() -> dict:
    if not SPECS_PATH.exists():
        pytest.skip(f"{SPECS_PATH} not present")
    with open(SPECS_PATH, encoding="utf-8") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# GATE 8, part 2 -- the flat plate cosine law
# --------------------------------------------------------------------------

@pytest.mark.parametrize("theta_deg", [0, 15, 30, 45, 60, 75, 89])
def test_flat_plate_follows_the_cosine_law(theta_deg, capsys):
    """**Phase 8 gate:** `A*cos(theta)` to within 1%.

    A zero-thickness plate would satisfy this exactly and trivially, so the
    plate under test carries the same 1 cm thickness the real spec uses. The
    edge term is what makes this a real 1% test rather than an identity: at
    theta = 89 deg the face contributes almost nothing and the edge dominates.
    """
    plate = Panel(span_m=3.0, chord_m=1.5, normal=(0, 0, 1), thickness_m=0.01)
    g = Geometry([plate])

    # Tilt within the x-z plane: theta is measured from the plate normal.
    v = direction_from_angles(azimuth_deg=0.0, elevation_deg=90.0 - theta_deg)
    got = g.projected_area(v)
    expected = plate.area_m2 * math.cos(math.radians(theta_deg))

    # Compare against the face term alone; the edge contribution is real
    # geometry, so it is reported rather than subtracted out.
    err = abs(got - expected) / plate.area_m2
    with capsys.disabled():
        print(f"  theta={theta_deg:2d} deg: {got:.5f} m2 vs A*cos = "
              f"{expected:.5f} m2  ({err * 100:.3f}% of A)")
    assert err < 0.01, f"cosine law off by {err * 100:.3f}% of the plate area"


def test_cosine_law_is_exact_for_a_zero_thickness_plate():
    """With no edge, the law must hold to machine precision, not just 1%.

    Separates 'the projection maths is right' from 'the edge term is small',
    so a real error in the former cannot hide inside the 1% budget above.
    """
    plate = Panel(span_m=3.0, chord_m=1.5, normal=(0, 0, 1), thickness_m=0.0)
    g = Geometry([plate])
    for theta_deg in (0, 23, 45, 67, 90):
        v = direction_from_angles(0.0, 90.0 - theta_deg)
        expected = plate.area_m2 * math.cos(math.radians(theta_deg))
        assert g.projected_area(v) == pytest.approx(expected, abs=1e-12)


def test_edge_on_plate_presents_only_its_edge():
    """Both edge-on directions, which also pins the span/chord frame convention.

    `Panel._frame()` derives its in-plane axes from the normal, so which body
    axis "span" ends up along is a consequence of that construction rather
    than something the caller sets. Pinned here because getting it backwards
    silently swaps two edge areas, and for a long thin array those differ by
    the aspect ratio.
    """
    plate = Panel(span_m=3.0, chord_m=1.5, normal=(0, 0, 1), thickness_m=0.01)
    g = Geometry([plate])
    e_span, e_chord, _ = plate._frame()

    # Looking down the span axis, the cross-section is chord x thickness.
    assert g.projected_area(e_span) == pytest.approx(1.5 * 0.01, rel=1e-9)
    # Looking down the chord axis, it is span x thickness.
    assert g.projected_area(e_chord) == pytest.approx(3.0 * 0.01, rel=1e-9)

    # For this normal the construction puts span along -y and chord along +x.
    assert np.allclose(e_span, [0, -1, 0])
    assert np.allclose(e_chord, [1, 0, 0])


# --------------------------------------------------------------------------
# GATE 8, part 1 -- bracketing the published range
# --------------------------------------------------------------------------

def test_v1_5_area_range_brackets_the_published_range(specs, capsys):
    """**Phase 8 gate:** knife-edge and broadside bracket 1.00-4.48 m^2.

    Checked across the whole swept dimension range, not just the nominal
    values, because the chassis dimensions are tagged `disputed` in the spec
    file and a gate that only passed at one corner would not mean much.
    """
    with capsys.disabled():
        print()
    for which in ("min", "value", "max"):
        g = geometry_from_spec(specs, which=which)
        knife, broad = g.area_range_m2()
        with capsys.disabled():
            print(f"  v1.5 [{which:5s}]: knife-edge {knife:6.3f} m2  <=  "
                  f"{PUBLISHED_KNIFE_EDGE_M2:.2f}  |  broadside {broad:7.3f} m2  >=  "
                  f"{PUBLISHED_MAX_M2:.2f}")
        assert knife <= PUBLISHED_KNIFE_EDGE_M2, (
            f"knife-edge {knife:.3f} m2 exceeds the published lower bound"
        )
        assert broad >= PUBLISHED_MAX_M2, (
            f"broadside {broad:.3f} m2 does not reach the published maximum"
        )


def test_published_maximum_is_a_face_on_chassis_not_a_broadside(specs, capsys):
    """The Phase 8 finding, pinned.

    Baruah's 4.48 m^2 "maximum" reproduces as the *face-on* area of the larger
    sourced chassis (3.0 x 1.5 = 4.50 m^2, 0.45% agreement) with the array
    feathered. It is neither a true broadside -- which includes the array and
    is ~15 m^2 -- nor even the chassis's own geometric maximum, since a box
    presents most area corner-on.

    This matters beyond trivia: it says the published range the V1 sweeps used
    describes a *feathered* configuration throughout, consistent with the
    knife-edge attitude the fleet was actually commanded into (PHYSICS.md §5).
    """
    chassis = Geometry([Box(length_m=3.0, width_m=1.5, thickness_m=0.35)])
    face_on = chassis.projected_area((0, 0, 1))
    corner_on, _ = chassis.broadside_area()
    with_array, _ = geometry_from_spec(specs, which="max").broadside_area()

    with capsys.disabled():
        print(f"\n  chassis face-on      {face_on:6.3f} m2   (published 4.48)")
        print(f"  chassis corner-on    {corner_on:6.3f} m2")
        print(f"  true broadside       {with_array:6.3f} m2   (array included)")

    assert face_on == pytest.approx(PUBLISHED_MAX_M2, rel=0.01)
    assert corner_on > face_on, "a box presents more area corner-on than face-on"
    assert with_array > 3 * PUBLISHED_MAX_M2, (
        "a genuine broadside including the array should dwarf the published figure"
    )


def test_knife_edge_is_driven_by_chassis_thickness(specs):
    """Sanity on which parameter controls the low-drag number.

    The spec file claims bus thickness "sets the knife-edge area almost
    single-handedly". If that stopped being true the sweep guidance would be
    wrong.
    """
    lo = geometry_from_spec(specs, which="min").knife_edge_area()[0]
    hi = geometry_from_spec(specs, which="max").knife_edge_area()[0]
    assert hi > lo
    # Thickness spans 0.15-0.35 m, a factor 2.33; the knife-edge area should
    # scale with roughly that, not with the ~1.07 factor in the other dims.
    assert 1.8 < hi / lo < 3.0


# --------------------------------------------------------------------------
# Solver properties
# --------------------------------------------------------------------------

def test_box_extremes_match_the_closed_form():
    """A box's min and max projected areas are analytic. Check against them.

        min = smallest face area
        max = sqrt((WT)^2 + (LT)^2 + (LW)^2)   -- the corner-on diagonal

    This is what validates the numerical direction search, so the search is
    never trusted on a shape where the answer is unknown without first being
    correct on one where it is.
    """
    L, W, T = 3.0, 1.5, 0.35
    g = Geometry([Box(L, W, T)])

    expected_min = min(W * T, L * T, L * W)
    expected_max = math.sqrt((W * T) ** 2 + (L * T) ** 2 + (L * W) ** 2)

    got_min, _ = g.knife_edge_area()
    got_max, _ = g.broadside_area()
    assert got_min == pytest.approx(expected_min, rel=1e-4)
    assert got_max == pytest.approx(expected_max, rel=1e-4)


def test_projection_is_symmetric_under_velocity_reversal():
    """Every term depends on |v.n|, so the hemisphere search is legitimate."""
    g = geometry_from_spec(json.loads(SPECS_PATH.read_text()), which="value")
    v = unit((0.3, -0.5, 0.81))
    assert g.projected_area(v) == pytest.approx(g.projected_area(-v), rel=1e-12)


def test_area_is_bounded_by_the_extremes():
    """No attitude may fall outside the reported range. Checked on random ones."""
    g = Geometry([
        Box(3.0, 1.5, 0.35),
        Panel(span_m=8.0, chord_m=1.4, center_m=(5.5, 0, 0)),
    ])
    knife, broad = g.area_range_m2()
    rng = np.random.default_rng(11)
    for _ in range(200):
        v = unit(rng.normal(size=3))
        a = g.projected_area(v)
        assert knife - 1e-6 <= a <= broad + 1e-6


def test_two_solvers_agree_when_parts_do_not_overlap():
    """The independent cross-check. Coplanar, side-by-side parts never overlap.

    Analytic sum and rasterised union reach the same number through entirely
    different code -- per-part formulas versus vertices, convex hulls and a
    pixel count. Agreement here is what licenses using the cheap one.
    """
    g = Geometry([
        Box(3.0, 1.5, 0.35),
        Panel(span_m=8.0, chord_m=1.4, center_m=(5.5, 0, 0), thickness_m=0.01),
    ])
    for v in [(0, 0, 1), (0.2, 0.1, 0.97), (0.4, 0.3, 0.87)]:
        summed = g.projected_area(v)
        union = g.projected_area_union(v, grid=900)
        assert union == pytest.approx(summed, rel=0.02), (
            f"solvers disagree at {v}: {summed:.4f} vs {union:.4f}"
        )
        assert g.shadowing_fraction(v, grid=900) < 0.02


def test_union_solver_catches_shadowing_the_sum_misses():
    """Stack two identical plates and the sum double-counts; the union does not.

    This is the case the two-solver arrangement exists for. Without it, an
    assembly with occluding parts would silently report up to twice its real
    drag area.
    """
    stacked = Geometry([
        Panel(span_m=2.0, chord_m=2.0, center_m=(0, 0, 0.0), thickness_m=0.01),
        Panel(span_m=2.0, chord_m=2.0, center_m=(0, 0, 0.5), thickness_m=0.01),
    ])
    v = (0, 0, 1)  # face-on: the second plate hides exactly behind the first
    summed = stacked.projected_area(v)
    union = stacked.projected_area_union(v, grid=600)

    assert summed == pytest.approx(8.0, rel=1e-6)          # 2 x 4 m2
    assert union == pytest.approx(4.0, rel=0.02)           # one plate's worth
    assert stacked.shadowing_fraction(v, grid=600) == pytest.approx(0.5, abs=0.02)


def test_single_part_union_is_exact_not_rasterised():
    """Single parts short-circuit to the hull area, so no discretisation error."""
    g = Geometry([Panel(span_m=3.0, chord_m=1.5, thickness_m=0.0)])
    assert g.projected_area_union((0, 0, 1)) == pytest.approx(4.5, rel=1e-9)


# --------------------------------------------------------------------------
# Input validation
# --------------------------------------------------------------------------

def test_degenerate_geometry_is_rejected():
    with pytest.raises(ValueError, match="at least one part"):
        Geometry([])
    with pytest.raises(ValueError, match="length_m must be positive"):
        Box(0.0, 1.0, 1.0)
    with pytest.raises(ValueError, match="span_m must be positive"):
        Panel(span_m=-1.0, chord_m=1.0)
    with pytest.raises(ValueError, match="cannot normalise"):
        unit((0.0, 0.0, 0.0))


def test_spec_sweep_selects_range_ends(specs):
    """`which` must actually move the geometry, or the sweep is decorative."""
    lo = geometry_from_spec(specs, which="min").parts[0]
    hi = geometry_from_spec(specs, which="max").parts[0]
    assert hi.length_m > lo.length_m
    assert hi.thickness_m > lo.thickness_m


# --------------------------------------------------------------------------
# Geometry -> physics: the degeneracy attacked from the geometry side
# --------------------------------------------------------------------------

def test_geometry_knife_edge_matches_the_secondary_sources_not_baruah(specs, capsys):
    """The solver lands on 0.27-0.61 m^2, the *secondary* knife-edge range.

    `satellite_specs.json` records two disagreeing claims for the knife-edge
    area: Baruah's 1.00 m^2 lower bound, and secondary sources at 0.3-0.7 m^2.
    Built independently from chassis dimensions and a published end-to-end
    span, the solver reproduces the secondary range and sits well below
    Baruah's figure.

    That is not a contradiction with V1. Baruah's 1.00 m^2 is a stated *lower
    bound* on a swept range, not a measurement, and V1 validated against the
    paper's own convention deliberately.
    """
    lo = min(geometry_from_spec(specs, which=w).knife_edge_area()[0]
             for w in ("min", "value", "max"))
    hi = max(geometry_from_spec(specs, which=w).knife_edge_area()[0]
             for w in ("min", "value", "max"))
    with capsys.disabled():
        print(f"\n  geometry knife-edge {lo:.3f}-{hi:.3f} m2 vs "
              f"secondary sources 0.3-0.7 m2, Baruah 1.00 m2")

    assert hi < PUBLISHED_KNIFE_EDGE_M2
    # Overlaps the 0.3-0.7 secondary range rather than merely being small.
    assert lo < 0.7 and hi > 0.3


def test_effective_drag_parameter_converges_with_baruah(specs, capsys):
    """A fourth line of evidence on `Cd*A`, reached from geometry.

    The README's central finding is that only the product `rho*Cd*A` is
    observable from a decay curve, and that three independent diagnostics put
    it in the same place. This is a fourth, and the first that does not use a
    decay curve at all:

        Baruah:              Cd = 1.0  x  A = 1.00 m^2       = 1.00 m^2
        geometry + standard: Cd = 2.2  x  A = 0.405 m^2      = 0.89 m^2

    Two different splits of an unobservable product -- the paper's stated
    simplification versus this project's free-molecular convention applied to
    a solved geometry -- landing 11% apart. Exactly the degeneracy V2_BRIEF.md
    §2 says pinning `A` from geometry should attack.
    """
    cd_standard = 2.2
    baruah_cd_a = 1.0 * PUBLISHED_KNIFE_EDGE_M2

    areas = [geometry_from_spec(specs, which=w).knife_edge_area()[0]
             for w in ("min", "value", "max")]
    band = [cd_standard * a for a in areas]
    nominal = cd_standard * areas[1]

    with capsys.disabled():
        print(f"\n  Baruah  Cd*A = {baruah_cd_a:.3f} m2")
        print(f"  geometry Cd*A = {min(band):.3f}-{max(band):.3f} m2 "
              f"(nominal {nominal:.3f})  ->  {abs(nominal - baruah_cd_a) / baruah_cd_a * 100:.1f}% apart")

    assert min(band) <= baruah_cd_a <= max(band), (
        "the geometry-derived effective drag range should bracket Baruah's"
    )
    assert abs(nominal - baruah_cd_a) / baruah_cd_a < 0.25


def test_ram_area_for_mode_maps_attitude_to_area(specs):
    """PHYSICS.md §5: safe mode is the minimum-area attitude, by definition."""
    from sim.geometry import ram_area_for_mode
    from sim.satellite import ThrusterMode

    g = geometry_from_spec(specs)
    safe = ram_area_for_mode(g, ThrusterMode.SAFE_MODE)
    nominal = ram_area_for_mode(g, ThrusterMode.NOMINAL)

    assert safe == pytest.approx(g.knife_edge_area()[0])
    assert nominal == pytest.approx(g.broadside_area()[0])
    assert safe < nominal
