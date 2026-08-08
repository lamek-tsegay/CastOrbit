"""Stage-1 mass closure. V2_BRIEF.md §5 and §7, Phase 9.

Gate 9 is `test_gate_9_held_out_dry_mass`. It reports every held-out
prediction whether it passes or not, because V2_BRIEF.md §7 says to report the
misses rather than tune the table, and a test that only printed on success
would make that impossible to honour.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim.mass_model import (
    MassEstimationError,
    estimate_dry_mass,
    held_out_set,
    load_reference_satellites,
    population_scatter,
    score_held_out,
    training_set,
)

REFERENCE_PATH = Path(__file__).resolve().parent.parent / "data" / "reference_satellites.json"


@pytest.fixture(scope="module")
def sats():
    if not REFERENCE_PATH.exists():
        pytest.skip(f"{REFERENCE_PATH} not present")
    return load_reference_satellites()


# --------------------------------------------------------------------------
# The table
# --------------------------------------------------------------------------

def test_every_entry_is_sourced(sats):
    """V2_BRIEF.md §5: tag every output with its source, as satellite_specs.json does."""
    for s in sats:
        assert s.source, f"{s.id} has no source"
        assert s.mass_confidence in {
            "published", "derived", "estimated", "disputed"
        }, f"{s.id} has an unrecognised confidence tag"


def test_table_spans_smallsat_to_large_geo(sats):
    """The range V2_BRIEF.md §5 asks for, asserted rather than assumed."""
    masses = [s.mass_kg for s in sats]
    assert min(masses) < 10.0, "no cubesat-class entry"
    assert max(masses) > 3000.0, "no large-GEO-class entry"
    assert max(masses) / min(masses) > 100.0

    assert {s.altitude_class for s in sats} >= {"LEO", "GEO"}
    assert {s.payload_class for s in sats} >= {
        "earth_observation", "communications"
    }


def test_mass_kind_is_always_explicit(sats):
    """Confusing launch mass with dry mass would corrupt the table beyond use.

    GOES-16 is 5192 kg at launch and 2857 kg dry -- a 1.8x error if mixed up.
    """
    for s in sats:
        assert s.mass_kind in {"dry", "launch"}, f"{s.id}: {s.mass_kind!r}"


def test_launch_mass_entries_are_excluded_from_interpolation(sats):
    """They are kept in the table as evidence but must never set a dry mass."""
    launch_only = [s for s in sats if not s.is_dry]
    assert launch_only, "the test is vacuous if every entry is a dry mass"
    for s in launch_only:
        assert not s.usable_for_power_interpolation


def test_exactly_three_held_out(sats):
    ho = held_out_set(sats)
    assert len(ho) == 3
    # Selection criteria stated in the data file: one per size class, each with
    # published dry mass and published power.
    for s in ho:
        assert s.is_dry, f"{s.id} held out but has no dry mass to score against"
        assert s.power_w, f"{s.id} held out but has no power to predict from"
    masses = sorted(s.mass_kg for s in ho)
    assert masses[0] < 500 < masses[1] < 2000 < masses[2], (
        f"held-out set does not span size classes: {masses}"
    )


def test_held_out_are_invisible_to_the_interpolator(sats):
    held = {s.id for s in held_out_set(sats)}
    assert held
    assert not (held & {s.id for s in training_set(sats)})


# --------------------------------------------------------------------------
# GATE 9
# --------------------------------------------------------------------------

def test_gate_9_held_out_dry_mass(sats, capsys):
    """**Phase 9 gate: FAILED, 1 of 3 within 25%.** Recorded, not tuned away.

    This test is the standing record of that failure. It asserts the outcome
    that was actually measured -- one resolvable case inside the bar, two
    refused -- so a later change that appears to "fix" Gate 9 has to come here
    and say so explicitly rather than passing quietly.

    The response to the failure was to bound the method rather than extend it;
    `tests/test_bounded.py` covers that. The table was not adjusted.
    """
    rows = score_held_out(sats)
    answered = [r for r in rows if r["resolvable"]]
    refused = [r for r in rows if not r["resolvable"]]
    n_within = sum(bool(r["within_25pct"]) for r in answered)

    with capsys.disabled():
        print("\n  GATE 9 -- held-out dry mass, 25% bar: FAIL (1/3)")
        for r in rows:
            lo, hi = r["interval_kg"]
            got = (
                f"{r['predicted_kg']:6.0f} kg ({r['error_frac'] * 100:+.1f}%)"
                if r["resolvable"] else "refused, no point estimate"
            )
            print(f"    {r['name'][:26]:26s} actual {r['actual_kg']:6.0f} kg  "
                  f"{got}  range {lo:.0f}-{hi:.0f} kg")

    # Every estimate stays tagged and sourced whether or not it resolves.
    for r in rows:
        assert r["estimate"]["tag"] == "estimated"
        assert len(r["estimate"]["neighbours"]) == 2
        assert r["estimate"]["provenance"]

    assert n_within == 1, "Gate 9 is recorded as 1/3; this is the record"
    assert len(refused) == 2


def test_the_two_misses_became_refusals_not_wrong_numbers(sats, capsys):
    """A miss is only acceptable if its cause is identified. Both are.

      * PROBA-V: a 42x-wide power bracket, and 3.6x kg/W scatter among
        comparable spacecraft. Two real spacecraft at 320 and 330 W differ 2.2x
        in mass, so power alone cannot resolve it -- at any table density.
      * GOES-16: 0.71 kg/W against neighbours at 0.50 and 0.45. An unusually
        heavy platform for its power, not an arithmetic error.

    Under the bound both now decline to answer instead of answering wrongly.
    """
    rows = {r["id"]: r for r in score_held_out(sats)}

    proba = rows["proba_v"]
    assert not proba["resolvable"]
    assert proba["predicted_kg"] is None
    assert any(
        "apart in power" in r for r in proba["estimate"]["refusal_reasons"]
    ), "the wide bracket must be given as a reason on the estimate itself"

    goes = rows["goes_16"]
    assert not goes["resolvable"]
    actual_ratio = goes["actual_kg"] / goes["power_w"]
    neighbour_ratios = [n["kg_per_w"] for n in goes["estimate"]["neighbours"]]
    with capsys.disabled():
        print(f"\n  GOES-16 is {actual_ratio:.3f} kg/W; its neighbours are "
              f"{neighbour_ratios[0]:.3f} and {neighbour_ratios[1]:.3f} kg/W")
    assert actual_ratio > max(neighbour_ratios), (
        "GOES-16 should be heavier per watt than both neighbours -- that is "
        "the whole explanation for the miss"
    )


def test_population_scatter_bounds_any_power_only_predictor(sats, capsys):
    """The ceiling on this approach, measured rather than asserted.

    PROBA-V (320 W, 140 kg) and Deimos-2 (330 W, 310 kg) are real spacecraft
    at essentially identical power whose masses differ by 2.2x. No predictor
    taking power alone can tell them apart, so a 25% bar is not reachable at
    the small end by *any* interpolation scheme -- which is the honest reading
    of the PROBA-V miss.
    """
    s = population_scatter("earth_observation", sats)
    pair = s["closest_power_pair"]

    with capsys.disabled():
        print(f"\n  earth_observation kg/W spans {s['kg_per_w_min']:.3f}-"
              f"{s['kg_per_w_max']:.3f} ({s['kg_per_w_spread']:.1f}x)")
        print(f"  closest pair in power: {pair['a']} and {pair['b']} at "
              f"{pair['power_w'][0]:.0f}/{pair['power_w'][1]:.0f} W "
              f"-> {pair['mass_ratio']:.2f}x apart in mass")

    assert pair["power_gap"] < 1.1, "these should be at effectively the same power"
    assert pair["mass_ratio"] > 2.0, (
        "if the scatter has collapsed, the PROBA-V miss needs re-explaining"
    )


# --------------------------------------------------------------------------
# Provenance -- the non-negotiable part
# --------------------------------------------------------------------------

def test_every_estimate_names_the_spacecraft_it_came_from():
    """V2_BRIEF.md §5: 'If dry mass came from interpolating two real
    satellites, the UI says so.'"""
    est = estimate_dry_mass(power_w=1500.0, payload_class="earth_observation")
    assert est.tag == "estimated"
    assert len(est.neighbours) == 2
    for n in est.neighbours:
        assert n["name"] and n["source"]
    assert "interpolated between" in est.provenance
    d = est.as_dict()
    assert d["provenance"] and d["method"] and d["tag"] == "estimated"


def test_estimate_is_bracketed_by_its_neighbours():
    est = estimate_dry_mass(power_w=1500.0, payload_class="earth_observation")
    lo, hi = sorted(n["mass_kg"] for n in est.neighbours)
    assert lo <= est.mass_kg <= hi


def test_more_power_never_means_less_mass():
    """Monotonicity. A design tool that got this backwards would be worse than useless.

    Checked on the interval rather than the point estimate, since most powers
    now decline to produce one -- and the interval is the thing a caller always
    gets, so it is the thing that has to behave.
    """
    intervals = [
        estimate_dry_mass(p, "earth_observation").interval_kg
        for p in (100.0, 500.0, 1500.0, 3000.0)
    ]
    assert [lo for lo, _ in intervals] == sorted(lo for lo, _ in intervals)
    assert [hi for _, hi in intervals] == sorted(hi for _, hi in intervals)


def test_estimates_that_do_resolve_are_still_monotone():
    """Where a point estimate survives the bound, it must order correctly too."""
    points = []
    for p in (900.0, 1200.0, 1700.0, 2400.0):
        est = estimate_dry_mass(p, "earth_observation")
        if est.resolvable:
            points.append(est.mass_kg)
    assert len(points) >= 2, "no resolvable cases left to check ordering on"
    assert points == sorted(points)


def test_extrapolation_is_flagged_not_silent():
    est = estimate_dry_mass(power_w=50000.0, payload_class="communications")
    assert est.extrapolated
    assert any("outside the reference range" in w for w in est.warnings)


def test_family_exclusion_prevents_predicting_from_a_twin(sats):
    """Without it, GOES-16 would be 'predicted' from GOES-17 and mean nothing."""
    est = estimate_dry_mass(
        power_w=4000.0, payload_class="earth_observation",
        satellites=sats, exclude_family="himawari",
    )
    assert all(n["id"] != "himawari_8" for n in est.neighbours)


def test_thin_payload_class_falls_back_loudly():
    est = estimate_dry_mass(power_w=800.0, payload_class="technology")
    assert any("widened to all classes" in w for w in est.warnings)


def test_thin_class_can_be_made_an_error_instead():
    with pytest.raises(MassEstimationError, match="need 2"):
        estimate_dry_mass(
            power_w=800.0, payload_class="technology", allow_class_fallback=False
        )


def test_rejects_nonsense_power():
    with pytest.raises(ValueError, match="power must be positive"):
        estimate_dry_mass(power_w=0.0, payload_class="earth_observation")


def test_no_mass_is_produced_without_a_table():
    """Deliberately an error, not a fallback guess. Nothing faked."""
    with pytest.raises(MassEstimationError):
        estimate_dry_mass(
            power_w=1000.0, payload_class="earth_observation", satellites=[]
        )


def test_reference_file_is_valid_json_with_documented_conventions():
    with open(REFERENCE_PATH, encoding="utf-8") as fh:
        payload = json.load(fh)
    for key in ("_mass_convention", "_power_convention", "_held_out", "_family_note"):
        assert key in payload, f"{key} missing -- the conventions must be written down"
    assert len(payload["satellites"]) >= 15
