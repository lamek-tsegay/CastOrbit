"""JSON export and the Phase 4 read-back gate.

ARCHITECTURE.md §5, Phase 4 gate: a JSON file exists that fully describes a
run, and a script can read it back and reproduce the plots -- using only the
file, not the in-memory objects that produced it. `replot_batch_from_json` and
`replot_sweeps_from_json` are exactly that second script; these tests exercise
them against freshly-written JSON rather than trusting a stale file in out/.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim.atmosphere import SpaceWeather
from sim.export import export_fleet_batch, load_json, replot_batch_from_json

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"


@pytest.fixture(scope="module")
def sw() -> SpaceWeather:
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return SpaceWeather.load(DATA)


@pytest.fixture(scope="module")
def batch_payload(sw, tmp_path_factory) -> dict:
    path = tmp_path_factory.mktemp("export") / "batch.json"
    export_fleet_batch(sw, path=path, dt=60.0)
    return load_json(path)


def test_batch_json_has_both_fleet_runs(batch_payload):
    """README.md's central finding needs both Cd conventions in one file."""
    labels = {r["label"]: r for r in batch_payload["runs"]}
    assert set(labels) == {"cd2.2", "cd1.0"}
    assert labels["cd2.2"]["config"]["cd"] == 2.2
    assert labels["cd1.0"]["config"]["cd"] == 1.0


def test_batch_json_matches_the_observed_event(batch_payload):
    """data/event_feb2022.json: 38 lost, 11 survived -- carried through as-is."""
    assert batch_payload["observed"] == {
        "lost": 38, "survived": 11, "source": "data/event_feb2022.json",
    }


def test_every_run_accounts_for_all_49_satellites(batch_payload):
    for run in batch_payload["runs"]:
        assert len(run["satellites"]) == 49
        assert sum(run["outcome_counts"].values()) == 49


def test_cd2_2_loses_the_whole_fleet(batch_payload):
    """Pins the headline fleet-reproduction result behind the JSON export."""
    cd22 = next(r for r in batch_payload["runs"] if r["label"] == "cd2.2")
    assert cd22["outcome_counts"]["REENTERED"] == 49


def test_cd1_0_brackets_the_observed_loss_count(batch_payload):
    cd10 = next(r for r in batch_payload["runs"] if r["label"] == "cd1.0")
    lost = cd10["outcome_counts"]["REENTERED"]
    assert 30 <= lost <= 42, f"{lost} lost is outside the range seen at Gate 3"


def test_every_satellite_is_tagged_with_cd_times_area(batch_payload):
    """The whole point of this file: Cd*A travels with each satellite.

    Not just present -- numerically consistent with that satellite's own cd
    and ram_area_m2, and inside the run's stated effective_drag_range_m2.
    """
    for run in batch_payload["runs"]:
        cd = run["config"]["cd"]
        lo, hi = run["config"]["effective_drag_range_m2"]
        for sat in run["satellites"]:
            p = sat["params"]
            assert p["cd_times_area_m2"] == pytest.approx(
                cd * p["ram_area_m2"], abs=1e-3
            )
            assert lo - 1e-6 <= p["cd_times_area_m2"] <= hi + 1e-6


def test_effective_drag_ranges_match_the_readme_finding(batch_payload):
    """README.md: Cd*A = 2.2*[1.00,4.48] and Baruah's own 1.0*[1.00,4.48]."""
    runs = {r["label"]: r for r in batch_payload["runs"]}
    assert runs["cd1.0"]["config"]["effective_drag_range_m2"] == [1.0, 4.48]
    lo, hi = runs["cd2.2"]["config"]["effective_drag_range_m2"]
    assert lo == pytest.approx(2.2, abs=1e-6)
    assert hi == pytest.approx(2.2 * 4.48, abs=1e-6)


def test_json_is_plain_data(batch_payload):
    """Round-trips through json.dumps with no numpy types, enums, or NaN."""
    text = json.dumps(batch_payload, allow_nan=False)
    assert json.loads(text) == batch_payload


def test_read_back_reproduces_the_altitude_chart(batch_payload, tmp_path):
    """Phase 4 gate: replot from the file alone, not the objects that made it."""
    path = tmp_path / "batch.json"
    path.write_text(json.dumps(batch_payload))

    out_png = replot_batch_from_json(path)

    assert out_png.exists()
    assert out_png.stat().st_size > 10_000, "suspiciously small for a 2-panel figure"


def test_replot_sweeps_from_json_uses_only_the_file(tmp_path):
    """Same gate for sweeps.json, against a minimal synthetic payload.

    Deliberately does not depend on a full sim.sweeps.main() run (minutes, not
    seconds) -- this checks that plot_from_payload only reads what's in the
    dict, matching what sim/sweeps.py's own docstring promises.
    """
    from sim.export import replot_sweeps_from_json

    payload = {
        "meta": {"n_satellites": 49, "actual_insertion_km": 210.0},
        "critical_altitude_km": {"storm": [181.0, 188.0], "quiet": [179.0, 185.0]},
        "sweep_insertion_altitude": {
            "x_km": [180, 190, 200],
            "survival": {
                "storm|1.00": [0.2, 0.6, 1.0], "storm|1.19": [0.1, 0.5, 0.9],
                "quiet|1.00": [0.3, 0.7, 1.0], "quiet|1.19": [0.2, 0.6, 0.95],
            },
        },
        "sweep_ram_area": {
            "x_m2": [1.0, 4.48, 6.0],
            "survival": {
                "storm|1.00": [1.0, 1.0, 0.9], "storm|1.19": [1.0, 0.95, 0.3],
                "quiet|1.00": [1.0, 1.0, 1.0], "quiet|1.19": [1.0, 1.0, 0.7],
            },
            "published_range_m2": [1.0, 4.48],
            "source_ranges_m2": [],
        },
        "sweep_safe_mode_timing": {
            "x_hours": [0, 24, 48],
            "survival": {
                "storm|1.00": [1.0, 0.8, 0.5], "storm|1.19": [1.0, 0.7, 0.3],
                "quiet|1.00": [1.0, 0.9, 0.7], "quiet|1.19": [1.0, 0.8, 0.5],
            },
        },
    }
    path = tmp_path / "sweeps.json"
    path.write_text(json.dumps(payload))

    paths = replot_sweeps_from_json(path)

    assert len(paths) == 3
    for p in paths:
        assert p.exists()
        assert p.stat().st_size > 5_000
