// Reshapes already-computed batch.json trajectories into the row-per-time
// shape recharts wants (one object per sample, one key per series). No value
// is computed here that wasn't already in the export -- every satellite's
// t_s array is a strict prefix of the run's shared, untruncated time grid
// (sim/export.py's `t_ds`), because it was built by slicing that same array
// at the satellite's own outcome time. That means index i always means the
// same absolute time across every satellite and the critical-altitude
// series in a given run, so this is a straight re-index, not an
// interpolation or a new number.
export function buildAltitudeRows(run) {
  const sats = run.satellites;
  let canonical = sats[0]?.trajectory.t_s ?? [];
  for (const s of sats) {
    if (s.trajectory.t_s.length > canonical.length) canonical = s.trajectory.t_s;
  }

  const critT = run.critical_altitude.times_s ?? null;
  const critH = run.critical_altitude.h_crit_km;
  const n = Math.max(canonical.length, critH?.length ?? 0);

  const rows = [];
  for (let i = 0; i < n; i++) {
    const row = {
      t_s: i < canonical.length ? canonical[i] : critT?.[i] ?? null,
      h_crit: i < (critH?.length ?? 0) ? critH[i] : null,
    };
    for (const s of sats) {
      row[`sat_${s.id}`] =
        i < s.trajectory.h_km.length ? s.trajectory.h_km[i] : null;
    }
    rows.push(row);
  }
  return rows;
}

// critical_altitude.times comes as ISO timestamps (sim/export.py); recharts
// needs a numeric x-axis to share with t_s. Converting ISO -> seconds-since-
// epoch is arithmetic on already-known values, not a new physical quantity.
export function withCriticalAltitudeSeconds(run, epochIso) {
  const epochMs = new Date(epochIso).getTime();
  return {
    ...run,
    critical_altitude: {
      ...run.critical_altitude,
      times_s: run.critical_altitude.times.map(
        (iso) => (new Date(iso).getTime() - epochMs) / 1000,
      ),
    },
  };
}
