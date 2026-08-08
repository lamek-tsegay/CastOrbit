// Display-only helpers: string/date formatting and reading already-computed
// sample arrays at a given time. Nothing here derives a new physical
// quantity -- every number consumed or returned already exists in the JSON.

// JS drops trailing zeros on plain numbers (1.0 -> "1"), which reads as a
// typo next to "Cd = 2.2". Every Cd in this project is quoted to one decimal
// (data/event_feb2022.json, satellite_specs.json), so format it that way.
export function formatCd(cd) {
  return Number(cd).toFixed(1);
}

export function formatUTC(isoString) {
  if (!isoString) return "—";
  const d = new Date(isoString);
  return d.toISOString().replace("T", " ").replace(/:\d{2}\.\d{3}Z$/, " UTC");
}

export function formatHours(seconds) {
  const h = seconds / 3600;
  return `${h >= 0 ? "+" : ""}${h.toFixed(2)} h`;
}

export function isoPlusSeconds(isoEpoch, seconds) {
  const d = new Date(isoEpoch);
  d.setUTCSeconds(d.getUTCSeconds() + seconds);
  return d.toISOString();
}

/**
 * Linear interpolation between the two trajectory samples bracketing `t`.
 *
 * `trajectory` is exactly one satellite's exported `trajectory` object
 * (t_s, h_km, rho, lat_deg, lon_deg -- all already computed in Python).
 * This only reads between two already-known samples for smooth playback; it
 * does not compute anything the export didn't already produce, the same way
 * a chart library interpolates between plotted points on a line.
 *
 * Returns null if `t` is before the trajectory starts (satellite not yet
 * deployed) or after its last sample (already resolved/truncated).
 */
export function sampleTrajectoryAt(trajectory, t) {
  const { t_s } = trajectory;
  if (!t_s || t_s.length === 0) return null;
  if (t < t_s[0]) return null;
  if (t >= t_s[t_s.length - 1]) {
    return pickSample(trajectory, t_s.length - 1);
  }

  // t_s is sorted ascending (Python export order); linear scan is fine at
  // ~200-600 samples per satellite.
  let i = 0;
  while (i < t_s.length - 1 && t_s[i + 1] < t) i++;

  const t0 = t_s[i];
  const t1 = t_s[i + 1];
  const frac = t1 > t0 ? (t - t0) / (t1 - t0) : 0;

  const s0 = pickSample(trajectory, i);
  const s1 = pickSample(trajectory, i + 1);
  const lerp = (a, b) => a + (b - a) * frac;

  return {
    t_s: t,
    h_km: lerp(s0.h_km, s1.h_km),
    rho: lerp(s0.rho, s1.rho),
    lat_deg: lerp(s0.lat_deg, s1.lat_deg),
    lon_deg: lerpAngle(s0.lon_deg, s1.lon_deg, frac),
  };
}

function pickSample(trajectory, i) {
  return {
    t_s: trajectory.t_s[i],
    h_km: trajectory.h_km[i],
    rho: trajectory.rho ? trajectory.rho[i] : undefined,
    lat_deg: trajectory.lat_deg ? trajectory.lat_deg[i] : undefined,
    lon_deg: trajectory.lon_deg ? trajectory.lon_deg[i] : undefined,
  };
}

// Longitude wraps at +/-180; a plain lerp across the seam would swing the
// long way around. Take the shorter arc instead.
function lerpAngle(a, b, frac) {
  if (a === undefined || b === undefined) return undefined;
  let diff = b - a;
  if (diff > 180) diff -= 360;
  if (diff < -180) diff += 360;
  let v = a + diff * frac;
  if (v > 180) v -= 360;
  if (v < -180) v += 360;
  return v;
}
