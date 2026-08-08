import { useEffect, useMemo, useRef, useState } from "react";
import Globe from "react-globe.gl";
import OutcomeLegend from "../components/OutcomeLegend";
import RunSelector from "../components/RunSelector";
import { usePlayback } from "../context/PlaybackContext";
import { formatCd, sampleTrajectoryAt } from "../lib/format";
import { outcomeColor } from "../lib/outcomeColors";

// react-globe.gl's pointAltitude is already a fraction of Earth's radius, so
// h_km / EARTH_RADIUS_KM alone gives a correctly-scaled bump (~0.03 for a
// 200 km orbit) with no exaggeration needed. This is a rendering choice
// (how tall to draw the marker), not a physics computation, and it is
// applied identically to every point regardless of outcome, run, or
// altitude, so it never changes which points are higher or lower than
// which others.
const DISPLAY_ALTITUDE_SCALE = 1;
const EARTH_RADIUS_KM = 6371;

export default function GlobeView({ batch }) {
  const { currentT } = usePlayback();
  const [selectedLabel, setSelectedLabel] = useState(null);
  const [clicked, setClicked] = useState(null);
  const [dims, setDims] = useState({ width: 600, height: 520 });
  const wrapRef = useRef(null);

  useEffect(() => {
    if (!wrapRef.current) return undefined;
    const el = wrapRef.current;
    const ro = new ResizeObserver((entries) => {
      const { width, height } = entries[0].contentRect;
      if (width > 0 && height > 0) setDims({ width, height });
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  const run = useMemo(
    () => batch.runs.find((r) => r.label === selectedLabel) ?? null,
    [batch, selectedLabel],
  );

  const points = useMemo(() => {
    if (!run) return [];
    return run.satellites
      .map((sat) => {
        const s = sampleTrajectoryAt(sat.trajectory, currentT);
        if (!s || s.lat_deg === undefined) return null;
        return {
          id: sat.id,
          lat: s.lat_deg,
          lng: s.lon_deg,
          altitude: (s.h_km / EARTH_RADIUS_KM) * DISPLAY_ALTITUDE_SCALE,
          h_km: s.h_km,
          color: outcomeColor(sat.outcome),
          outcome: sat.outcome,
          outcomeTime: sat.outcome_time,
          cause: sat.cause,
          params: sat.params,
        };
      })
      .filter(Boolean);
  }, [run, currentT]);

  return (
    <>
      <div className="view-header">
        <span className="view-title">Globe</span>
        <span className="view-note">
          Ground-track position (latitude/longitude) is illustrative display
          geometry only -- computed from each satellite's own simulated
          altitude plus an illustrative orbital-plane spread, not measured or
          tracked. Lat/lon modelling is outside the validated physics (see
          PHYSICS.md §1 and sim/groundtrack.py). Altitude, outcome, and cause
          of loss are real simulation output.
        </span>
      </div>

      <RunSelector
        runs={batch.runs}
        selected={selectedLabel}
        onSelect={(label) => {
          setSelectedLabel(label);
          setClicked(null);
        }}
        nSatellites={batch.meta.n_satellites}
      />

      {!run ? (
        <div className="panel" style={{ textAlign: "center", padding: "48px 16px" }}>
          Choose a Cd convention above to render the fleet.
        </div>
      ) : (
        <div className="globe-layout">
          <div className="globe-canvas-wrap" ref={wrapRef}>
            <Globe
              width={dims.width}
              height={dims.height}
              backgroundColor="#000000"
              globeImageUrl="/textures/earth-dark.jpg"
              showAtmosphere={true}
              atmosphereColor="#4fa3ff"
              atmosphereAltitude={0.12}
              pointsData={points}
              pointLat="lat"
              pointLng="lng"
              pointAltitude="altitude"
              pointColor="color"
              pointRadius={0.35}
              pointLabel={(p) => `#${p.id} · ${p.outcome} · ${p.h_km.toFixed(1)} km`}
              onPointClick={(p) => setClicked(p)}
            />
          </div>

          <div className="panel globe-detail">
            <div className="panel-title">
              {clicked ? `Satellite #${clicked.id}` : "Click a satellite"}
            </div>
            {clicked ? (
              <dl>
                <dt>Outcome</dt>
                <dd style={{ color: clicked.color }}>{clicked.outcome}</dd>
                <dt>Current altitude</dt>
                <dd>{clicked.h_km.toFixed(2)} km</dd>
                <dt>Mass</dt>
                <dd>{clicked.params.mass_kg.toFixed(1)} kg</dd>
                <dt>Ram area</dt>
                <dd>{clicked.params.ram_area_m2.toFixed(3)} m²</dd>
                <dt>Cd</dt>
                <dd>{clicked.params.cd}</dd>
                <dt>Cd·A (effective drag)</dt>
                <dd>{clicked.params.cd_times_area_m2.toFixed(3)} m²</dd>
                <dt>Deploy offset</dt>
                <dd>{clicked.params.deploy_time_s.toFixed(0)} s</dd>
                <dt>Cause</dt>
                <dd className="cause">{clicked.cause}</dd>
              </dl>
            ) : (
              <p style={{ fontSize: 12, color: "var(--text-dim)" }}>
                Click any point on the globe for its parameters and cause of
                loss.
              </p>
            )}
            <div className="panel-title" style={{ marginTop: 16 }}>
              Run: Cd = {formatCd(run.config.cd)}
            </div>
            <OutcomeLegend counts={run.outcome_counts} />
          </div>
        </div>
      )}
    </>
  );
}
