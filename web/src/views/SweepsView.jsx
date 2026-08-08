import { useMemo, useState } from "react";
import {
  Area,
  CartesianGrid,
  Line,
  LineChart,
  ReferenceArea,
  ReferenceLine,
  ResponsiveContainer,
  XAxis,
  YAxis,
} from "recharts";

const STORM_COLOR = "#d1495b";
const QUIET_COLOR = "#4fa3ff";

// Reshapes sweeps.json's "storm|1.00"-keyed survival dicts into recharts
// rows, and folds the x1.00/x1.19 pair into a two-value range per condition
// so <Area> can shade the density-uncertainty band between them. No survival
// fraction is computed here -- both endpoints of every band are already in
// the JSON; this only pairs them up.
function buildRows(xKey, x, survival) {
  return x.map((xv, i) => {
    const s100 = survival[`storm|1.00`]?.[i];
    const s119 = survival[`storm|1.19`]?.[i];
    const q100 = survival[`quiet|1.00`]?.[i];
    const q119 = survival[`quiet|1.19`]?.[i];
    return {
      [xKey]: xv,
      storm_1: s100,
      storm_19: s119,
      storm_band: s100 != null && s119 != null ? [Math.min(s100, s119), Math.max(s100, s119)] : null,
      quiet_1: q100,
      quiet_19: q119,
      quiet_band: q100 != null && q119 != null ? [Math.min(q100, q119), Math.max(q100, q119)] : null,
    };
  });
}

export default function SweepsView({ sweeps }) {
  const [tab, setTab] = useState("altitude");

  return (
    <>
      <div className="view-header">
        <span className="view-title">Sweeps</span>
        <span className="view-note">
          Shaded band spans density scale 1.00 (NRLMSIS as reported) to 1.19
          (README.md's central finding). Solid = 1.00, dashed = 1.19.{" "}
          <strong className="band-caveat">Bounds, not scenarios</strong> — the
          two edges are the same systematic bias assumed absent and assumed
          present, so the truth lies between them. Neither edge is a forecast,
          and neither is a best or worst case.
        </span>
      </div>

      <div className="nav-tabs" style={{ alignSelf: "flex-start" }}>
        {[
          ["altitude", "Insertion altitude"],
          ["ram_area", "Ram area"],
          ["safe_mode", "Safe-mode timing"],
        ].map(([id, label]) => (
          <button
            key={id}
            className={`nav-tab${tab === id ? " active" : ""}`}
            onClick={() => setTab(id)}
          >
            {label}
          </button>
        ))}
      </div>

      {tab === "altitude" && <InsertionAltitudeSweep sweeps={sweeps} />}
      {tab === "ram_area" && <RamAreaSweep sweeps={sweeps} />}
      {tab === "safe_mode" && <SafeModeSweep sweeps={sweeps} />}
    </>
  );
}

function InsertionAltitudeSweep({ sweeps }) {
  const blk = sweeps.sweep_insertion_altitude;
  const rows = useMemo(() => buildRows("x_km", blk.x_km, blk.survival), [blk]);
  const [stormLo, stormHi] = sweeps.critical_altitude_km.storm;
  const [quietLo, quietHi] = sweeps.critical_altitude_km.quiet;

  return (
    <div className="panel">
      <div className="panel-title">Survival vs insertion altitude</div>
      <ResponsiveContainer width="100%" height={420}>
        <LineChart data={rows} margin={{ top: 8, right: 16, left: 0, bottom: 8 }}>
          <CartesianGrid strokeDasharray="2 4" stroke="#262b36" />
          <XAxis dataKey="x_km" type="number" domain={[175, 260]} allowDataOverflow stroke="#7a8194" fontSize={11}
            label={{ value: "insertion altitude (km)", position: "insideBottom", offset: -4, fill: "#7a8194", fontSize: 11 }} />
          <YAxis domain={[0, 1]} stroke="#7a8194" fontSize={11}
            label={{ value: "survival fraction", angle: -90, position: "insideLeft", fill: "#7a8194", fontSize: 11 }} />
          <ReferenceArea x1={stormLo} x2={stormHi} fill={STORM_COLOR} fillOpacity={0.12} />
          <ReferenceArea x1={quietLo} x2={quietHi} fill={QUIET_COLOR} fillOpacity={0.12} />
          <ReferenceLine x={sweeps.meta.actual_insertion_km} stroke="#f2f4f8" strokeDasharray="3 3"
            label={{ value: "actual insertion", position: "top", fill: "#f2f4f8", fontSize: 10 }} />
          <Area dataKey="storm_band" stroke="none" fill={STORM_COLOR} fillOpacity={0.18} isAnimationActive={false} />
          <Area dataKey="quiet_band" stroke="none" fill={QUIET_COLOR} fillOpacity={0.18} isAnimationActive={false} />
          <Line dataKey="storm_1" stroke={STORM_COLOR} strokeWidth={2} dot={false} isAnimationActive={false} name="storm, x1.00" />
          <Line dataKey="storm_19" stroke={STORM_COLOR} strokeWidth={1.25} strokeDasharray="4 3" dot={false} isAnimationActive={false} name="storm, x1.19" />
          <Line dataKey="quiet_1" stroke={QUIET_COLOR} strokeWidth={2} dot={false} isAnimationActive={false} name="quiet, x1.00" />
          <Line dataKey="quiet_19" stroke={QUIET_COLOR} strokeWidth={1.25} strokeDasharray="4 3" dot={false} isAnimationActive={false} name="quiet, x1.19" />
        </LineChart>
      </ResponsiveContainer>
      <SweepLegend />
      <p className="view-note">
        Analytic critical altitude (independent of this Monte Carlo sweep):
        storm {stormLo.toFixed(1)}–{stormHi.toFixed(1)} km, quiet {quietLo.toFixed(1)}–{quietHi.toFixed(1)} km.
      </p>
    </div>
  );
}

function RamAreaSweep({ sweeps }) {
  const blk = sweeps.sweep_ram_area;
  const rows = useMemo(() => buildRows("x_m2", blk.x_m2, blk.survival), [blk]);
  const [pubLo, pubHi] = blk.published_range_m2;

  return (
    <div className="panel">
      <div className="panel-title">Survival vs ram area (210 km insertion)</div>
      <ResponsiveContainer width="100%" height={420}>
        <LineChart data={rows} margin={{ top: 8, right: 16, left: 0, bottom: 8 }}>
          <CartesianGrid strokeDasharray="2 4" stroke="#262b36" />
          <XAxis dataKey="x_m2" type="number" stroke="#7a8194" fontSize={11}
            label={{ value: "nominal ram area (m²)", position: "insideBottom", offset: -4, fill: "#7a8194", fontSize: 11 }} />
          <YAxis domain={[0, 1]} stroke="#7a8194" fontSize={11}
            label={{ value: "survival fraction", angle: -90, position: "insideLeft", fill: "#7a8194", fontSize: 11 }} />
          <ReferenceArea x1={pubLo} x2={pubHi} fill="#8d99ae" fillOpacity={0.14}
            label={{ value: "Baruah range", position: "insideTop", fill: "#7a8194", fontSize: 10 }} />
          {blk.source_ranges_m2.map(([[a, b], color, label], i) => (
            <ReferenceArea key={i} x1={a} x2={b} fill={color} fillOpacity={0.14}
              label={{ value: label, position: "insideTop", fill: color, fontSize: 10 }} />
          ))}
          <Area dataKey="storm_band" stroke="none" fill={STORM_COLOR} fillOpacity={0.18} isAnimationActive={false} />
          <Area dataKey="quiet_band" stroke="none" fill={QUIET_COLOR} fillOpacity={0.18} isAnimationActive={false} />
          <Line dataKey="storm_1" stroke={STORM_COLOR} strokeWidth={2} dot={false} isAnimationActive={false} />
          <Line dataKey="storm_19" stroke={STORM_COLOR} strokeWidth={1.25} strokeDasharray="4 3" dot={false} isAnimationActive={false} />
          <Line dataKey="quiet_1" stroke={QUIET_COLOR} strokeWidth={2} dot={false} isAnimationActive={false} />
          <Line dataKey="quiet_19" stroke={QUIET_COLOR} strokeWidth={1.25} strokeDasharray="4 3" dot={false} isAnimationActive={false} />
        </LineChart>
      </ResponsiveContainer>
      <SweepLegend />
    </div>
  );
}

function SafeModeSweep({ sweeps }) {
  const blk = sweeps.sweep_safe_mode_timing;
  const rows = useMemo(() => buildRows("x_hours", blk.x_hours, blk.survival), [blk]);

  return (
    <div className="panel">
      <div className="panel-title">Survival vs safe-mode exit time (210 km)</div>
      <ResponsiveContainer width="100%" height={420}>
        <LineChart data={rows} margin={{ top: 8, right: 16, left: 0, bottom: 8 }}>
          <CartesianGrid strokeDasharray="2 4" stroke="#262b36" />
          <XAxis dataKey="x_hours" type="number" stroke="#7a8194" fontSize={11}
            label={{ value: "safe-mode exit time (hours)", position: "insideBottom", offset: -4, fill: "#7a8194", fontSize: 11 }} />
          <YAxis domain={[0, 1]} stroke="#7a8194" fontSize={11}
            label={{ value: "survival fraction", angle: -90, position: "insideLeft", fill: "#7a8194", fontSize: 11 }} />
          <Area dataKey="storm_band" stroke="none" fill={STORM_COLOR} fillOpacity={0.18} isAnimationActive={false} />
          <Area dataKey="quiet_band" stroke="none" fill={QUIET_COLOR} fillOpacity={0.18} isAnimationActive={false} />
          <Line dataKey="storm_1" stroke={STORM_COLOR} strokeWidth={2} dot={false} isAnimationActive={false} />
          <Line dataKey="storm_19" stroke={STORM_COLOR} strokeWidth={1.25} strokeDasharray="4 3" dot={false} isAnimationActive={false} />
          <Line dataKey="quiet_1" stroke={QUIET_COLOR} strokeWidth={2} dot={false} isAnimationActive={false} />
          <Line dataKey="quiet_19" stroke={QUIET_COLOR} strokeWidth={1.25} strokeDasharray="4 3" dot={false} isAnimationActive={false} />
        </LineChart>
      </ResponsiveContainer>
      <SweepLegend />
    </div>
  );
}

// Every chart in this view carries a shaded band, so the legend is where the
// "bounds, not scenarios" reading belongs -- it sits next to the swatch rather
// than only in the header, and follows the band into every tab.
function SweepLegend() {
  return (
    <>
      <div className="legend">
        <span className="legend-item"><span className="legend-swatch" style={{ background: STORM_COLOR }} /> Storm</span>
        <span className="legend-item"><span className="legend-swatch" style={{ background: QUIET_COLOR }} /> Quiet</span>
        <span className="legend-item">— solid: density ×1.00</span>
        <span className="legend-item">- - dashed: density ×1.19</span>
      </div>
      <p className="band-caveat-note">
        Shading shows <strong>bounds, not scenarios.</strong> The band's edges
        are the same density assumption applied and withheld — the answer lies
        between them. An edge is not a forecast.
      </p>
    </>
  );
}
