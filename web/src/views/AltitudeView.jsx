import { useMemo } from "react";
import {
  CartesianGrid,
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  XAxis,
  YAxis,
} from "recharts";
import OutcomeLegend from "../components/OutcomeLegend";
import { usePlayback } from "../context/PlaybackContext";
import { buildAltitudeRows, withCriticalAltitudeSeconds } from "../lib/chartData";
import { formatCd } from "../lib/format";
import { outcomeColor } from "../lib/outcomeColors";

// Both Cd runs, always shown side by side -- the hard rule is that batch.json
// carries both and the UI must never default to one silently. There is
// nothing to toggle here: both panels render every time.
export default function AltitudeView({ batch }) {
  const { currentT } = usePlayback();

  return (
    <>
      <div className="view-header">
        <span className="view-title">Altitude</span>
        <span className="view-note">
          All {batch.meta.n_satellites} trajectories per run, both Cd
          conventions. Dashed line is the recovery boundary (critical
          altitude): the altitude at which the rated thrust would balance
          drag if orbit-raising resumed -- safe mode holds F = 0 throughout,
          so it is a counterfactual, not what actually happened.
        </span>
      </div>
      <div className="altitude-grid">
        {batch.runs.map((run) => (
          <RunPanel key={run.label} run={run} epochIso={batch.meta.epoch} currentT={currentT} />
        ))}
      </div>
    </>
  );
}

function RunPanel({ run, epochIso, currentT }) {
  const runWithSeconds = useMemo(
    () => withCriticalAltitudeSeconds(run, epochIso),
    [run, epochIso],
  );
  const rows = useMemo(() => buildAltitudeRows(runWithSeconds), [runWithSeconds]);
  const colorById = useMemo(() => {
    const m = new Map();
    for (const sat of run.satellites) m.set(sat.id, outcomeColor(sat.outcome));
    return m;
  }, [run]);

  return (
    <div className="panel">
      <div className="panel-title">
        Cd = {formatCd(run.config.cd)} — Cd·A ∈ [{run.config.effective_drag_range_m2[0].toFixed(2)},{" "}
        {run.config.effective_drag_range_m2[1].toFixed(2)}] m²
      </div>
      <ResponsiveContainer width="100%" height={340}>
        <LineChart data={rows} margin={{ top: 4, right: 12, left: 0, bottom: 4 }}>
          <CartesianGrid strokeDasharray="2 4" stroke="#262b36" />
          <XAxis
            dataKey="t_s"
            type="number"
            tickFormatter={(v) => (v / 3600).toFixed(0)}
            stroke="#7a8194"
            fontSize={11}
            label={{ value: "hours after epoch", position: "insideBottom", offset: -4, fill: "#7a8194", fontSize: 11 }}
          />
          <YAxis
            domain={[90, 220]}
            stroke="#7a8194"
            fontSize={11}
            label={{ value: "altitude (km)", angle: -90, position: "insideLeft", fill: "#7a8194", fontSize: 11 }}
          />
          <ReferenceLine y={100} stroke="#c8ccd6" strokeDasharray="2 2" />
          <ReferenceLine x={currentT} stroke="#4fa3ff" strokeWidth={1.5} />
          {run.satellites.map((sat) => (
            <Line
              key={sat.id}
              dataKey={`sat_${sat.id}`}
              stroke={colorById.get(sat.id)}
              dot={false}
              strokeWidth={1}
              isAnimationActive={false}
              connectNulls={false}
            />
          ))}
          <Line
            dataKey="h_crit"
            stroke="#f2f4f8"
            strokeDasharray="5 3"
            strokeWidth={1.75}
            dot={false}
            isAnimationActive={false}
            connectNulls
          />
        </LineChart>
      </ResponsiveContainer>
      <OutcomeLegend counts={run.outcome_counts} />
    </div>
  );
}
