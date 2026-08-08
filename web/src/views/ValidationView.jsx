import { formatCd } from "../lib/format";

// The Baruah comparison table, discrepancies shown plainly. Every number
// here is read directly from sweeps.json's "validation" key or batch.json's
// "observed"/"runs" -- nothing is recomputed in the browser.
export default function ValidationView({ sweeps, batch }) {
  const v = sweeps.validation;

  return (
    <>
      <div className="view-header">
        <span className="view-title">Validation</span>
        <span className="view-note">
          PHYSICS.md §8's four tests, plus the Swarm C secondary check and the
          fleet loss-count comparison. Nothing tuned to match; discrepancies
          are reported as computed.
        </span>
      </div>

      <div className="panel">
        <div className="panel-title">Tests 1–3 (analytic)</div>
        <table className="data-table">
          <thead>
            <tr><th>Test</th><th>Measured</th><th>Limit</th><th>Result</th></tr>
          </thead>
          <tbody>
            <AnalyticRow t={v.analytic_tests.test_1_energy_conservation}
              name="1 — Energy conservation" valueKey="relative_change_in_a" />
            <AnalyticRow t={v.analytic_tests.test_2_thrust_spiral}
              name="2 — Thrust spiral" valueKey="relative_error" />
            <AnalyticRow t={v.analytic_tests.test_3_critical_density_fixed_point}
              name="3 — Critical density fixed point" valueKey="relative_rate" />
          </tbody>
        </table>
      </div>

      <div className="panel">
        <div className="panel-title">Test 4 — Baruah et al. (2024) reproduction</div>
        <p className="view-note" style={{ marginBottom: 10 }}>
          Cd = {formatCd(v.test_4_baruah_reproduction.cd)}, mass = {v.test_4_baruah_reproduction.mass_kg} kg,
          insertion = {v.test_4_baruah_reproduction.insertion_altitude_km} km, thrust = 0 N (safe mode).
          Acceptance: {v.test_4_baruah_reproduction.cases[0].acceptance_pct}% ("roughly 20%", PHYSICS.md §8).
        </p>
        <table className="data-table">
          <thead>
            <tr>
              <th>Ram area</th><th>Target</th><th>Outcome</th>
              <th>Altitude at reference</th><th>Reentry error</th><th>Decay error</th><th>Result</th>
            </tr>
          </thead>
          <tbody>
            {v.test_4_baruah_reproduction.cases.map((c) => {
              const errPct = c.reentry_error_pct ?? c.decay_error_pct;
              const pass = Math.abs(errPct) < c.acceptance_pct;
              return (
                <tr key={c.ram_area_m2}>
                  <td className="label">{c.ram_area_m2.toFixed(2)} m²</td>
                  <td>{c.target_altitude_km.toFixed(2)} km</td>
                  <td>{c.outcome}</td>
                  <td>{c.altitude_at_reference_km != null ? `${c.altitude_at_reference_km.toFixed(2)} km` : "—"}</td>
                  <td>{c.reentry_error_pct != null ? `${c.reentry_error_pct.toFixed(1)}%` : "—"}</td>
                  <td style={c.reentry_error_pct != null ? { color: "var(--text-dim)" } : undefined}>
                    {c.decay_error_pct != null ? `${c.decay_error_pct.toFixed(1)}%` : "—"}
                    {c.reentry_error_pct != null && (
                      <span style={{ fontSize: 10 }}> (not the tested metric here — see README.md)</span>
                    )}
                  </td>
                  <td className={pass ? "pass" : "fail"}>{pass ? "PASS" : "FAIL"}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
        <p className="view-note" style={{ marginTop: 10 }}>
          Implied uniform density multiplier needed to match each case exactly:{" "}
          {Object.entries(v.test_4_baruah_reproduction.implied_density_multiplier)
            .map(([area, k]) => `${area} m² → ×${k.toFixed(3)}`)
            .join("  ·  ")}
          . Two cases 4.48× apart in drag area agree to within a few percent — see README.md's central finding.
        </p>
      </div>

      <div className="panel">
        <div className="panel-title">
          Fleet loss count (2022-02-03 → {batch.meta.window_end.slice(0, 10)})
        </div>
        <table className="data-table">
          <thead>
            <tr><th>Source</th><th>Lost</th><th>Survived</th></tr>
          </thead>
          <tbody>
            <tr>
              <td className="label">Observed ({batch.observed.source})</td>
              <td>{batch.observed.lost}</td>
              <td>{batch.observed.survived}</td>
            </tr>
            {batch.runs.map((r) => {
              const lost = r.outcome_counts.REENTERED ?? 0;
              return (
                <tr key={r.label}>
                  <td className="label">Simulated, Cd = {formatCd(r.config.cd)} ({r.description})</td>
                  <td>{lost}</td>
                  <td>{batch.meta.n_satellites - lost}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <div className="panel">
        <div className="panel-title">
          Swarm C (secondary) <span className="badge badge-weakest">flagged weakest</span>
        </div>
        <table className="data-table">
          <thead>
            <tr><th>Quantity</th><th>Value</th></tr>
          </thead>
          <tbody>
            <tr><td className="label">Altitude / mass / ram area / Cd</td>
              <td>{v.swarm_c_secondary.altitude_km} km, {v.swarm_c_secondary.mass_kg} kg,{" "}
                {v.swarm_c_secondary.ram_area_m2} m², Cd {formatCd(v.swarm_c_secondary.cd)}</td></tr>
            <tr><td className="label">CastOrbit decay</td><td>{v.swarm_c_secondary.decay_m.toFixed(2)} m</td></tr>
            <tr><td className="label">Paper's modelled decay</td>
              <td>{v.swarm_c_secondary.paper_modelled_decay_m.toFixed(2)} m
                (×{v.swarm_c_secondary.implied_multiplier_vs_paper.toFixed(3)})</td></tr>
            <tr><td className="label">Observed decay</td>
              <td>{v.swarm_c_secondary.observed_decay_m.toFixed(2)} m
                (×{v.swarm_c_secondary.implied_multiplier_vs_observed.toFixed(3)})</td></tr>
          </tbody>
        </table>
        <p className="view-note" style={{ marginTop: 10 }}>{v.swarm_c_secondary.flagged_reason}</p>
      </div>
    </>
  );
}

function AnalyticRow({ t, name, valueKey }) {
  return (
    <tr>
      <td className="label">{name}</td>
      <td>{t[valueKey].toExponential(2)}</td>
      <td>{t.limit.toExponential(0)}</td>
      <td className={t.passed ? "pass" : "fail"}>{t.passed ? "PASS" : "FAIL"}</td>
    </tr>
  );
}
