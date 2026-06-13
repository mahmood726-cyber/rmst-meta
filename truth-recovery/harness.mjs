// ============================================================
// harness.mjs -- Truth-recovery yardstick for rmst-meta.
//
// rmst-meta pools RMST differences with a DL/REML random-effects model, a primary
// z-interval CI, an HKSJ-modified CI, and a prediction interval. Using the app's
// OWN dlMeta (engine.mjs, verbatim) on a known (mu, tau^2) + selection DGP we
// measure two things against the truth:
//   (1) CI: does the primary z-interval cover the true pooled RMST difference,
//       and does the HKSJ card recover what it loses under heterogeneity?
//   (2) PI: the shipped prediction interval uses t_{k-2}. Cochrane v6.5 / metafor
//       v4 use t_{k-1}. We measure coverage of a NEW study's TRUE effect under
//       both, instead of assuming the "standard" df is better.
//
// Truth-first: every number printed comes from seeded simulation here.
// Run:  node truth-recovery/harness.mjs --reps 2000
// ============================================================

import { dlMeta, tQuantile } from './engine.mjs';
import { generate, makeRng, SCENARIOS } from './dgp.mjs';

const BASE_SEED = 20260613;
function randn(rng) { let u1 = rng(), u2 = rng(); if (u1 < 1e-12) u1 = 1e-12; return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2); }
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

export function runCell(mu, tau2, k, scenario, reps, rng, rngNew) {
  const acc = { 'z-CI': { c: 0 }, 'HKSJ-CI': { c: 0 }, 'PI t_{k-2} (shipped)': { c: 0 }, 'PI t_{k-1} (Cochrane)': { c: 0 } };
  let n = 0, biasSum = 0, nPI = 0;
  for (let r = 0; r < reps; r++) {
    const { yi, vi } = generate(mu, tau2, k, scenario, rng);
    const se = vi.map(Math.sqrt);
    let m; try { m = dlMeta(yi, se); } catch { continue; }
    if (!m || !isFinite(m.est)) continue;
    n++;
    biasSum += m.est - mu;
    if (m.ci[0] <= mu && mu <= m.ci[1]) acc['z-CI'].c++;
    if (m.hksj && m.hksj.ci[0] <= mu && mu <= m.hksj.ci[1]) acc['HKSJ-CI'].c++;
    // new study's TRUE effect ~ N(mu, tau2) (unconditional RE realisation)
    const thetaNew = mu + Math.sqrt(tau2) * randn(rngNew);
    if (m.pi && k >= 3) {
      nPI++;
      if (m.pi[0] <= thetaNew && thetaNew <= m.pi[1]) acc['PI t_{k-2} (shipped)'].c++;
      const half = tQuantile(0.975, k - 1) * Math.sqrt(m.tau2 + m.se * m.se);
      const lo = m.est - half, hi = m.est + half;
      if (lo <= thetaNew && thetaNew <= hi) acc['PI t_{k-1} (Cochrane)'].c++;
    }
  }
  return {
    n, bias: +(biasSum / n).toFixed(4),
    zCI: +(acc['z-CI'].c / n).toFixed(4),
    hksjCI: +(acc['HKSJ-CI'].c / n).toFixed(4),
    piK2: nPI ? +(acc['PI t_{k-2} (shipped)'].c / nPI).toFixed(4) : null,
    piK1: nPI ? +(acc['PI t_{k-1} (Cochrane)'].c / nPI).toFixed(4) : null,
  };
}

export function runGrid({ reps = 2000, ks = [5, 10, 20], mu = 0.3, tau2 = 0.05, scenarios = SCENARIOS } = {}) {
  const rng = makeRng(BASE_SEED), rngNew = makeRng(BASE_SEED ^ 0x9e3779b9);
  const grid = [];
  for (const scen of scenarios) for (const k of ks) grid.push({ scen, k, results: runCell(mu, tau2, k, scen, reps, rng, rngNew) });
  return grid;
}

const isMain = process.argv[1]?.endsWith('harness.mjs');
if (isMain) {
  const i = process.argv.indexOf('--reps');
  const reps = i >= 0 ? Number(process.argv[i + 1]) : 2000;
  const t0 = Date.now();
  const grid = runGrid({ reps });
  console.log(`\n# Truth-recovery yardstick -- rmst-meta`);
  console.log(`reps=${reps}/cell  mu=0.3 tau2=0.05  seed=${BASE_SEED}\n`);
  console.log('## CI coverage of true pooled RMST diff: primary z-interval vs HKSJ card\n');
  console.log('scenario       k    z-CI    HKSJ-CI');
  for (const c of grid) console.log(c.scen.padEnd(14), String(c.k).padStart(2), String(c.results.zCI).padStart(7), String(c.results.hksjCI).padStart(8));
  console.log('\n## PI coverage of a NEW study true effect: shipped t_{k-2} vs Cochrane t_{k-1}\n');
  console.log('scenario       k    t_{k-2}(shipped)   t_{k-1}(Cochrane)');
  for (const c of grid) console.log(c.scen.padEnd(14), String(c.k).padStart(2), String(c.results.piK2).padStart(13), String(c.results.piK1).padStart(17));
  console.log(`\n(${(Date.now() - t0) / 1000}s)`);
}
