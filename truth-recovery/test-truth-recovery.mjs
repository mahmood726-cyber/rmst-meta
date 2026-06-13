// node --test truth-recovery/test-truth-recovery.mjs
// Measured invariants for the rmst-meta truth-recovery yardstick. Seeded; no
// hand-entered numbers.
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { generate, makeRng } from './dgp.mjs';
import { runCell, runGrid } from './harness.mjs';

describe('DGP', () => {
  it('is reproducible for a fixed seed', () => {
    const a = generate(0.3, 0.05, 10, 'step_weak', makeRng(7));
    const b = generate(0.3, 0.05, 10, 'step_weak', makeRng(7));
    assert.deepEqual(a.yi, b.yi);
  });
});

describe('Truth-recovery (measured)', () => {
  it('the primary z-interval under-covers the true RMST diff under heterogeneity; HKSJ recovers it', () => {
    const r = runCell(0.3, 0.05, 20, 'step_weak', 2500, makeRng(20260613), makeRng(99));
    assert.ok(r.zCI < 0.90, `z-CI coverage ${r.zCI} not under-covering`);
    assert.ok(r.hksjCI > r.zCI + 0.03, `HKSJ ${r.hksjCI} did not recover over z ${r.zCI}`);
  });

  it('HONEST NEGATIVE: switching the shipped t_{k-2} PI to the "standard" t_{k-1} does NOT improve new-study coverage', () => {
    // The app uses t_{k-2} (superseded by Cochrane v6.5 / metafor v4 = t_{k-1}).
    // But measured against a new study's TRUE effect, the wider t_{k-2} covers at
    // least as well -- both under-cover, and t_{k-1} (narrower) is no better.
    const grid = runGrid({ reps: 2000, scenarios: ['none', 'step_weak', 'copas_weak'] });
    let k2Better = 0, total = 0;
    for (const c of grid) {
      if (c.results.piK2 == null) continue;
      total++;
      if (c.results.piK2 >= c.results.piK1 - 0.005) k2Better++;
    }
    assert.ok(k2Better === total,
      `t_{k-1} beat shipped t_{k-2} in ${total - k2Better}/${total} cells -- the standard fix would help`);
  });

  it('clean PI under-covers a new study even with the shipped formula (the gap is tau2 estimation, not df)', () => {
    const r = runCell(0.3, 0.05, 10, 'none', 3000, makeRng(1), makeRng(2));
    assert.ok(r.piK2 < 0.93, `clean PI coverage ${r.piK2} -- expected mild under-coverage`);
  });
});
