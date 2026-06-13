# Truth-recovery yardstick — rmst-meta

**Verdict: HKSJ recommendation (measured) + an instructive HONEST NEGATIVE on the
prediction-interval df. The primary z-interval under-covers the true pooled RMST
difference under heterogeneity (HKSJ fixes it); but the "standard" t_{k-1} PER df
would NOT improve the prediction interval — the shipped (superseded) t_{k-2} is
marginally better for the new-study target.**

## Method
rmst-meta pools RMST differences with a DL/REML random-effects model, a primary
z-interval, an HKSJ-modified CI card, and a prediction interval. Using the app's
OWN `dlMeta` (engine.mjs, verbatim) on a known `(mu=0.3, tau2=0.05)` + selection
DGP (step / Copas), we measure coverage of the true pooled effect and of a NEW
study's TRUE effect. k∈{5,10,20}, 3000 reps/cell.

## Results

### CI coverage of the true pooled RMST difference
| scenario     | k  | z-interval (primary) | HKSJ card |
|--------------|----|---------------------:|----------:|
| none         | 5  | 0.895 | 0.961 |
| step_weak    | 20 | 0.820 | 0.866 |
| copas_weak   | 10 | 0.869 | 0.920 |
| step_strong  | 10 | 0.315 | 0.459 |

### PI coverage of a NEW study's true effect
| scenario   | k  | t_{k-2} (shipped) | t_{k-1} (Cochrane "standard") |
|------------|----|------------------:|------------------------------:|
| none       | 5  | 0.912 | 0.889 |
| step_weak  | 5  | 0.894 | 0.868 |
| copas_weak | 5  | 0.894 | 0.865 |
| none       | 10 | 0.869 | 0.866 |

## Findings (all measured)
1. **The primary displayed CI (z-interval) under-covers under heterogeneity.** It
   drops to 0.82–0.87 under weak selection and is already ~0.90 clean at k=5. The
   **HKSJ card recovers +4–9pp** toward nominal (0.96 clean, 0.87–0.95 under weak
   selection). → **recommend HKSJ as the primary CI**, not a secondary card
   (consistent with the proportionma / dta / htmlpairwise findings). The HKSJ here
   already uses the IQWiG variance floor `max(varHKSJ, seRE²)`, which is correct.
2. **HONEST NEGATIVE on the PI df.** The shipped PI uses `t_{k-2}`, which is
   *superseded* (Cochrane Handbook v6.5 and metafor v4 use `t_{k-1}`). The obvious
   "standardise it" change would be `t_{k-2} → t_{k-1}`. **Measured, that change
   does NOT help** — against a new study's true effect, the wider `t_{k-2}` covers
   *slightly better* in every cell (e.g. 0.912 vs 0.889 clean k=5); `t_{k-1}` is
   narrower and marginally worse. So do **not** make this change for coverage
   reasons. (It may still be worth aligning with Cochrane for *consistency*, but
   that is a convention choice, not a truth-recovery improvement.)
3. **The real PI gap is τ² estimation, not the df.** Both formulas under-cover the
   new-study target by ~5–8pp even with no selection, because DL `τ²` is downward
   biased at small k — neither df choice fixes that. (Same root cause as the
   wave-1 htmlpairwise PI finding.)
4. **Strong selection wrecks every IV interval** (z and HKSJ both collapse, down
   to 0.08–0.13 at step_strong k=20) — expected; no inverse-variance interval
   models selection (matches meta-stats-core / ubcma).

## What did NOT transfer
RMST differences are continuous effect sizes, so the wave-1 pairwise selection DGP
transferred directly; NPE/conformal machinery does not apply to this formula
library. No runtime dependency added; engine unchanged.

## Reproduce
```
node truth-recovery/harness.mjs --reps 3000
node --test truth-recovery/test-truth-recovery.mjs
```
