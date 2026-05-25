# RMST Meta

A browser-based tool for meta-analysis of restricted mean survival time (RMST) differences across trials, avoiding the proportional-hazards assumption inherently required by standard hazard-ratio synthesis.

**Live dashboard:** <https://mahmood726-cyber.github.io/rmstmeta/>

## What it does

- Pools pre-computed RMST differences across studies (DerSimonian-Laird random-effects).
- Optional Hartung-Knapp-Sidik-Jonkman adjustment and prediction intervals.
- Leave-one-out sensitivity, funnel plot, Egger's test.
- ~1200 lines of dependency-free JavaScript — works offline.

## Run

Open `index.html` (or `rmst-meta.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q          # Python smoke + parity tests
Rscript validate.R           # numerical parity against R survRM2 (requires R)
```

## Repo layout

| Path | Purpose |
|---|---|
| `rmst-meta.html` | the dashboard (main artifact) |
| `index.html` | landing page |
| `tests/` | pytest smoke + R-parity tests |
| `validate.R` | R reference implementation for parity |
| `manuscript_stat_med.md` | Stat Med submission manuscript |
| `e156-submission/` | E156 micro-paper bundle |
| `E156-PROTOCOL.md` | project metadata (E156 entry #148) |

## Citation

See `manuscript_stat_med.md` for full methodological reference. The tool implements the RMST framework of Uno et al. 2014 (`doi:10.1200/JCO.2014.55.2208`).

## License

See `LICENSE` (MIT).
