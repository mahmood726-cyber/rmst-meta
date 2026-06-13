# RMST Meta

[![ci](https://github.com/mahmood726-cyber/rmst-meta/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/rmst-meta/actions/workflows/ci.yml) [![codeql](https://github.com/mahmood726-cyber/rmst-meta/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/rmst-meta/actions/workflows/codeql.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)

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
Rscript validate.R           # numerical parity against R metafor (requires R on PATH)
```

If `Rscript` is not on `PATH`, set `RSCRIPT=/full/path/to/Rscript` before running the parity gate.

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
