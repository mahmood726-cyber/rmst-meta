"""R parity gate — compare RMSTmeta pooling against metafor::rma().

Covers two built-in examples hardcoded in rmst-meta.html:
  - NSCLC immunotherapy (5 trials from CheckMate-017/057, KEYNOTE-010/024, OAK)
  - SGLT2i heart failure (3 trials: DAPA-HF, EMPEROR-Reduced, SOLOIST-WHF)

The SGLT2i case at k=3 is a stress test (HKSJ becomes t_{k-1}; PI requires
k ≥ 3). Skipped cleanly when R or scipy is unavailable.
"""
from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
from pathlib import Path

import pytest

try:
    from scipy.optimize import minimize_scalar
    from scipy.stats import norm
    from scipy.stats import t as student_t
except ImportError:
    minimize_scalar = None
    student_t = None
    norm = None

RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript") or shutil.which("Rscript.exe")
HERE = Path(__file__).resolve().parent
R_DIR = HERE.parent / "validation" / "r_parity"
R_SCRIPT = R_DIR / "validate.R"
DATASETS = {
    "nsclc": (R_DIR / "nsclc_data.csv", R_DIR / "nsclc_r_results.csv"),
    "sglt2": (R_DIR / "sglt2_data.csv", R_DIR / "sglt2_r_results.csv"),
}

TOL = 1e-4
TOL_REML_TAU2 = 1e-3  # REML tau² wobblier under iterative optimisation


pytestmark = pytest.mark.skipif(
    not RSCRIPT or minimize_scalar is None,
    reason="Rscript or scipy not available",
)


def _load(csv_path: Path):
    out = []
    with csv_path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            diff = float(row["rmst1"]) - float(row["rmst2"])
            se = math.hypot(float(row["se1"]), float(row["se2"]))
            out.append((diff, se * se))
    return out


def _dl(yi, vi):
    k = len(yi)
    wi = [1.0 / v for v in vi]
    sw = sum(wi)
    mu_fe = sum(w * y for w, y in zip(wi, yi)) / sw
    q = sum(w * (y - mu_fe) ** 2 for w, y in zip(wi, yi))
    sw2 = sum(w * w for w in wi)
    c = sw - sw2 / sw
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    wre = [1.0 / (v + tau2) for v in vi]
    swre = sum(wre)
    mu = sum(w * y for w, y in zip(wre, yi)) / swre
    se = math.sqrt(1.0 / swre)
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    return mu, se, tau2, i2, q, wre, swre


def _reml(yi, vi):
    def neg_ll(tau2: float) -> float:
        if tau2 < 0:
            return float("inf")
        w = [1.0 / (v + tau2) for v in vi]
        sw = sum(w)
        mu = sum(wj * yj for wj, yj in zip(w, yi)) / sw
        ll = -0.5 * (
            sum(math.log(v + tau2) for v in vi)
            + sum(wj * (yj - mu) ** 2 for wj, yj in zip(w, yi))
            + math.log(sw)
        )
        return -ll

    res = minimize_scalar(neg_ll, bounds=(0.0, 100.0), method="bounded",
                          options={"xatol": 1e-10})
    tau2 = max(0.0, res.x)
    wre = [1.0 / (v + tau2) for v in vi]
    swre = sum(wre)
    mu = sum(w * y for w, y in zip(wre, yi)) / swre
    se = math.sqrt(1.0 / swre)
    return mu, se, tau2, wre, swre


def _hksj(yi, wre, swre, mu):
    """Hartung-Knapp-Sidik-Jonkman standard error.

    SE_HKSJ = sqrt(Σ w_re_i (y_i - μ)² / ((k-1) Σ w_re_i)).
    metafor does NOT floor this at the regular RE SE — older convention did,
    but current metafor releases compute the raw HKSJ SE. Floor would cause
    spurious disagreement, so we compute the raw value.
    """
    k = len(yi)
    q_h = sum(w * (y - mu) ** 2 for w, y in zip(wre, yi))
    var_h = q_h / ((k - 1) * swre)
    return math.sqrt(var_h)


@pytest.fixture(scope="module")
def r_results():
    res = subprocess.run(
        [RSCRIPT, str(R_SCRIPT)],
        cwd=str(R_DIR),
        capture_output=True, text=True,
    )
    if res.returncode != 0:
        pytest.skip(f"Rscript failed: {res.stderr.strip()[:200]}")
    out = {}
    for name, (_, results_csv) in DATASETS.items():
        out[name] = {}
        with results_csv.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                method = row["method"]
                out[name][method] = {k: (float(v) if v not in ("", "NA") else None)
                                     for k, v in row.items() if k != "method"}
    return out


@pytest.fixture(scope="module", params=["nsclc", "sglt2"])
def dataset(request):
    data_csv, _ = DATASETS[request.param]
    return request.param, _load(data_csv)


def test_dl_est_se_tau2_Q(r_results, dataset):
    name, studies = dataset
    yi = [y for y, _ in studies]; vi = [v for _, v in studies]
    mu, se, tau2, _, q, _, _ = _dl(yi, vi)
    r = r_results[name]["DL"]
    assert mu == pytest.approx(r["est"], abs=TOL), f"[{name}] DL est"
    assert se == pytest.approx(r["se"], abs=TOL),  f"[{name}] DL se"
    assert tau2 == pytest.approx(r["tau2"], abs=TOL), f"[{name}] DL tau2"
    assert q == pytest.approx(r["Q"], abs=TOL), f"[{name}] DL Q"


def test_dl_i_squared(r_results, dataset):
    name, studies = dataset
    yi = [y for y, _ in studies]; vi = [v for _, v in studies]
    _, _, _, i2, _, _, _ = _dl(yi, vi)
    assert i2 == pytest.approx(r_results[name]["DL"]["I2"], abs=1.0), f"[{name}] I²"


def test_reml_est_se_tau2(r_results, dataset):
    name, studies = dataset
    yi = [y for y, _ in studies]; vi = [v for _, v in studies]
    mu, se, tau2, _, _ = _reml(yi, vi)
    r = r_results[name]["REML"]
    assert mu == pytest.approx(r["est"], abs=TOL), f"[{name}] REML est"
    assert se == pytest.approx(r["se"], abs=TOL), f"[{name}] REML se"
    assert tau2 == pytest.approx(r["tau2"], abs=TOL_REML_TAU2), f"[{name}] REML tau²"


def test_hksj_se_parity(r_results, dataset):
    """HKSJ SE against metafor knha."""
    name, studies = dataset
    yi = [y for y, _ in studies]; vi = [v for _, v in studies]
    _, _, _, wre, swre = _reml(yi, vi)
    mu_reml, _, _, _, _ = _reml(yi, vi)
    se_h = _hksj(yi, wre, swre, mu_reml)
    r_hksj = r_results[name]["REML"].get("hksj_se")
    if r_hksj is None:
        pytest.skip(f"[{name}] metafor HKSJ unavailable")
    assert se_h == pytest.approx(r_hksj, abs=TOL), f"[{name}] HKSJ se"


def test_prediction_interval_parity(r_results, dataset):
    """Prediction interval from REML fit, t_{k-2} distribution."""
    name, studies = dataset
    yi = [y for y, _ in studies]; vi = [v for _, v in studies]
    k = len(yi)
    if k < 3:
        pytest.skip(f"[{name}] PI undefined for k<3")
    mu, se, tau2, _, _ = _reml(yi, vi)
    # metafor: when τ² = 0, predict() returns PI = CI because there is no
    # between-study variance for a new study to inherit. When τ² > 0,
    # metafor uses t_{k-1} degrees of freedom (not t_{k-2} as in
    # Higgins-Thompson 2009 — the current metafor release matches the
    # IntHout et al. 2016 convention).
    # scipy's bounded optimiser returns ~1e-13 fp dust rather than exact
    # zero when the likelihood is maximised at the boundary; treat anything
    # below 1e-9 as zero to match metafor's τ²=0 branch.
    if tau2 < 1e-9:
        z = norm.ppf(0.975)
        half = z * se
    else:
        t_val = student_t.ppf(0.975, df=k - 1)
        half = t_val * math.sqrt(tau2 + se * se)
    pi_lb, pi_ub = mu - half, mu + half
    r = r_results[name]["REML"]
    if r.get("pi_lb") is None or r.get("pi_ub") is None:
        pytest.skip(f"[{name}] metafor PI unavailable")
    assert pi_lb == pytest.approx(r["pi_lb"], abs=TOL), f"[{name}] PI lower"
    assert pi_ub == pytest.approx(r["pi_ub"], abs=TOL), f"[{name}] PI upper"


def test_k_equals_2_dl_edge_case():
    """DL with k=2 must not divide-by-zero; tau² should clamp to 0 when Q<k-1."""
    yi = [0.5, 1.0]
    vi = [0.04, 0.05]
    mu, se, tau2, i2, q, _, _ = _dl(yi, vi)
    assert math.isfinite(mu) and math.isfinite(se)
    assert tau2 >= 0.0
    assert i2 >= 0.0


def test_rmst_difference_se_formula():
    """Guard: SE(RMST_diff) = sqrt(se1² + se2²) by independence."""
    y, v = _load(DATASETS["nsclc"][0])[0]
    assert y == pytest.approx(11.2 - 9.4, abs=1e-9)
    assert math.sqrt(v) == pytest.approx(math.hypot(0.6, 0.5), abs=1e-9)
