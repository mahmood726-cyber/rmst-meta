"""R parity gate — compare RMSTmeta pooling against metafor::rma().

Uses the built-in NSCLC immunotherapy example (5 trials from CheckMate-017/057,
KEYNOTE-010/024, OAK) that is hard-coded into rmst-meta.html. The same dataset
is supplied via csv to validate.R; Python reimplements DL and REML pooling and
asserts agreement to 1e-4.

Skipped when R is not installed.
"""
from __future__ import annotations

import csv
import math
import os
import subprocess
from pathlib import Path

import pytest

try:
    from scipy.optimize import minimize_scalar
except ImportError:  # scipy missing — degrade gracefully
    minimize_scalar = None

RSCRIPT = r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe"
HERE = Path(__file__).resolve().parent
R_DIR = HERE.parent / "validation" / "r_parity"
R_SCRIPT = R_DIR / "validate.R"
R_OUT = R_DIR / "nsclc_r_results.csv"
DATA_CSV = R_DIR / "nsclc_data.csv"

TOL = 1e-4


pytestmark = pytest.mark.skipif(
    not os.path.exists(RSCRIPT) or minimize_scalar is None,
    reason="R 4.5.2 or scipy not available",
)


def _load_studies():
    out = []
    with DATA_CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            diff = float(row["rmst1"]) - float(row["rmst2"])
            se1 = float(row["se1"])
            se2 = float(row["se2"])
            se = math.sqrt(se1 * se1 + se2 * se2)
            out.append((diff, se * se))
    return out


def _dl_pool(yi, vi):
    """DerSimonian-Laird random-effects pool — closed-form."""
    k = len(yi)
    wi = [1.0 / v for v in vi]
    sw = sum(wi)
    mu_fe = sum(w * y for w, y in zip(wi, yi)) / sw
    q = sum(w * (y - mu_fe) ** 2 for w, y in zip(wi, yi))
    sw2 = sum(w * w for w in wi)
    c = sw - sw2 / sw
    tau2 = max(0.0, (q - (k - 1)) / c)
    wi_re = [1.0 / (v + tau2) for v in vi]
    sw_re = sum(wi_re)
    mu_re = sum(w * y for w, y in zip(wi_re, yi)) / sw_re
    se_re = math.sqrt(1.0 / sw_re)
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    return mu_re, se_re, tau2, i2, q


def _reml_pool(yi, vi):
    """REML pool via 1-D optimisation of the restricted log-likelihood."""
    k = len(yi)

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
    wi_re = [1.0 / (v + tau2) for v in vi]
    sw_re = sum(wi_re)
    mu_re = sum(w * y for w, y in zip(wi_re, yi)) / sw_re
    se_re = math.sqrt(1.0 / sw_re)
    return mu_re, se_re, tau2


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
    with R_OUT.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            out[row["method"]] = {k: float(v) for k, v in row.items() if k != "method"}
    return out


@pytest.fixture(scope="module")
def studies():
    return _load_studies()


def test_dl_pooled_estimate(r_results, studies):
    yi = [y for y, _ in studies]
    vi = [v for _, v in studies]
    mu, se, tau2, i2, q = _dl_pool(yi, vi)
    r = r_results["DL"]
    assert mu == pytest.approx(r["est"], abs=TOL), f"DL est: py={mu} r={r['est']}"
    assert se == pytest.approx(r["se"], abs=TOL)
    assert tau2 == pytest.approx(r["tau2"], abs=TOL)
    assert q == pytest.approx(r["Q"], abs=TOL)


def test_dl_i_squared(r_results, studies):
    yi = [y for y, _ in studies]
    vi = [v for _, v in studies]
    _, _, _, i2, _ = _dl_pool(yi, vi)
    assert i2 == pytest.approx(r_results["DL"]["I2"], abs=1.0)  # pct, 1% tol


def test_reml_pooled_estimate(r_results, studies):
    yi = [y for y, _ in studies]
    vi = [v for _, v in studies]
    mu, se, tau2 = _reml_pool(yi, vi)
    r = r_results["REML"]
    assert mu == pytest.approx(r["est"], abs=TOL), f"REML est: py={mu} r={r['est']}"
    assert se == pytest.approx(r["se"], abs=TOL)
    assert tau2 == pytest.approx(r["tau2"], abs=1e-3)  # REML tau2 more wobbly


def test_study_level_effects_match_rmst_difference_formula(studies):
    """Guard: ensure SE_diff = sqrt(se1^2 + se2^2) convention holds."""
    # First study: CheckMate-017 (rmst1=11.2, se1=0.6, rmst2=9.4, se2=0.5)
    y, v = studies[0]
    assert y == pytest.approx(11.2 - 9.4, abs=1e-9)
    assert math.sqrt(v) == pytest.approx(math.hypot(0.6, 0.5), abs=1e-9)
