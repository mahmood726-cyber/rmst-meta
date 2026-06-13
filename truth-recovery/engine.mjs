// engine.mjs -- pure RMST meta-analysis core EXTRACTED VERBATIM from
// rmst-meta.html (606-930): normalCDF/Quantile, tQuantile, dlMeta (DL+HKSJ+PI),
// chi2CDF, tCDF, remlMeta. SAME math the app ships.

function normalCDF(x) {
  if (x < -8) return 0;
  if (x > 8) return 1;
  var a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
  var p = 0.3275911;
  var sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.SQRT2;
  var t = 1.0 / (1.0 + p * x);
  var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return 0.5 * (1.0 + sign * y);
}

/* Normal quantile (rational approximation, Beasley-Springer-Moro) */
function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;
  // Rational approximation
  var a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  var b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  var c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00, 2.938163982698783e+00
  ];
  var d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];
  var pLow = 0.02425, pHigh = 1 - pLow;
  var q, r;
  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

/* t-distribution quantile — exact table for small df, Cornish-Fisher for larger */
function tQuantile(p, df) {
  if (df <= 0) return normalQuantile(p);
  if (df >= 300) return normalQuantile(p);
  // Exact t_{0.975,df} lookup for df 1..30 (scipy.stats.t.ppf(0.975, df))
  var T975 = [0,
    12.7062,4.3027,3.1824,2.7764,2.5706,2.4469,2.3646,2.3060,2.2622,2.2281,
    2.2010,2.1788,2.1604,2.1448,2.1314,2.1199,2.1098,2.1009,2.0930,2.0860,
    2.0796,2.0739,2.0687,2.0639,2.0595,2.0555,2.0518,2.0484,2.0452,2.0423
  ];
  // For p=0.975 or p=0.025, use exact table
  if (df <= 30 && df === Math.floor(df)) {
    var z = normalQuantile(p);
    // For 0.975 and 0.025 specifically
    if (Math.abs(p - 0.975) < 1e-9) return T975[df];
    if (Math.abs(p - 0.025) < 1e-9) return -T975[df];
    // For other p values with small df, use exact formulas for df 1-2 and improved expansion otherwise
    if (df === 1) return Math.tan(Math.PI * (p - 0.5));
    if (df === 2) {
      var alpha2 = 2 * Math.min(p, 1 - p);
      var val = Math.sqrt(2.0 / alpha2 - 2.0);
      return p < 0.5 ? -val : val;
    }
    // Cornish-Fisher 4-term expansion for other quantiles
    var g1 = (z * z * z + z) / (4 * df);
    var g2 = (5 * z * z * z * z * z + 16 * z * z * z + 3 * z) / (96 * df * df);
    var g3 = (3 * Math.pow(z,7) + 19 * Math.pow(z,5) + 17 * z * z * z - 15 * z) / (384 * df * df * df);
    return z + g1 + g2 + g3;
  }
  // df > 30 or non-integer: Cornish-Fisher expansion
  var z2 = normalQuantile(p);
  var g1b = (z2 * z2 * z2 + z2) / (4 * df);
  var g2b = (5 * z2 * z2 * z2 * z2 * z2 + 16 * z2 * z2 * z2 + 3 * z2) / (96 * df * df);
  return z2 + g1b + g2b;
}

/* ---------- DerSimonian-Laird meta-analysis ---------- */
function dlMeta(diffs, ses) {
  var k = diffs.length;
  if (k === 0) return null;
  if (k === 1) {
    var z975 = normalQuantile(0.975);
    return { k:1, est:diffs[0], se:ses[0], ci:[diffs[0]-z975*ses[0], diffs[0]+z975*ses[0]], tau2:0, I2:0, Q:0, pQ:1, pi:null, hksj:null, weights:[1], diffs:diffs, ses:ses };
  }
  var wi = ses.map(function(s) { return 1 / (s * s); });
  var sumW = wi.reduce(function(a, b) { return a + b; }, 0);
  var sumWd = 0;
  for (var i = 0; i < k; i++) sumWd += wi[i] * diffs[i];
  var muFE = sumWd / sumW;

  // Q statistic
  var Q = 0;
  for (var i2 = 0; i2 < k; i2++) Q += wi[i2] * (diffs[i2] - muFE) * (diffs[i2] - muFE);

  var sumW2 = wi.reduce(function(a, w) { return a + w * w; }, 0);
  var C = sumW - sumW2 / sumW;
  var tau2 = Math.max(0, (Q - (k - 1)) / C);

  // RE weights
  var wiStar = ses.map(function(s) { return 1 / (s * s + tau2); });
  var sumWStar = wiStar.reduce(function(a, b) { return a + b; }, 0);
  var sumWdStar = 0;
  for (var i3 = 0; i3 < k; i3++) sumWdStar += wiStar[i3] * diffs[i3];
  var muRE = sumWdStar / sumWStar;
  var seRE = Math.sqrt(1 / sumWStar);

  var z975 = normalQuantile(0.975);
  var ci = [muRE - z975 * seRE, muRE + z975 * seRE];

  // I-squared
  var I2 = k > 1 ? Math.max(0, (Q - (k - 1)) / Q) * 100 : 0;

  // Prediction interval (requires k >= 3)
  var pi = null;
  if (k >= 3) {
    var tVal = tQuantile(0.975, k - 2);
    var piHalf = tVal * Math.sqrt(tau2 + seRE * seRE);
    pi = [muRE - piHalf, muRE + piHalf];
  }

  // HKSJ adjustment
  var hksj = null;
  if (k >= 2) {
    var qHksj = 0;
    for (var ih = 0; ih < k; ih++) {
      qHksj += wiStar[ih] * (diffs[ih] - muRE) * (diffs[ih] - muRE);
    }
    var varHksj = qHksj / ((k - 1) * sumWStar);
    var seHksj = Math.sqrt(Math.max(varHksj, seRE * seRE));
    var tValHksj = tQuantile(0.975, k - 1);
    hksj = {
      est: muRE,
      se: seHksj,
      ci: [muRE - tValHksj * seHksj, muRE + tValHksj * seHksj]
    };
  }

  return {
    k: k,
    est: muRE,
    se: seRE,
    ci: ci,
    tau2: tau2,
    I2: I2,
    Q: Q,
    pQ: 1 - chi2CDF(Q, k - 1),
    pi: pi,
    hksj: hksj,
    weights: wiStar,
    diffs: diffs,
    ses: ses
  };
}

/* Chi-squared CDF via regularized gamma (Wilson-Hilferty) */
function chi2CDF(x, df) {
  if (df <= 0 || x <= 0) return 0;
  // Wilson-Hilferty approximation
  var z = Math.pow(x / df, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df));
  z = z / Math.sqrt(2.0 / (9.0 * df));
  return normalCDF(z);
}

/* t-distribution CDF via bisection on tQuantile */
function tCDF(x, df) {
  if (df <= 0) return normalCDF(x);
  if (df >= 300) return normalCDF(x);
  // Bisection: find p such that tQuantile(p, df) = x
  var lo = 1e-12, hi = 1 - 1e-12, mid;
  for (var i = 0; i < 40; i++) {
    mid = (lo + hi) / 2;
    if (tQuantile(mid, df) < x) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

/* ---------- REML meta-analysis (Fisher scoring, Viechtbauer 2005) ---------- */
function remlMeta(diffs, ses) {
  var k = diffs.length;
  if (k === 0) return null;
  if (k === 1) {
    var z975s = normalQuantile(0.975);
    return { k:1, est:diffs[0], se:ses[0], ci:[diffs[0]-z975s*ses[0], diffs[0]+z975s*ses[0]], tau2:0, I2:0, Q:0, pQ:1, pi:null, hksj:null, weights:[1/(ses[0]*ses[0])], diffs:diffs, ses:ses, converged:true, iterations:0 };
  }

  // Fixed-effect Q (invariant, same as DL)
  var vi = ses.map(function(s) { return s * s; });
  var wiFE = vi.map(function(v) { return 1 / v; });
  var sumWFE = wiFE.reduce(function(a, b) { return a + b; }, 0);
  var muFE = 0;
  for (var i = 0; i < k; i++) muFE += wiFE[i] * diffs[i];
  muFE /= sumWFE;
  var Q = 0;
  for (var iq = 0; iq < k; iq++) Q += wiFE[iq] * (diffs[iq] - muFE) * (diffs[iq] - muFE);

  // DL tau2 as warm start
  var sumW2FE = wiFE.reduce(function(a, w) { return a + w * w; }, 0);
  var C = sumWFE - sumW2FE / sumWFE;
  var tau2 = Math.max(0, (Q - (k - 1)) / C);

  // Fisher scoring iteration
  var converged = false;
  var maxIter = 100;
  var iter;
  for (iter = 0; iter < maxIter; iter++) {
    var wiStar = vi.map(function(v) { return 1 / (v + tau2); });
    var sumW = wiStar.reduce(function(a, b) { return a + b; }, 0);
    var sumW2 = 0, sumW3 = 0;
    for (var iw = 0; iw < k; iw++) {
      sumW2 += wiStar[iw] * wiStar[iw];
      sumW3 += wiStar[iw] * wiStar[iw] * wiStar[iw];
    }
    var muStar = 0;
    for (var im = 0; im < k; im++) muStar += wiStar[im] * diffs[im];
    muStar /= sumW;

    var P = 0;
    for (var ip = 0; ip < k; ip++) P += wiStar[ip] * wiStar[ip] * (diffs[ip] - muStar) * (diffs[ip] - muStar);

    // REML score: proportional to (P - A) where A = tr(W) - tr(W^2)/tr(W)
    var A = sumW - sumW2 / sumW;
    var score = P - A;

    // Fisher information (expected): sum(wi*^2) - 2*sum(wi*^3)/sum(wi*) + (sum(wi*^2)/sum(wi*))^2
    var info = sumW2 - 2 * sumW3 / sumW + (sumW2 / sumW) * (sumW2 / sumW);

    if (info <= 0) break; // degenerate, stop

    var tau2New = Math.max(0, tau2 + score / info);
    if (Math.abs(tau2New - tau2) < 1e-8) {
      tau2 = tau2New;
      converged = true;
      break;
    }
    tau2 = tau2New;
  }

  // Self-check: score residual at convergence (only for interior solutions, not boundary tau2=0)
  if (converged && tau2 > 0) {
    var wiCheck = vi.map(function(v) { return 1 / (v + tau2); });
    var sumWCheck = wiCheck.reduce(function(a, b) { return a + b; }, 0);
    var muCheck = 0;
    for (var ic = 0; ic < k; ic++) muCheck += wiCheck[ic] * diffs[ic];
    muCheck /= sumWCheck;
    var qRE = 0;
    for (var iqc = 0; iqc < k; iqc++) qRE += wiCheck[iqc] * (diffs[iqc] - muCheck) * (diffs[iqc] - muCheck);
    var Acheck = sumWCheck - wiCheck.reduce(function(a, w) { return a + w * w; }, 0) / sumWCheck;
    if (Math.abs(qRE - Acheck) > 1e-4) {
      console.warn('REML score residual:', Math.abs(qRE - Acheck));
    }
  }

  // Compute final estimates with converged tau2
  var wiRE = vi.map(function(v) { return 1 / (v + tau2); });
  var sumWRE = wiRE.reduce(function(a, b) { return a + b; }, 0);
  var muRE = 0;
  for (var ir = 0; ir < k; ir++) muRE += wiRE[ir] * diffs[ir];
  muRE /= sumWRE;
  var seRE = Math.sqrt(1 / sumWRE);

  var z975 = normalQuantile(0.975);
  var ci = [muRE - z975 * seRE, muRE + z975 * seRE];
  var I2 = k > 1 ? Math.max(0, (Q - (k - 1)) / Q) * 100 : 0;

  // Prediction interval (k >= 3)
  var pi = null;
  if (k >= 3) {
    var tVal = tQuantile(0.975, k - 2);
    var piHalf = tVal * Math.sqrt(tau2 + seRE * seRE);
    pi = [muRE - piHalf, muRE + piHalf];
  }

  // HKSJ adjustment (k >= 2)
  var hksj = null;
  if (k >= 2) {
    var qHksj = 0;
    for (var ih = 0; ih < k; ih++) {
      qHksj += wiRE[ih] * (diffs[ih] - muRE) * (diffs[ih] - muRE);
    }
    var varHksj = qHksj / ((k - 1) * sumWRE);
    var seHksj = Math.sqrt(Math.max(varHksj, seRE * seRE));
    var tValHksj = tQuantile(0.975, k - 1);
    hksj = {
      est: muRE,
      se: seHksj,
      ci: [muRE - tValHksj * seHksj, muRE + tValHksj * seHksj]
    };
  }

  return {
    k: k,
    est: muRE,
    se: seRE,
    ci: ci,
    tau2: tau2,
    I2: I2,
    Q: Q,
    pQ: 1 - chi2CDF(Q, k - 1),
    pi: pi,
    hksj: hksj,
    weights: wiRE,
    diffs: diffs,
    ses: ses,
    converged: converged,
    iterations: iter + 1
  };
}

/* ---------- Egger's test for funnel plot asymmetry (OLS, Egger 1997) ---------- */

export { dlMeta, remlMeta, tQuantile, normalQuantile, normalCDF };
