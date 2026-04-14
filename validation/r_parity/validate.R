# R parity reference for RMSTmeta.
# Reads nsclc_data.csv and sglt2_data.csv, computes RMST-difference
# random-effects pooling via metafor::rma() under DL and REML, plus HKSJ
# and prediction interval for the REML fit. Writes nsclc_r_results.csv
# and sglt2_r_results.csv consumed by tests/test_r_parity.py.

suppressPackageStartupMessages({
  library(metafor)
})

# Locate script dir via --file=... arg (set by Rscript), fall back to cwd.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_dir <- dirname(sub("^--file=", "", file_arg))
} else {
  script_dir <- getwd()
}

run_one <- function(data_csv, out_csv) {
  d <- read.csv(file.path(script_dir, data_csv), stringsAsFactors = FALSE)
  d$diff    <- d$rmst1 - d$rmst2
  d$se_diff <- sqrt(d$se1^2 + d$se2^2)
  d$vi      <- d$se_diff^2

  rows <- list()
  for (method in c("DL", "REML")) {
    fit <- rma(yi = d$diff, vi = d$vi, method = method)
    # HKSJ-adjusted fit (method="knha") for the same pool
    fit_knha <- tryCatch(
      rma(yi = d$diff, vi = d$vi, method = method, test = "knha"),
      error = function(e) NULL
    )
    # Prediction interval (REML-like; DL is a fallback)
    pi <- tryCatch(predict(fit, level = 95), error = function(e) NULL)

    row <- data.frame(
      method = method,
      k      = fit$k,
      est    = as.numeric(fit$b),
      se     = fit$se,
      ci_lb  = fit$ci.lb,
      ci_ub  = fit$ci.ub,
      tau2   = fit$tau2,
      I2     = fit$I2,
      Q      = fit$QE,
      pQ     = fit$QEp,
      hksj_se = if (!is.null(fit_knha)) fit_knha$se else NA_real_,
      hksj_ci_lb = if (!is.null(fit_knha)) fit_knha$ci.lb else NA_real_,
      hksj_ci_ub = if (!is.null(fit_knha)) fit_knha$ci.ub else NA_real_,
      pi_lb = if (!is.null(pi)) pi$pi.lb else NA_real_,
      pi_ub = if (!is.null(pi)) pi$pi.ub else NA_real_
    )
    rows[[method]] <- row
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(script_dir, out_csv), row.names = FALSE)
}

run_one("nsclc_data.csv", "nsclc_r_results.csv")
run_one("sglt2_data.csv", "sglt2_r_results.csv")
cat("wrote nsclc_r_results.csv and sglt2_r_results.csv\n")
