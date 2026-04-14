# R parity reference for RMSTmeta.
# Reads nsclc_data.csv, computes RMST-difference random-effects pooling via
# metafor::rma() under DL and REML, writes pooled summaries to
# nsclc_r_results.csv. Called from tests/test_r_parity.py.

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

data_path <- file.path(script_dir, "nsclc_data.csv")
out_path  <- file.path(script_dir, "nsclc_r_results.csv")

d <- read.csv(data_path, stringsAsFactors = FALSE)
d$diff    <- d$rmst1 - d$rmst2
d$se_diff <- sqrt(d$se1^2 + d$se2^2)
d$vi      <- d$se_diff^2

rows <- list()
for (method in c("DL", "REML")) {
  fit <- rma(yi = d$diff, vi = d$vi, method = method)
  rows[[method]] <- data.frame(
    method = method,
    k      = fit$k,
    est    = as.numeric(fit$b),
    se     = fit$se,
    ci_lb  = fit$ci.lb,
    ci_ub  = fit$ci.ub,
    tau2   = fit$tau2,
    I2     = fit$I2,
    Q      = fit$QE,
    pQ     = fit$QEp
  )
}

out <- do.call(rbind, rows)
write.csv(out, out_path, row.names = FALSE)
cat("wrote", out_path, "\n")
