# Validate RMST Meta v1.1 against metafor
library(meta)
library(metafor)

# NSCLC data (two-arm -> differences)
rmst_exp  <- c(11.2, 12.3, 13.4, 12.8, 15.1)
se_exp    <- c(0.6, 0.5, 0.7, 0.4, 0.8)
rmst_ctrl <- c(9.4, 9.6, 10.1, 10.5, 11.3)
se_ctrl   <- c(0.5, 0.4, 0.6, 0.4, 0.7)

rmst_diff <- rmst_exp - rmst_ctrl
se_diff   <- sqrt(se_exp^2 + se_ctrl^2)

cat("=== NSCLC Dataset ===\n")
cat("Diffs:", rmst_diff, "\n")
cat("SEs:", round(se_diff, 4), "\n\n")

# DL
m_dl <- rma(yi = rmst_diff, sei = se_diff, method = "DL")
cat("DL tau2:", m_dl$tau2, "\n")
cat("DL est:", m_dl$beta[1], "CI:", m_dl$ci.lb, m_dl$ci.ub, "\n\n")

# REML
m_reml <- rma(yi = rmst_diff, sei = se_diff, method = "REML")
cat("REML tau2:", m_reml$tau2, "\n")
cat("REML est:", m_reml$beta[1], "CI:", m_reml$ci.lb, m_reml$ci.ub, "\n\n")

# Egger's test (OLS)
rt <- regtest(m_reml, model = "lm", predictor = "sei")
cat("Egger intercept:", coef(rt$fit)[1], "\n")
cat("Egger t:", rt$zval, "p:", rt$pval, "\n\n")

# SGLT2 data
rmst_exp2  <- c(20.1, 14.2, 8.1)
se_exp2    <- c(0.3, 0.2, 0.3)
rmst_ctrl2 <- c(19.4, 13.6, 7.5)
se_ctrl2   <- c(0.3, 0.2, 0.3)

rmst_diff2 <- rmst_exp2 - rmst_ctrl2
se_diff2   <- sqrt(se_exp2^2 + se_ctrl2^2)

cat("=== SGLT2 Dataset ===\n")
m_reml2 <- rma(yi = rmst_diff2, sei = se_diff2, method = "REML")
cat("REML tau2:", m_reml2$tau2, "\n")
cat("REML est:", m_reml2$beta[1], "CI:", m_reml2$ci.lb, m_reml2$ci.ub, "\n")

m_dl2 <- rma(yi = rmst_diff2, sei = se_diff2, method = "DL")
cat("DL tau2:", m_dl2$tau2, "\n")
cat("DL est:", m_dl2$beta[1], "CI:", m_dl2$ci.lb, m_dl2$ci.ub, "\n")
