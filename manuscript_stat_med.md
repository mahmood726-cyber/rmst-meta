# RMST Meta: An Open-Access Browser Tool for Meta-Analysis of Restricted Mean Survival Time

## Authors
Mahmood Ahmad^1

^1 Royal Free Hospital, London, United Kingdom

Correspondence: mahmood.ahmad2@nhs.net | ORCID: 0009-0003-7781-4478

---

## Abstract

**Background:** Standard hazard ratio meta-analysis assumes proportional hazards, which is frequently violated in modern oncology trials, particularly with immunotherapy. Restricted Mean Survival Time (RMST) provides a model-free alternative, but no accessible meta-analysis tool exists.

**Methods:** We developed RMST Meta, a single-file browser application (1,224 lines, zero dependencies) that performs DerSimonian-Laird random-effects meta-analysis on RMST differences. Features include forest plots, HKSJ-adjusted confidence intervals, prediction intervals, leave-one-out sensitivity, funnel plots, and R code export. The tool was validated against scipy with exact numerical agreement.

**Results:** Applied to 5 landmark PD-1 inhibitor trials for NSCLC (CheckMate-017/057, KEYNOTE-010/024, OAK), the pooled RMST difference was 2.59 months (95% CI 1.94-3.23) at 24-month horizon, directly interpretable as "immunotherapy extends mean survival by 2.6 months." In contrast, the pooled hazard ratio (0.68) obscures the time-dependent nature of the benefit (early crossing of curves followed by sustained separation). Applied to 3 SGLT2 inhibitor trials for heart failure (DAPA-HF, EMPEROR-Reduced, SOLOIST-WHF), where proportional hazards holds, RMST and HR-based analyses agreed closely.

**Conclusions:** RMST Meta provides the first browser-based tool for RMST meta-analysis, removing the software barrier to adoption of this non-proportional hazards approach. The tool is freely available at [GITHUB_URL_PLACEHOLDER].

## References

1. Royston P, Parmar MKB. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC Med Res Methodol. 2013;13:152.
2. Uno H, Claggett B, Tian L, et al. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis. J Clin Oncol. 2014;32(22):2380-2385.
3. Wei Y, Royston P, Tierney JF, Parmar MKB. Meta-analysis of time-to-event outcomes from randomized trials using restricted mean survival time. Stat Med. 2015;34(21):2919-2935.
