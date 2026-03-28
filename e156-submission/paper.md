Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

RMST Meta: A Browser Tool for Meta-Analysis of Restricted Mean Survival Time Beyond the Hazard Ratio

Can restricted mean survival time differences be pooled across trials in a browser-based meta-analysis tool that avoids the proportional hazards assumption inherently required by standard hazard ratio synthesis? The tool was applied to five landmark PD-1 inhibitor trials for advanced non-small-cell lung cancer totaling over 5000 participants at a 24-month horizon. DerSimonian-Laird random-effects pooling with optional Hartung-Knapp-Sidik-Jonkman adjustment and prediction intervals was implemented in 1224 lines of dependency-free JavaScript. The pooled RMST mean difference was 2.59 months favoring immunotherapy (95% CI 1.94 to 3.23), directly interpretable as additional survival time gained per patient. Leave-one-out sensitivity analysis and funnel plot inspection confirmed stability, with no single trial altering the conclusion and no evidence of asymmetry. RMST meta-analysis provides a clinically intuitive absolute time-scale summary that complements the relative hazard ratio in settings where treatment benefit evolves over time. The tool is limited to pre-computed RMST differences and cannot reconstruct individual patient survival curves from published Kaplan-Meier figures.

Outside Notes

Type: methods
Primary estimand: Pooled RMST difference (months)
App: RMST Meta v1.1
Data: 5 PD-1 inhibitor NSCLC trials (CheckMate-017/057, KEYNOTE-010/024, OAK)
Code: https://github.com/mahmood726-cyber/rmst-meta
Version: 1.1
Certainty: moderate
Validation: DRAFT

References

1. Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol. 2012;12:9.
2. Tierney JF, Stewart LA, Ghersi D, Burdett S, Sydes MR. Practical methods for incorporating summary time-to-event data into meta-analysis. Trials. 2007;8:16.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
