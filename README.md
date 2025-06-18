
## Estimating and presenting nonlinear associations and interactions with Restricted Cubic Splines (RCS)

Last updated: 2025-06-18

------------------------------------------------------------------------

- [Slides](https://andreabellavia.github.io/RCSplines/continuous/2025_01_29_catalyst_Bellavia_.pdf)
  for the Harvard Catalyst seminar on January
  29, 2025

### 1. RCS to model a continuous variable in logistic and Cox regression

Supplementary material for the manuscript [*“Estimating and presenting
nonlinear associations with Restricted Cubic Splines”*](https://academic.oup.com/ije/article/54/4/dyaf088/8166023?searchresult=1)

- [Dataset](https://github.com/andreabellavia/RCSplines/blob/main/continuous/data_example.xlsx)
  used for the examples

#### R material

- [Code](https://github.com/andreabellavia/RCSplines/blob/main/continuous/R/rcs_logistic_cox.qmd)
  to reproduce the logistic and Cox regression examples

- [Output](https://andreabellavia.github.io/RCSplines/continuous/R/rcs_logistic_cox.html)
  for the examples in html format

- [Functions](https://github.com/andreabellavia/RCSplines/blob/main/continuous/R/R_functions.R)
  used for the examples

#### SAS material

- Macros to plot ORs from a logistic model: [unadjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/OR_spline.sas)
  / [adjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/OR_splines_adj.sas)

- Macros to plot HRs from a Cox model: [unadjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/HR_splines.sas)
  / [adjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/HR_splines_adj.sas)

- Macros to plot predicted probabilities after logistic model:
  [unadjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/logistic_splines_evprob.sas)
  / [adjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/logistic_splines_evprob_adj.sas)

- Macros to plot absolute risks after Cox model: [unadjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/cox_splines_eventprob.sas)
  / [adjusted
  models](https://github.com/andreabellavia/RCSplines/blob/main/continuous/SAS/cox_splines_eventprob_adj.sas)

#### Stata material

- [Code](https://github.com/andreabellavia/RCSplines/blob/main/continuous/Stata/rcs_logistic.do)
  to reproduce the logistic regression example

- [Code](https://github.com/andreabellavia/RCSplines/blob/main/continuous/Stata/rcs_cox.do)
  to reproduce the Cox regression example

### 2. RCS to model a nonlinear interaction between a binary and a continuous variable

#### R package `interactionRCS`

- Package on [CRAN](https://cran.r-project.org/package=interactionRCS)

- [Vignette](https://cran.r-project.org/web/packages/interactionRCS/vignettes/vignette.html)

#### Absolute Risks and Absolute Risk Differences from a Cox model with nonlinear interactions

- Supplemenraty material for the
  [paper](https://academic.oup.com/aje/article/193/8/1155/7678921)
  *“Estimating and presenting hazard ratios and absolute risks from a
  Cox model with complex nonlinear interactions”* published in the
  American Journal of Epidemiology (2024)

- [R
  material](https://timi.org/wp-content/uploads/2023/11/R-material.zip)

- [SAS
  material](https://timi.org/wp-content/uploads/2023/09/SAS-macros_rev.txt)

### 3. Generalized linear model (GLM) for estimating event risks using Pseudo-values

#### Pseudo-values derivation and comparison with Kaplan-Meier estimates

- [SAS
  Macro](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/Deriving_Pseudo_Values_Macro.sas)

- [ENAR 2025
  Poster](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/ENAR_poster_final.pdf)

#### Pseudo-values GLM with RCS to estimate nonlinear Absolute Risks and Absolute Risk Differences

- [SAS
  Macro](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/Macro_Splines_Pseudo.sas)
  to model nonlinear effects using pseudo-values

- [Macro](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/Macro_Splines_Cox.sas)
  to model nonlinear effects using Cox regression

- [Simulation
  Example](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/SAS_Simulation_Example.sas)

- [Simulation
  Data](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/simdatafinal.sas7bdat)

- [JSM 2023
  Poster](https://github.com/andreabellavia/RCSplines/blob/main/continuous/pseudovalues/JSM_poster_final_HR.pdf)
