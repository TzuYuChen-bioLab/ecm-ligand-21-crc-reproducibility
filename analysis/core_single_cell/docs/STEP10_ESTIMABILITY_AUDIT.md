# GSE39582 Step 10 estimability audit

Date: 2026-07-16

## Observed pre-fix failure

The supplied Step 10 result bundle contained repeated `coxph` warnings that variables 4–6 could have infinite
coefficients and one model that exhausted its iterations. In the primary formula
`score_z + age + sex + stage + mmr`, columns 4–6 are the three TNM-stage indicators.

Validation primary complete cases contained 98 patients and 30 recurrence events:

| Stage | Patients | Events |
|---|---:|---:|
| I | 9 | 0 |
| II | 51 | 10 |
| III | 37 | 19 |
| IV | 1 | 1 |

The mutation-complete Validation subset contained 89 patients/27 events; dMMR contained 6 patients/0 events.
These zero-event coefficient levels produce monotone likelihood. Increasing `iter.max` is not a valid correction.

## Results that remained numerically interpretable before refitting

- Discovery prespecified primary score: HR 1.243, 95% CI 1.022–1.513, P=0.0297.
- Validation univariable score: HR 1.252, 95% CI 0.888–1.765, P=0.1996.
- Pooled ordinary score-by-cohort interaction: P=0.789; no detectable heterogeneity, with limited Validation power.
- Discovery score nonlinearity LRT: P=0.668.
- The score did not violate PH in the univariable models. Sex violated PH in Discovery (P=0.00115), Validation
  (P=0.0199), and the pooled model (P=0.0263).
- Maximum absolute score DFBETA was below the descriptive `2/sqrt(n)` heuristic in every supplied model.

The pre-fix Validation adjusted C-index, delta C-index, full-model Wald statistics, PH table, and adjusted
nonlinearity test must not be reported because they came from a non-estimable ordinary Cox model.

## Implemented v2.1.0 decision

1. Preserve the Discovery prespecified primary model and label it as primary.
2. Attempt the ordinary Validation primary model only behind an estimability gate; it is expected to be recorded
   as `not_estimable`, with no coefficient table released.
3. Fit `score_z + age + mmr + strata(sex, stage)` in both splits as a transparent sensitivity analysis. Stage
   stratification avoids estimating separated stage coefficients, and sex stratification addresses its PH violation.
4. Fit a similarly stratified pooled score-by-cohort interaction sensitivity.
5. Do not fit Validation site/CMS/mutation secondary models because 27–30 events cannot support their parameter
   burden and sparse levels. Record the decision and counts instead.
6. Publish model-status, event-level, PH, nonlinearity, DFBETA and incremental-value outputs only when their
   underlying fit passes the same estimability gate.

Independent Efron-partial-likelihood checks predict a stratified score HR of approximately 1.235 in Discovery
(95% CI 1.013–1.504; P≈0.0365) and 1.108 in Validation (95% CI 0.787–1.559; P≈0.558). These are verification
targets, not final manuscript values; the R rerun is authoritative.
