# Validation report

Date: 2026-07-16

## Passed static checks

- All 17 R source files have balanced quoted strings and `()`, `[]`, `{}` delimiters.
- No numbered entry script contains a hard-coded `significance-filtered edge with zero.
- GSE39582 processing contains hard assertions for 566 tumours and 443 discovery reference tumours.
- GSE14333 censoring is explicitly inverted (`event = 1 - censor`).
- The meta-analysis helper stops if GSE14333 and GSE17536 coexist.
- Cox model helpers include `cox.zph`, dfbeta, C-index and score nonlinearity outputs.
- Cox convergence/infinite-coefficient warnings are captured and invalidate coefficient publication.
- Zero-event coefficient levels are audited before results are released.
- Step 10 contains a sex/stage-stratified sensitivity model and does not force Validation secondary models.
- DFBETA output uses true coefficient names and sample identifiers and includes an influence summary.
- No iteration-limit workaround is used to make a separated Cox model appear converged.
- Mutation parsing recognizes the exact single-letter value `M` as mutated.
- Bulk scoring requires exactly 21 locked signature genes.

Run the included check with:

```bash
python corrected_pipeline/tests/static_guardrails.py
```

An executable synthetic R test is also included:

```bash
Rscript corrected_pipeline/tests/test_step10_estimability.R
```

It constructs a zero-event factor level, verifies that the coefficient model is rejected, verifies that the
corresponding stratified sensitivity is estimable, and checks DFBETA labels.

## Not executable in this delivery environment

`R`/`Rscript`, the raw GEO expression matrices, the Seurat objects and the pinned CellChat/NicheNet resources were not all available here. The supplied pre-fix Step 10 tables were audited numerically and reproduced independently with Efron partial likelihood, but the revised R code itself still requires execution on the analysis workstation. This report therefore does not claim successful numerical execution of steps 0–13. Run the static test, the synthetic R test, then rerun Step 10 and review `GSE39582_model_estimability_status.csv` before Step 11.
