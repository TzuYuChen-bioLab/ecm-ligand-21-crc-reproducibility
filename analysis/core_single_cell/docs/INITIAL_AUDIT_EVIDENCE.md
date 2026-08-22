# Initial audit evidence from the supplied tables/scripts

This document records why the corrections were triggered. These are audits of the supplied outputs, not results from the new pipeline.

## Single-cell replicate structure

- Original epithelial comparison: 17,469 tumour-tissue epithelial cells versus 1,070 normal-tissue epithelial cells.
- Tumour epithelial cells came from 23 patient IDs; normal epithelial cells came from 10 patient IDs.
- The old DEG table contains many adjusted P values printed as zero, consistent with an inflated cell-level effective sample size; it does not prove that every biological effect is false, but the P values are not patient-level inference.

## Stromal tissue composition hidden by the CAF label

| Original group | Normal cells | Tumour cells | Normal fraction |
|---|---:|---:|---:|
| CAF_Stromal1 | 972 | 73 | 93.0% |
| CAF_Stromal2 | 236 | 158 | 59.9% |
| CAF_Stromal3 | 744 | 124 | 85.7% |
| CAF_Myofibroblasts | 9 | 1,146 | 0.8% |

Therefore, the generic `CAF_` prefix encoded a biological interpretation that was not supported uniformly across tissue conditions. The corrected scripts use neutral stromal subtype names.

## GSE39582 population leakage

- Supplied score/survival table contained 585 rows.
- The old RFS complete-case subset contained 536 rows and 147 events.
- Of those, 17 were non-tumour mucosa records and 2 were coded as RFS events.
- Removing them left 519 RFS-complete tumours and 145 events in that particular old complete-case analysis.
- The corrected processing population is defined earlier: 443 discovery tumours + 123 validation tumours = 566 tumours before any probe selection or score construction.

## GSE14333 parsing and overlap

- GSE14333 has 290 primary colorectal cancers.
- Its clinical metadata were embedded in one semicolon-separated field, for example: location, Dukes stage, age, sex, DFS time and DFS censoring in a single string.
- Complete parsing yielded 226 valid DFS records and 50 recurrence events (`DFS_Cens=0` means an event; `DFS_Cens=1` means censored).
- Exact clinical/endpoint fingerprint matching produced 131 candidate match rows, representing 129 unique GSE14333 samples and 129 unique GSE17536 samples.
- Across matched rows, the supplied ECM score correlation was approximately 0.934.

This is strong data-level evidence that GSE14333 and GSE17536 overlap. The corrected pipeline re-runs this audit and enforces mutual exclusion in every meta-analysis.

## Missing survival diagnostics

The supplied survival scripts called `coxph`, but the result directories did not contain systematic per-model `cox.zph` tests/plots, dfbeta influence files, continuous-score nonlinearity tests, or clinical-only versus clinical+score likelihood-ratio comparisons. The corrected helper writes all of these for every estimable primary model and labels low-event models as not estimable rather than silently omitting the limitation.
