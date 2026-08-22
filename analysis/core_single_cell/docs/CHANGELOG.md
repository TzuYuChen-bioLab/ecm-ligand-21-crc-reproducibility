# v2.1.0 change log

## 2.1.0 — 2026-07-16

| Script | Main correction | New guardrail |
|---|---|---|
| helpers_bulk | Capture Cox warnings; audit zero-event factor levels and extreme/non-finite estimates | Invalid Cox fits cannot publish coefficients, PH, DFBETA or C-index as valid results |
| 10 | Keep discovery primary model; add sex/stage-stratified sensitivity in both splits and pooled interaction | Validation secondary models are recorded as not attempted; iteration limits are never used to hide separation |
| 10/helpers_bulk | Label influence output with true coefficient and GSM sample identifiers | DFBETA summary includes a transparent `2/sqrt(n)` heuristic flag |
| 12 | Plot Discovery primary together with explicitly labelled Discovery/Validation stratified sensitivities | Non-estimable Validation ordinary Cox output cannot re-enter the publication figure |
| tests | Add static convergence guardrails and an executable synthetic zero-event test | The synthetic coefficient model must fail while the stratified model must pass |

## 2.0.0 — 2026-07-12

| Script | Main correction | New guardrail |
|---|---|---|
| 0 | Portable root and input inventory | No fixed Windows path |
| 1 | Duplicate gene rows are summed | No invented `make.unique()` symbols |
| 2–4 | Metadata alignment and explicit tissue provenance | No unproven malignant label |
| 5 | Paired patients selected using bilateral epithelial coverage | Neutral stromal labels |
| 6 | CellChat run separately by condition | No missing-to-zero probability imputation |
| 7 | Patient-paired epithelial pseudobulk | No cell-level inferential replicate |
| 8 | Fixed 21-gene equal-weight signature | No outcome-driven membership/weight selection |
| 9 | Tissue filtering before probe selection/scoring | GSE39582 must equal 566 tumours |
| 10 | Prespecified GSE39582 primary model | PH, nonlinearity, influence and incremental checks |
| 11 | Full semicolon GEO parser and overlap detection | GSE14333/GSE17536 cannot coexist in a meta-analysis |
| 12–13 | Figures read corrected outputs only | Figure source data and interpretation inventory |

No original file is overwritten. All computed outputs use `corrected` directories.
