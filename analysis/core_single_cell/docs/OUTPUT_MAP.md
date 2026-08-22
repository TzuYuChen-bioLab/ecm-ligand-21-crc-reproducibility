# Output map and manuscript replacement rules

| Manuscript item | Use this corrected source | Do not reuse |
|---|---|---|
| Epithelial tumour-vs-normal DE | `GSE132465_epithelial_paired_pseudobulk_DEG.csv` | Seurat cell-level `FindMarkers` P values |
| Tumour-vs-normal communication | `GSE132465_CellChat_donor_supported_comparison.csv` | Pooled mixed-object receiver comparison |
| Mechanistic ligand ranking | `GSE132465_NicheNet_ligand_activities_paired_pseudobulk.csv` | NicheNet based on cell-level DEG |
| Signature | `CRC_ECM_LIGAND_21.csv` | Outcome-optimized modules/weights |
| GSE39582 population | 566 primary tumours | 585 samples including 19 mucosa |
| Main adjusted GSE39582 HR | Discovery row labelled `model_role=PRIMARY`; Validation adjustment is the explicitly labelled sex/stage-stratified sensitivity | Validation ordinary full-adjustment output with infinite stage coefficients; minimum-P or “best” model |
| External forest/meta | mutually exclusive primary and swap sensitivity | GSE14333 + GSE17536 together |
| Diagnostics | `GSE39582_model_estimability_status.csv`, event-level audit, `diagnostics/`, incremental and nonlinearity CSVs | coefficient-only output or any row marked `not_estimable`/`not_attempted` |

Before submission, regenerate all figures and replace every sample count, effect estimate, confidence interval and P value in the manuscript and supplement.
