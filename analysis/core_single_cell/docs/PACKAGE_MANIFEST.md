# R package manifest

Use a project-local `renv` lockfile after installation. The pipeline checks packages but never installs or upgrades them automatically.

## Core single-cell

- Seurat, Matrix, data.table, future, ggplot2, dplyr, tidyr, tibble, patchwork
- CellChat (pin the same major/minor release for all conditions)
- edgeR
- nichenetr plus versioned `ligand_target_matrix.rds` and `lr_network.rds`

## Bulk and survival

- GEOquery, Biobase, limma, AnnotationDbi, hgu133plus2.db
- survival, broom, metafor, splines (base R), scales

## Recommended reproducibility sequence

```r
renv::init(bare = TRUE)
# Install the packages above explicitly, including a pinned CellChat/nichenetr source.
renv::snapshot()
```

Record `sessionInfo()` from every completed step. The scripts already write it under `11_logs`.

## Step 10 validation sequence

```bash
python corrected_pipeline/tests/static_guardrails.py
Rscript corrected_pipeline/tests/test_step10_estimability.R
```

The second test requires only `survival` and `broom`; it does not use the CRC expression matrices.
