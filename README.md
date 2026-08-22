# ECM_LIGAND_21 CRC reproducibility release

This repository supports the analyses reported in the manuscript **“Functional Decomposition and Cross-Context Evaluation of a 21-Gene Extracellular Matrix Score in Colorectal Cancer.”** It is organized for direct use in GitHub and for versioned archival through Zenodo, and contains no user-specific absolute paths.

## What is included

- Complete supplied analysis scripts, organized by analytical role rather than draft history.
- The processed, analysis-level pseudobulk inputs supplied for GSE132465 and GSE144735.
- Analysis-ready GSE39582 recurrence tables and figure-source values supplied with this project.
- Frozen gene-set and claim-control resources.
- Session-information and run-status records.
- Figure-source data and reconstruction audits for the main and supplementary figures.
- The submission-ready supplementary tables workbook.
- Dataset accession links for large public raw inputs that were not present in the uploaded archive.

## Public-release mapping

The original delivery archives have been unpacked and normalized into one browsable repository:

| Public repository path | Organized source material |
| --- | --- |
| `analysis/`, `data/`, `results/`, and reproducibility documentation | Reproducibility code and analysis-level data |
| `source_data/` | Figure-source archive, with user-specific absolute paths replaced by package-relative paths |
| `supplementary/` | `Supplementary_Tables_IJMS_SUBMISSION_READY.xlsx` |

This repository was published as GitHub release `v1.0.0` and archived through the GitHub–Zenodo integration at https://doi.org/10.5281/zenodo.22052359                                                                                          . Separate copies of the original ZIP archives are not attached to the release because doing so would duplicate the same content and make the authoritative version ambiguous.

The archive does **not** claim that every cell-level GEO raw matrix, the complete GSE39582 series matrix, or raw mass-spectrometry file is redistributed here. Those public inputs remain available at their original repositories and are listed in `DATASET_ACCESSIONS_AND_LINKS.md`. The supplied 0.96-MB GSE39582 cache failed the analysis script's 10-MB early integrity guard and was therefore excluded rather than mislabeled as usable raw data; see `DATA_INVENTORY.md`.

## Directory layout

```text
analysis/
  core_single_cell/          GSE132465 single-cell, CellChat, NicheNet, bulk, and recurrence pipeline
  external_triangulation/    GSE144735, GSE92945, spatial, perturbational, proteomic, and external analyses
  module_recurrence/         Four-module GSE39582 Cox, bootstrap, Validation, and CAF analyses
  editorial_robustness/      Exact small-sample, small-k, spatial-shift, and 19-vs-21 sensitivities
  biological_coherence/      Outcome-blind pseudobulk specificity and matched-random audit
data/available/               Verified analysis inputs supplied with this project
source_data/                  Figure-source values, reconstruction files, and audits
supplementary/                Submission-ready supplementary tables workbook
results/generated/            Created at runtime; not pre-populated
```

## Reproduction sequence

1. Obtain the public inputs listed in `DATASET_ACCESSIONS_AND_LINKS.md` and place them under the dataset paths expected by `analysis/core_single_cell/scripts/config.R` and `analysis/external_triangulation/config.R`.
2. Run `analysis/core_single_cell/run_pipeline.R` for the core single-cell/communication workflow.
3. Run `analysis/external_triangulation/scripts/00_package_check.R` through `12_CBC_reproducibility_transportability_specificity_summary.R` in numerical order.
4. Run `analysis/biological_coherence/scripts/02_run_provenance_audit.R` with the two included pseudobulk input folders.
5. Download the complete GSE39582 series matrix to `data/available/GSE39582/GSE39582_series_matrix.txt.gz`; then run `analysis/module_recurrence/scripts/RUN_GSE39582_MODULE_COX_AND_BOOTSTRAP_DELTA_CINDEX.R`, followed by the CAF/Validation continuation script.
6. Run the scripts in `analysis/editorial_robustness/` for exact and small-k sensitivities.
7. Use the files under `source_data/` to audit or redraw the corresponding main and supplementary figures.

## Path configuration

No personal path is embedded. The core pipeline reads `CRC_PROJECT_ROOT`; the triangulation pipeline reads `CRC_TRIANGULATION_ROOT`; the module pipeline uses its own package-relative GSE39582 input by default and accepts `CBC_REPRO_ROOT` for generated outputs.

## Integrity and software

`FILE_MANIFEST.csv` and `SHA256SUMS.txt` cover all packaged files. Module-specific session records are preserved because the uploaded analyses were not produced under one common lockfile. Review `SOFTWARE_ENVIRONMENTS.md` before rerunning and archive any newly resolved package environment separately.

## Licensing

This repository uses scoped mixed licensing:

- Original R and Python source code is licensed under the MIT License; see `LICENSE-CODE`.
- Author-generated documentation, annotations, tables, figures, supplementary materials, and the authors' selection and arrangement of derived results are licensed under CC BY 4.0 to the extent that the authors hold the applicable rights; see `LICENSE-DATA`.
- Underlying third-party and third-party-derived data are not relicensed. Their original source terms continue to apply; see `LICENSE`, `THIRD_PARTY_DATA_NOTICE.md`, and the dataset inventory.

Zenodo records this repository as software under the MIT License. File-specific CC BY 4.0 terms and third-party-data exclusions remain governed by `LICENSE-DATA`, `LICENSE`, and `THIRD_PARTY_DATA_NOTICE.md`.


## Repository metadata

``CITATION.cff` provides citation metadata for this release. The development repository is available at https://github.com/TzuYuChen-bioLab/ecm-ligand-21-crc-reproducibility. The version-specific archival release is available at https://doi.org/10.5281/zenodo.22052359                                                                                                                                                                  . The final article DOI will be added when available.
