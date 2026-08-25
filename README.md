# Functional Decomposition of an Author-Defined 21-Gene Extracellular Matrix Composite Characterizes Compartment-Associated Transcriptional Patterns in Colorectal Cancer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22052359.svg)](https://doi.org/10.5281/zenodo.22052359)

**Version-specific archive:** [Zenodo v1.0.0 — 10.5281/zenodo.22052359](https://doi.org/10.5281/zenodo.22052359) 

**GitHub release:** [v1.0.0](https://github.com/TzuYuChen-bioLab/ecm-ligand-21-crc-reproducibility/releases/tag/v1.0.0)

## Overview

This repository contains the reproducibility materials supporting the manuscript:

> **Functional Decomposition of an Author-Defined 21-Gene Extracellular Matrix Composite Characterizes Compartment-Associated Transcriptional Patterns in Colorectal Cancer**

The study evaluates `ECM_LIGAND_21` as an author-defined, fixed, equally weighted 21-gene extracellular matrix research composite. Panel membership, gene weights, and module assignments were fixed before recurrence modeling.

The analyses examine:

* the cellular and functional composition of the composite;
* paired tumor–normal patterns in colorectal cancer transcriptomic cohorts;
* fibroblast specificity relative to expression- and detection-matched gene panels;
* spatial covariation between neighborhood extracellular matrix scores and epithelial target programs;
* ligand–receptor, perturbational, and proteomic context;
* the contribution of the 19-gene structural core and the `PTN`/`MDK` extension;
* recurrence associations and their robustness to simpler extracellular matrix representations and correlated stromal features.

The repository is intended to support reproducibility of the reported analyses. It does not present `ECM_LIGAND_21` as a clinically validated biomarker.

## Current release

| Item                   | Value                                                                                                                       |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Release                | `v1.0.0`                                                                                                                    |
| Release date           | 2026-08-22                                                                                                                  |
| Resource type          | Software and analysis-level reproducibility materials                                                                       |
| Version-specific DOI   | [10.5281/zenodo.22052359
        
                                            ](https://doi.org/10.5281/zenodo.22052359
        
        )                      |
| Development repository | [TzuYuChen-bioLab/ecm-ligand-21-crc-reproducibility](https://github.com/TzuYuChen-bioLab/ecm-ligand-21-crc-reproducibility) |
| R version              | R 4.5.0                                                                                                                     |
| Python version         | Python 3.12.13                                                                                                              |

## Panel definition

### ECM_LIGAND_21

`ECM_LIGAND_21` contains the following 21 genes:

```text
COL1A1, COL1A2, COL4A1, COL4A2, COL4A5,
COL6A1, COL6A2, COL6A3, FN1, LAMA4,
LAMA5, LAMB1, LAMB2, LAMC1, THBS1,
THBS2, HSPG2, TNC, TNXB, PTN, MDK
```

The genes are organized into four non-overlapping functional modules:

| Module | Functional category                    | Genes                                                                              |
| ------ | -------------------------------------- | ---------------------------------------------------------------------------------- |
| M1     | Interstitial and pericellular collagen | `COL1A1`, `COL1A2`, `COL6A1`, `COL6A2`, `COL6A3`                                   |
| M2     | Basement-membrane components           | `COL4A1`, `COL4A2`, `COL4A5`, `LAMA4`, `LAMA5`, `LAMB1`, `LAMB2`, `LAMC1`, `HSPG2` |
| M3     | Matrix-organizing glycoproteins        | `FN1`, `THBS1`, `THBS2`, `TNC`, `TNXB`                                             |
| M4     | Matrisome-associated secreted factors  | `PTN`, `MDK`                                                                       |

`ECM_CORE_19` consists of M1–M3. `ECM_LIGAND_21` consists of `ECM_CORE_19` plus the M4 genes `PTN` and `MDK`.

Human Matrisome Project/NABA annotations were used to assign the functional categories applied in the study. MatrisomeDB 2.0 was used as an additional annotation cross-reference and was not used to select or reselect panel members.

## Repository contents

| Path or file                      | Description                                                                    |
| --------------------------------- | ------------------------------------------------------------------------------ |
| `analysis/`                       | R and Python scripts used for the reported analyses                            |
| `data/available/`                 | Processed, analysis-level inputs that can be redistributed                     |
| `source_data/`                    | Figure-source values, reconstruction materials, and related audits             |
| `supplementary/`                  | Submission-ready supplementary tables and supporting materials                 |
| `CITATION.cff`                    | Machine-readable citation metadata                                             |
| `DATASET_ACCESSIONS_AND_LINKS.md` | Public dataset accessions, source links, and usage notes                       |
| `DATA_INVENTORY.md`               | Inventory of included and externally accessed data objects                     |
| `FILE_MANIFEST.csv`               | File-level manifest and content descriptions                                   |
| `SOFTWARE_ENVIRONMENTS.md`        | R, Python, package, and software-environment information                       |
| `SHA256SUMS.txt`                  | SHA-256 checksums for integrity verification                                   |
| `THIRD_PARTY_DATA_NOTICE.md`      | Terms and restrictions applying to third-party data                            |
| `LICENSE`                         | General repository licensing notice                                            |
| `LICENSE-CODE`                    | License covering eligible author-generated source code                         |
| `LICENSE-DATA`                    | License covering eligible author-generated documentation and derived materials |

## Included materials

The versioned release includes:

* reproducibility code for the reported analyses;
* processed, analysis-level pseudobulk inputs for `GSE132465` and `GSE144735`;
* analysis-ready recurrence tables and figure-source values for `GSE39582`;
* fixed gene-set definitions and module assignments;
* sample, donor, and dataset mapping audits;
* random seeds and frozen computational resources;
* session-information and software-environment records;
* matched-panel outputs and balance diagnostics;
* model, bootstrap, overlap, and donor-support diagnostics;
* figure-source data and reconstruction audits;
* submission-ready supplementary tables;
* dataset accession links for large public inputs that are not redistributed.

## Public data sources

Public transcriptomic, spatial, perturbational, proteomic, and clinical data used in the study remain available from their original repositories.

### Transcriptomic and spatial data

* `GSE132465`
* `GSE144735`
* `GSE92945`
* `GSE280315`, within SuperSeries `GSE280318`
* `GSE160686`
* `GSE162561`
* `GSE155343`
* `GSE39582`
* `GSE38832`
* `GSE14333`
* `GSE17536`
* `GSE17537`
* `GSE33113`

These datasets are available through the [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/).

### Proteomic data

The decellularized extracellular matrix measurements were obtained from the author-normalized supplementary workbook associated with Lee et al. (2025).

Associated mass-spectrometry data and search-result files are deposited through ProteomeXchange and MassIVE under:

* `PXD037824`
* `MSV000090604`

Dataset-specific links, roles, and usage notes are provided in [`DATASET_ACCESSIONS_AND_LINKS.md`](DATASET_ACCESSIONS_AND_LINKS.md).

## Getting started

### 1. Clone the repository

```bash
git clone https://github.com/TzuYuChen-bioLab/ecm-ligand-21-crc-reproducibility.git
cd ecm-ligand-21-crc-reproducibility
```

Alternatively, download the archived `v1.0.0` release from:

```text
https://doi.org/10.5281/zenodo.22052359
        
        
        
        
```

### 2. Review the repository inventory

Before running an analysis, review:

* [`DATA_INVENTORY.md`](DATA_INVENTORY.md)
* [`FILE_MANIFEST.csv`](FILE_MANIFEST.csv)
* [`DATASET_ACCESSIONS_AND_LINKS.md`](DATASET_ACCESSIONS_AND_LINKS.md)
* [`SOFTWARE_ENVIRONMENTS.md`](SOFTWARE_ENVIRONMENTS.md)

These files identify the analysis-level inputs, externally accessed datasets, software requirements, and expected outputs.

### 3. Verify file integrity

On Linux or macOS, run:

```bash
sha256sum -c SHA256SUMS.txt
```

On Windows PowerShell, individual files can be checked with:

```powershell
Get-FileHash -Algorithm SHA256 <file_path>
```

Compare the resulting values with `SHA256SUMS.txt`.

### 4. Prepare the software environment

The analyses were conducted using:

```text
R 4.5.0
Python 3.12.13
```

Exact package versions and session-information records are provided in [`SOFTWARE_ENVIRONMENTS.md`](SOFTWARE_ENVIRONMENTS.md) and the archived environment files.

Run scripts from the repository root so that repository-relative paths resolve correctly. The public release contains no user-specific local absolute paths.

### 5. Obtain externally hosted inputs when required

Large third-party raw files are not redistributed in this repository. Download them from their original repositories using the accession numbers and links documented in [`DATASET_ACCESSIONS_AND_LINKS.md`](DATASET_ACCESSIONS_AND_LINKS.md).

Do not replace the analysis-level files with differently processed versions unless the purpose is an explicit sensitivity or extension analysis.

### 6. Run the analyses

Scripts are organized under `analysis/` by analytical evidence layer. Follow the script headers, file manifest, and documented input–output mappings.

The release is primarily designed to reproduce the reported analyses from the included analysis-level inputs. Complete reconstruction from every third-party raw dataset may require substantial storage, computation, and dataset-specific preprocessing.

## Reproducibility safeguards

The reproducibility package records the main analytical safeguards used in the study:

* fixed 21-gene membership and equal gene weights;
* fixed four-module assignment;
* explicit `ECM_CORE_19` comparison;
* patient-level paired pseudobulk reporting where applicable;
* outcome-blind expression- and detection-matched gene panels;
* separation of patient-, specimen-, and experimental-model reporting units;
* preservation of modality-specific interpretation boundaries;
* bootstrap confidence intervals and multiplicity adjustments;
* recurrence-model, overlap, and matched-panel diagnostics;
* random seeds and frozen computational resources;
* figure-source and sample-mapping audits.

Evidence from different modalities was not combined into a cross-modal meta-analysis. Within-specimen spatial associations and computational ligand–receptor candidates were not interpreted as patient-level replication or established signaling mechanisms.

## Summary of analytical interpretation

The reported analyses support the following bounded interpretation:

* `ECM_LIGAND_21` showed greater fibroblast specificity than expression- and detection-matched panels in two paired cohorts.
* M1 and M2 had the largest descriptive point estimates in fibroblast and epithelial contrasts, respectively, in both paired cohorts, although direct module superiority was not established.
* The full-score neighborhood measure covaried positively with epithelial target programs within each of three spatial specimens.
* `ECM_CORE_19` closely tracked the 21-gene score in sample rankings and recurrence estimates.
* The validation, discrimination, CAF-adjusted, and matched-panel analyses did not establish incremental prognostic information for `ECM_LIGAND_21`.

The composite should therefore be interpreted as an author-defined research construct for component-resolved biological investigation, rather than as a validated or clinically transferable biomarker.

## Data availability and redistribution

Author-generated analysis code, processed analysis-level data, derived results, figure-source data, and supporting documentation are included when redistribution is permitted.

Large third-party raw datasets are not redistributed. They remain available from their original repositories and are subject to their original licenses, access conditions, and terms of use.

The presence of a public accession number does not transfer ownership or relicensing rights to this repository. See [`THIRD_PARTY_DATA_NOTICE.md`](THIRD_PARTY_DATA_NOTICE.md) for details.

## Licensing

The repository uses separate licenses for different material types:

* Eligible author-generated R and Python source code is licensed under the [MIT License](LICENSE-CODE).
* Eligible author-generated documentation, annotations, tables, figures, supplementary materials, and derived-result arrangements are licensed under [Creative Commons Attribution 4.0 International](LICENSE-DATA).
* Third-party data and materials are not relicensed and remain subject to their original repository or publisher terms.

Review [`LICENSE`](LICENSE), [`LICENSE-CODE`](LICENSE-CODE), [`LICENSE-DATA`](LICENSE-DATA), and [`THIRD_PARTY_DATA_NOTICE.md`](THIRD_PARTY_DATA_NOTICE.md) before reuse or redistribution.

## Citation

If you use the code, analysis-level data, or derived materials in this repository, cite the version-specific Zenodo release:

> Chen, T.-Y., Cai, Y., & Cai, T. (2026). *Functional Decomposition of an Author-Defined 21-Gene Extracellular Matrix Composite Characterizes Compartment-Associated Transcriptional Patterns in Colorectal Cancer* (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.22052359                                    

BibTeX:

```bibtex
@software{Chen_2026_ECM_LIGAND_21,
  author    = {Chen, Tzu-Yu and Cai, Yannan and Cai, Ting},
  title     = {Functional Decomposition of an Author-Defined 21-Gene Extracellular Matrix Composite Characterizes Compartment-Associated Transcriptional Patterns in Colorectal Cancer},
  year      = {2026},
  version   = {1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22052359},
  url       = {https://doi.org/10.5281/zenodo.22052359}
}
```

Machine-readable citation metadata are available in [`CITATION.cff`](CITATION.cff).

When the associated journal article becomes available, please cite both the article and the versioned reproducibility release.

## Authors

* Tzu-Yu Chen
* Yannan Cai
* Ting Cai

## Correspondence

**Ting Cai**
Ningbo No. 2 Hospital
Ningbo, Zhejiang, China

Email:

* [caiting@nbu.edu.cn](mailto:caiting@nbu.edu.cn)
* [caiting@ucas.ac.cn](mailto:caiting@ucas.ac.cn)

## Questions and issues

For questions about the repository, reproducibility package, or file organization, open a GitHub Issue or contact the corresponding author.

When reporting a reproducibility issue, include:

* the repository release or commit identifier;
* the operating system;
* the R or Python version;
* the script or analysis step involved;
* the complete error message;
* whether the versioned Zenodo release or the current GitHub branch was used.
