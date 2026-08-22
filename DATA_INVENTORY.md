# Data inventory and integrity status

## Included and verified as supplied

- GSE132465: pseudobulk count matrix, pseudobulk metadata, sample audit, and sample mapping used for the outcome-blind biological-coherence analysis.
- GSE144735: pseudobulk count matrix, pseudobulk metadata, sample audit, and sample mapping used for the outcome-blind biological-coherence analysis.
- Frozen ECM, epithelial-target, CAF-proxy, and claim-control resources used by the packaged scripts.
- Analysis-ready recurrence and sensitivity tables used by the packaged Python analyses.

These are processed or analysis-level inputs, not substitutes for the authoritative cell-level raw deposits.

## Pseudobulk integrity audit

The two compressed count matrices were read cell by cell as non-negative integers, and their sample identifiers were compared with the corresponding metadata in the same order.

| Dataset | Genes | Pseudobulk samples | Metadata rows | Count values checked | Sample order matches metadata | Duplicate gene/sample IDs | Invalid or negative counts |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| GSE132465 | 25,655 | 164 | 164 | 4,207,420 | Yes | 0 / 0 | 0 |
| GSE144735 | 30,191 | 101 | 101 | 3,049,291 | Yes | 0 / 0 | 0 |

## Deliberately excluded after integrity screening

- `GSE39582_series_matrix.txt.gz`: the uploaded cache was 959,599 bytes. The supplied module script defines 10 MiB as an early corruption threshold for the complete 585-sample matrix; therefore this cache was excluded. Download the complete public series matrix from GSE39582 and place it at `data/available/GSE39582/GSE39582_series_matrix.txt.gz` before rerunning the module pipeline.

## Referenced rather than redistributed

Cell-level matrices, spatial HDF5/parquet files, perturbational datasets, and raw proteomics deposits not present in the uploaded files are listed by accession in `DATASET_ACCESSIONS_AND_LINKS.md`. Their absence is declared explicitly so that package completeness is not overstated.
