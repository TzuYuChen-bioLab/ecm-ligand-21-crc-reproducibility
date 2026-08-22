# Software environment records

- Core single-cell pipeline: `analysis/core_single_cell/VALIDATION_REPORT.md`, the checks in `analysis/core_single_cell/scripts/0.R`, and the session-writing functions in the pipeline.
- External triangulation: `analysis/external_triangulation/environment/R_sessionInfo.txt` and the status files in `analysis/external_triangulation/run_status/`.
- Module recurrence: `analysis/module_recurrence/environment/R_sessionInfo.txt`.
- Biological-coherence audit: `analysis/biological_coherence/environment/R_sessionInfo.txt`.
- Editorial robustness: the Python scripts record software information with their outputs where implemented.

The archive preserves analysis-specific environment records rather than claiming that one lockfile generated every reported result. Reusers should record a fresh `sessionInfo()` or Python environment when rerunning.
