# Differential peak calling test history (local notes)

## 2026-01-10 — pytest failures resolved (current context)

- **ROTS crash with <2 replicates** (`rowSums(is.na(...))` error in `test_differential_chipbinner`).  
  **Fix:** disable ROTS when a condition has <2 samples and fall back to simple tests (or p=1) in `bin/chipbinner_rots.R`.

- **`DIFFERENTIAL_SUMMARY_MERGE` unbound variable** due to `$` in staged summary filename.  
  **Fix:** stage summaries with `stageAs: 'summaries/*'` and build `summary_files.list` via `find summaries -type f` in `modules/local/differential_summary.nf`.

- **`ANNOTATE_REGIONS` missing python in container** (command not found).  
  **Fix:** disable container for `ANNOTATE_REGIONS` so conda provides python (`conf/modules.config`).

- **`ANNOTATE_REGIONS` publishDir error** (`Unexpected path value` from publishDir closure).  
  **Fix:** replace publishDir closure with a static map + `enabled` closure in `modules/local/annotate_regions.nf`.
