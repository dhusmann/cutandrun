# Differential peak calling test history

## 2026-01-09 — full pytest suite (Slurm)

- **Failure:** `LexerNoViableAltException` / Groovy parsing errors during `span_compare` (seen in Nextflow stderr).  
  **Cause:** unescaped `$` tokens in AWK and regex anchors inside `modules/local/span_compare.nf`.  
  **Fix:** escape `$` in AWK fields and regex anchors (e.g., `\$2`, `\$i`, `/pattern\\$/`).

- **Failure:** Nextflow could not start because `java` was not found in the pytest job environment.  
  **Fix:** set `JAVA_HOME`, `NXF_JAVA_HOME`, and `JAVA_CMD` to `/share/software/user/open/java/17.0.4` in the pytest script.

- **Failure:** pytest workflow missing `PROFILE` environment variable, causing Nextflow config/profile resolution to fail in tests.  
  **Fix:** export `PROFILE=singularity` in the pytest script.

- **Failure:** Nextflow picked up `~/.nextflow/config` with `process.executor=slurm`, leading to nested Slurm execution during pytest runs.  
  **Fix:** set `NXF_HOME` to a scratch directory and `NXF_EXECUTOR=local` in the pytest script to isolate from user config and force local executor.

- **Failure:** `CHIPBINNER_COUNTS` input file name collision (same BAM/BAI filenames for target and control inputs) in `test_differential_chipbinner`.  
  **Fix:** stage control BAM/BAI into a dedicated subdirectory via `stageAs: { "controls/${it.getName()}" }` in `modules/local/chipbinner_counts.nf`.

## 2026-01-10 — differential tag + full pytest suite (Slurm)

- **Failure:** `ROTS` error (`rowSums(is.na(...))` expects matrix) when a condition has <2 replicates in `test_differential_chipbinner` (job `13175581`, run dir `/scratch/users/dhusmann/runs/pytest-suite/20260109_210345`).  
  **Fix:** disable ROTS when <2 samples per condition and fall back to simple tests (or p=1 if <2) in `bin/chipbinner_rots.R`.

- **Failure:** `DIFFERENTIAL_SUMMARY_MERGE` unbound variable due to `$` in staged summary filename (job `13175959`, run dir `/scratch/users/dhusmann/runs/pytest-suite/20260109_212945`).  
  **Fix:** stage summaries with `stageAs: 'summaries/*'` and build `summary_files.list` via `find summaries -type f` in `modules/local/differential_summary.nf`.

- **Failure:** `ANNOTATE_REGIONS` container missing `python` (job `13177002`, run dir `/scratch/users/dhusmann/runs/pytest-suite/20260109_222413`).  
  **Fix:** disable container for `ANNOTATE_REGIONS` so conda provides python (`conf/modules.config`).

- **Failure:** `ANNOTATE_REGIONS` publishDir error (`Unexpected path value` from publishDir closure) (job `13178611`, run dir `/scratch/users/dhusmann/runs/pytest-suite/20260109_225233`).  
  **Fix:** replace publishDir closure with a static map + `enabled` closure in `modules/local/annotate_regions.nf`.

- **Pass:** `test_differential_chipbinner` tag rerun succeeded (job `13179209`, run dir `/scratch/users/dhusmann/runs/pytest-suite/20260109_231631`).
- **Pass:** full pytest suite succeeded (job `13179771`, run dir `/scratch/users/dhusmann/runs/pytest-suite/20260109_234549`).
