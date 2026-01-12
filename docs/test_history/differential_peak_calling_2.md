# differential_peak_calling_2 test history

## 2026-01-11 — pytest-suite (Slurm)
- JobID: 13242556
- Run dir: /scratch/users/dhusmann/runs/pytest-suite/20260111_125922
- Failure: tests/test_07_differential.yml::test_chipbinner_use_input_missing_fail expected stderr "input BAMs are missing"; heredoc parsing broke in the workflow (log.err showed `cat: unrecognized option '--samples'`).
- Fix: switched heredocs to `printf` in chipbinner_use_input_missing_{fail,allow_partial} tests to avoid heredoc parsing issues; commit b3341a6.
- Debug rerun (tag test_chipbinner_use_input): JobID 13250774, run dir /scratch/users/dhusmann/runs/pytest-suite/20260111_183621 (COMPLETED 0:0).

## 2026-01-12 — pytest-suite (Slurm)
- JobID: 13250786
- Run dir: /scratch/users/dhusmann/runs/pytest-suite/20260111_183752
- Result: COMPLETED 0:0
