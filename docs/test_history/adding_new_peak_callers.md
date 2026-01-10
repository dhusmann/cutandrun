## 2026-01-09 — pytest-suite (Slurm)
- JobID/run dir: 13170828 — /scratch/users/dhusmann/nextflow-work/runs/pytest-suite/20260109_183946
- Command: pytest-suite (pytest-workflow) --maxfail=1 --kwdof --symlink --ga --color=yes (via pytest_suite.sh)
- Result: 323 passed, 0 failed (2:05:41), SLURM COMPLETED 0:0
- Failures observed (prior attempt 13116667): Nextflow compile error "Variable ch_seacr_pairs/ch_macs_pairs already defined" in subworkflows/local/peak_calling_extended.nf; fix in 1491b92 (remove def on ch_*_inputs)
