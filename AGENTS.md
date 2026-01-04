# AGENTS.md

## SLURM job monitoring (Codex sessions)

**What went wrong previously**
- Monitoring spawned many short-lived `squeue` polls across separate exec sessions; once the 60-session cap was hit, polling stopped and the job completion was missed.

**Standard monitoring process (use sherlock-slurm-operator helpers)**
1) Submit via the helper so a watcher is created automatically:
   - `~/.codex/skills/sherlock-slurm-operator/scripts/slurm_run.sh -n <jobname> -p hns -c <cpus> -m <mem> -t <time> -- <command>`
   - Capture the printed JobID and run directory.
2) The watcher writes compact status lines to:
   - `<run_dir>/status/latest.txt` (single line)
   - `<run_dir>/status/watch.tsv` (append-only)
3) Monitor without creating new exec sessions:
   - Preferred (blocking, quiet): `~/.codex/skills/sherlock-slurm-operator/scripts/slurm_wait.sh <jobid> 60 <run_dir>`
   - Or periodically read one line: `cat <run_dir>/status/latest.txt`
4) After completion, confirm with:
   - `sacct -j <jobid> -o JobIDRaw,State,ExitCode,Elapsed -P -n`
   - If not `COMPLETED 0:0`, inspect the last ~50 lines of `logs/<jobname>_<jobid>.err`.
