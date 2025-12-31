# AGENTS.md

## SLURM job monitoring (Codex sessions)

**What went wrong previously**
- Monitoring used many short-lived `squeue` polls across separate exec sessions; once the 60-session limit was hit, the polling stopped and the job completion was missed.

**Use this single-session loop instead**
1) Submit the job and capture the JobID:
   - `sbatch /path/to/script.sh` → save the JobID
2) Monitor in one long-running session (do not start new sessions for each poll):
   - Example:
     ```bash
     JOBID=123456
     LOG=/scratch/users/dhusmann/nextflow-work/logs/pytest_full_${JOBID}.out
     ERR=/scratch/users/dhusmann/nextflow-work/logs/pytest_full_${JOBID}.err
     while squeue -j "$JOBID" -h | grep -q .; do
       date
       squeue -j "$JOBID" -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
       [ -f "$LOG" ] && tail -n 5 "$LOG"
       [ -f "$ERR" ] && tail -n 5 "$ERR"
       sleep 300
     done
     echo "Job $JOBID no longer in queue"
     sacct -j "$JOBID" -X --format=JobID,State,ExitCode,Elapsed
     ```
3) If `squeue` shows nothing, always confirm completion with `sacct` and then inspect the log files.
