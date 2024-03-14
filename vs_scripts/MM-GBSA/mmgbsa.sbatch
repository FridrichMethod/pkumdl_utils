#!/bin/bash
#SBATCH --job-name=mmgbsa
#SBATCH --output=mmgbsa_%A_%a.out
#SBATCH --error=mmgbsa_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhaoyangli@stanford.edu
#SBATCH --partition=normal
#SBATCH --array=0-9
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2-00:00:00

n=$SLURM_ARRAY_TASK_ID

module load chemistry
module load schrodinger/2022-3
echo "Successfully loaded schrodinger/2022-3"

echo "Prime: Starting for $SLURM_JOB_NAME #$SLURM_ARRAY_TASK_ID"
srun prime_mmgbsa "lig_${n}.maegz" -report_prime_log yes -csv_output yes -out_type COMPLEX -flexdist 5.0 -job_type REAL_MIN -num_output_struct 1 -prime_opt OPLS_VERSION=S-OPLS -HOST "$HOSTNAME:$SLURM_CPUS_PER_TASK" -WAIT -TMPDIR "lig_${n}_tmp" -NJOBS 120 -RETRIES 6
echo "Prime: Completed for $SLURM_JOB_NAME #$SLURM_ARRAY_TASK_ID"
