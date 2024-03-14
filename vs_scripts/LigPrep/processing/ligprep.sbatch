#!/bin/bash
#SBATCH --job-name=ligprep
#SBATCH --output=ligprep_%A_%a.out
#SBATCH --error=ligprep_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhaoyangli@stanford.edu  # Your email address
#SBATCH --partition=normal
#SBATCH --array=0-39  # Job array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32  # CPU cores requested for your job
#SBATCH --mem-per-cpu=6G  # Memory requested for each CPU core
#SBATCH --time=2-00:00:00  # Time limit d-hh:mm:ss

n=$SLURM_ARRAY_TASK_ID

# Load required modules
module load chemistry
module load schrodinger/2022-3
echo "Successfully loaded schrodinger/2022-3"

echo "Ligprep: Starting for $SLURM_JOB_NAME #$SLURM_ARRAY_TASK_ID"
# Detailed explanation of the ligprep command should be found in the Schrodinger documentation
# Only -NSTRUCTS argument should be changed according to the total number of ligands.
srun ligprep -f ligfilter.lff -ma 200 -epikx -ph 7.4 -pht 2.0 -g -bff 16 -NSTRUCTS 200 -RETRIES 2 -HOST "$HOSTNAME:$SLURM_CPUS_PER_TASK" -WAIT -TMPDIR "lig_${n}_tmp" -ismi "lig_${n}.smi" -osd "lig_${n}.sdf"
echo "Ligprep: Completed for $SLURM_JOB_NAME #$SLURM_ARRAY_TASK_ID"
