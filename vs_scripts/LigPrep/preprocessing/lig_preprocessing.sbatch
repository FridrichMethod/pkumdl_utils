#!/bin/bash
#SBATCH --job-name=lig_preprocessing
#SBATCH --output=lig_preprocessing_%j.out
#SBATCH --error=lig_preprocessing_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhaoyangli@stanford.edu  # Your email address
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12  # CPU cores requested for your job
#SBATCH --mem-per-cpu=4G  # Memory requested for each CPU core
#SBATCH --time=0-06:00:00  # Time limit d-hh:mm:ss

# Load required modules
module load python/3.9.0
module load chemistry
module load py-rdkit/2022.09.1_py39
echo "Successfully loaded py-rdkit/2022.09.1_py39"

echo "Lig_preprocessing: Starting for $SLURM_JOB_NAME"
# -i: input directory, -o: output directory, -p: properties
# -n: number of molecules per file, -d: deduplicate, -s: sort, -t: include title
srun python3 sdf_to_smi.py -i compound_libraries -o lig -p IDNUMBER idnumber "ID Number" -n 100000 -d -s -t
echo "Lig_preprocessing: Completed for $SLURM_JOB_NAME"
