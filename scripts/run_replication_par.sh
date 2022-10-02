#!/bin/bash
#SBATCH --job-name=gamlssrep
#SBATCH --output=/dev/null
#SBATCH --ntasks=1
#SBATCH --mail-user=skhodursky@rockefeller.edu
#SBATCH --array=1-50
#SBATCH --cpus-per-task=4

eval "$(conda shell.bash hook)"
conda activate ood

i=${SLURM_ARRAY_TASK_ID}
echo $i
Rscript replication_permutation_hpc.R $i
