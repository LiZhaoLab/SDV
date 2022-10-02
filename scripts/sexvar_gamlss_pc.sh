#!/bin/bash
#SBATCH --job-name=sbvar
#SBATCH --output=/dev/null
#SBATCH --ntasks=1
#SBATCH --mail-user=skhodursky@rockefeller.edu
#SBATCH --array=43
#SBATCH --cpus-per-task=12

#conda init bash
#source ~/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate ood 

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p tissues.txt)

Rscript  sexvar_gamlss_pc.R $LINE
