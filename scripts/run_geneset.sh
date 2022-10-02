#!/bin/bash
#SBATCH --job-name=geneset
#SBATCH --output=geneset.txt
#SBATCH --ntasks=1
#SBATCH --mail-user=skhodursky@rockefeller.edu
#SBATCH -N 1
#SBATCH -n 24

#source ~/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

conda activate ood


Rscript geneset_hpc_new.R
