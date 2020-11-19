#!/bin/bash
#SBATCH --job-name=jasa-rev
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --partition=panda
#SBATCH --array=1-1000
echo "$SLURM_ARRAY_TASK_ID"

source ~/.bashrc
spack load -r /bxc56dm
R CMD BATCH '1-simulate.R'
