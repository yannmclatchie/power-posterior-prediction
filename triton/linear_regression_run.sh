#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1
#SBATCH --output=slurm/%A_%a.out

# load the required modules
module load r

# run experiment with TVD and elpd
srun Rscript ./R/linear-regression/linear-regression-elpd.R $1 $2 $3
