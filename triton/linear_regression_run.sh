#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1
#SBATCH --output=slurm/%A_%a.out

# load the required modules
module load r

# run experiment with TVD and elpd
#srun Rscript ./R/linear-regression/linear-regression-elpd.R $1 $2 $3
srun Rscript ./R/linear-regression/linear-regression-tvd.R $1 $2 $3
