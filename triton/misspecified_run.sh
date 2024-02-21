#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1
#SBATCH --output=slurm/%A_%a.out

# load the required modules
module load r

# run experiment with TVD and elpd
srun Rscript ./R/normal-location-misspecified/normal-location-misspecified-elpd.R $1 $2 $3
srun Rscript ./R/normal-location-misspecified/normal-location-misspecified-tvd.R $1 $2 $3
