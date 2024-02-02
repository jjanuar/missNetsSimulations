#!/bin/bash

#SBATCH --job-name=20231024_simMissSubset
#SBATCH -p cascade
#SBATCH --time=01:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=16384
#SBATCH --array=1-3

module load foss/2022a R/4.2.1

echo "start job"

cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20231024_simMissSubset.R 

echo "end job"