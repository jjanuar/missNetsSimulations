#!/bin/bash

#SBATCH --job-name=simMissErgm
#SBATCH -p physical
#SBATCH --time=01:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=4096

#### Pass all of your environment variables through to the script
#SBATCH --array=1-3
#SBATCH --export=ALL,ergmCoef=3,ergmSpec=1
#SBATCH -o /home/jjanuar/covNets/Output/20230814_simMissNets/Slurms/jon%j


module load r/4.0.0
module load r/4.1.0

echo "start job"

cd /home/jjanuar/covNets
Rscript --vanilla Scripts/20230816_simMissErgm.R 

echo "job complete"