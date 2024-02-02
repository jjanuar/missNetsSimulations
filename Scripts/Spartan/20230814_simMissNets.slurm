#!/bin/bash

#SBATCH --job-name=simMissTest
#SBATCH -p physical
#SBATCH --time=01:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=4096

#### Pass all of your environment variables through to the script
#SBATCH --export=ALL,ergmCoef=2,latentBeta=0.4
#SBATCH --array=1-3
#SBATCH --mail-user="j.januar@student.unimelb.edu.au"
#SBATCH --mail-type=BEGIN,END,FAIL

module load r/4.0.0
module load r/4.1.0

echo "start job"

cd /home/jjanuar/covNets
Rscript --vanilla Scripts/20230814_simMissNets.R 

echo "job complete"