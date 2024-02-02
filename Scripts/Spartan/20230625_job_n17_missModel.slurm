#!/bin/bash

#SBATCH --job-name=20230625_n17_miss_models
#SBATCH -p physical
#SBATCH --time=01:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=1024



#SBATCH --mail-user="jjanuar@student.unimelb.edu.au"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-30
#SBATCH -o /home/jjanuar/covNets/Output/20230625_n17_missModels/Slurms/jon%j

missNetInd=$1

module load r/4.0.0

cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20230625_n17_missModels.R $missNetInd