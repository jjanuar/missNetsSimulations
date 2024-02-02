#!/bin/bash

#SBATCH --job-name=20230726_missNetsTrueModels
#SBATCH -p cascade
#SBATCH --time=01:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=1024


#SBATCH --mail-user="j.januar@student.unimelb.edu.au"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-5
#SBATCH -o /home/jjanuar/covNets/Output/20230726_missNetsTrueModels/Slurms/jon%j

netInd=$1

#SBATCH --export=netInd=$netInd

module load foss/2022a R/4.2.1

echo "start job"

cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20230726_missNetsTrueModels.R $netInd

echo "end job"