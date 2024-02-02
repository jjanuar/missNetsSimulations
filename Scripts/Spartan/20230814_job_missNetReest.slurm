#!/bin/bash

#SBATCH --job-name=20230814_missNetReest
#SBATCH -p cascade
#SBATCH --time=06:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=16384
#SBATCH --array=1-50
#SBATCH -o /home/jjanuar/covNets/Output/20230825_simMissReest/Slurms/jon%j
#SBATCH --export=ALL,missSaveInd=2,propMissVal=0.6

netInd=$1
missModInd=$2


module load foss/2022a R/4.2.1

echo "start job"


cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20230814_simMissReest.R $netInd $missModInd

echo "end job"