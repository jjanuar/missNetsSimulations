#!/bin/bash

#SBATCH --job-name=20230825_missErgmReest
#SBATCH -p cascade
#SBATCH --time=02:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=4096
#SBATCH --array=1-30
#SBATCH -o /home/jjanuar/covNets/Output/20230825_simMissReest/Slurms/jon%j
#SBATCH --export=ALL,missSaveInd=2,propMissVal=0.1,ergmCoef=1

# note on exported environmental variables
# missSaveInd
# An index for how the missing data are saved in R
# 1 = NA ("Peter")
# 2 = 0 ("Todd", because nobody likes Todd.)

# propMissVal
# The proportion of missingness.
# only three options for now (0.1, 0.35, 0.6)
# requires simulations to generate this from a separate file

# ergmCoef
# An index for the chosen set of ergm coefficients to have generated the missingness.
# Implies specific and different dependence assumptions in the missingness.
# Need to open the script to view specifics.

# NetInd is specified in the batch SLURM that runs this SLURM
# it goes from 1 to 6 corresponding to the specific network chosen.

netInd=$1


module load foss/2022a R/4.2.1

echo "start job"


cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20230825_simMissErgmReest.R $netInd

echo "end job"