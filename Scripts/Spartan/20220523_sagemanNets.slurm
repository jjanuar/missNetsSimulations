#!/bin/bash

#SBATCH --job-name=20220523_sagemanNets
#SBATCH -p physical
#SBATCH --time=02:00:00

#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=2000

#SBATCH --mail-user="jjanuar@student.unimelb.edu.au"
#SBATCH --mail-type=BEGIN,END,FAIL


module load r/4.0.0

cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20220523_sagemanERGMs.R