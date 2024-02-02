#!/bin/bash

#SBATCH --job-name=20220616_sagemanERGMs_Todd_Peter
#SBATCH -p physical
#SBATCH --time=08:00:00

#SBATCH --nodes=2
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=4096

#SBATCH --mail-user="jjanuar@student.unimelb.edu.au"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-10


module load r/4.0.0

cd /home/jjanuar/covNets

Rscript --vanilla Scripts/20220616_sagemanERGMs_Todd_Peter.R