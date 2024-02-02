#!/bin/bash

for netInd in {1..6}
do

	sbatch /home/jjanuar/covNets/Scripts/Slurms/20230825_job_missErgmReest.slurm $netInd
			
done