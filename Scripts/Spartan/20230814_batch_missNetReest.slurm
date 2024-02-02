#!/bin/bash

for netInd in {3 5}
do
	for missModelInd in {1..3}
	do
	
		sbatch /home/jjanuar/covNets/Scripts/Slurms/20230814_job_missNetReest.slurm $netInd $missModelInd
	
	done	
done