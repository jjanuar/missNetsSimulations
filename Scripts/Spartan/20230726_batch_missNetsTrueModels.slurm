#!/bin/bash

for netInd in 6
do

	sbatch /home/jjanuar/covNets/Scripts/Slurms/20230726_job_missNetsTrueModels.slurm $netInd 
	
			
done