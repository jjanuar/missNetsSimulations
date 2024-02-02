#!/bin/bash


for missModInd in {1..1}
do

	sbatch /home/jjanuar/covNets/Scripts/Slurms/20230625_job_n17_missModel.slurm $missModInd 
	
			
done