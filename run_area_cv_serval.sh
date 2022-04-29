#!/bin/bash
#SBATCH --job-name=irsf_area_cv    # name of the job
#SBATCH --partition=defq,intel     # partition to be used (defq, gpu or intel)
#SBATCH --time=24:00:00            # walltime (up to 96 hours)
#SBATCH --nodes=1                  # number of nodes
#SBATCH --ntasks-per-node=1        # number of tasks (i.e. parallel processes) to be started
#SBATCH --cpus-per-task=1          # number of cpus required to run the script
#SBATCH --mem-per-cpu=16G	   # memory required for process
#SBATCH --array=1-15%15    	   # set number of total simulations and number that can run simultaneously	  


module load gcc

export LD_LIBRARY_PATH="/home/alston92/software/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/gdal-3.3.0/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/proj-8.0.1/lib:$LD_LIBRARY_PATH"

ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/terra/libs/terra.so
ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/rgdal/libs/rgdal.so

module load R

cd /home/alston92/proj/irsf   # where executable and data is located

list=(/home/alston92/proj/irsf/data/serval*.csv)

date
echo "Initiating script"


if [ -f results/area_cv_summary_serval.csv ]; then
	echo "Results file already exists! continuing..."
else
	echo "creating results file area_cv_summary.csv"
	echo "aid,ind_file,cor_irsf,cor_crsf,kld_irsf,kld_crsf,kld_r2" > results/area_cv_summary_serval.csv
fi

if [ -f results/area_cv_data_serval.csv ]; then
        echo "Results file already exists! continuing..."
else
        echo "creating results file area_cv_data.csv"
        echo "aid,prob_irsf,prob_crsf,emp_count" > results/area_cv_data_serval.csv
fi


Rscript area_cross_validation_serval.R ${list[SLURM_ARRAY_TASK_ID]}     # name of script
echo "Script complete"
date
