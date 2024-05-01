#!/bin/bash
#SBATCH -t 2-12:30:00
#SBATCH --job-name=NWSimImp
#SBATCH -N 2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --array=1-16

# load any software environment module required for app (e.g. matlab, gcc, cuda)
module load gcc/12.2
module load R/4.2.2

# run my job (e.g. matlab, python)
Rscript ../NWSimulationGen.R ${SLURM_ARRAY_TASK_ID}