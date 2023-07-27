#!/bin/bash
#SBATCH -N1 --ntasks-per-node=8
#SBATCH -t 1:00:00
#SBATCH -p normal_q
#SBATCH -A Precipit

module purge
module load gcc/8.2.0
module load apps site/tinkercliffs/easybuild/setup
module load R/4.1.0-foss-2021a

#export R_LIBS="$HOME/R/gcc/3.6/" # dir for gcc version of package

export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

i=${1:-1}
#j=${2:-0}

R CMD BATCH --no-restore "--args seed=$i" bssanovaHydro.R Rout/bssanovaHydro$i.Rout

