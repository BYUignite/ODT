#!/bin/bash

# Submit this script with: sbatch thisFileName caseName
# This is a parallel (ODT MPI) version that runs a parallel job

#SBATCH --time=00:10:00                    # walltime
#SBATCH --ntasks=1                         # number of processor cores (i.e. tasks) 288
#SBATCH -J "coldJet"                       # slurm job name
#SBATCH --gid=fslg_crfrg
#SBATCH --mem-per-cpu=1000M                 # memory per CPU core

module load python/3/5


###############################################################################
echo "the start time is"
date
###############################################################################

caseName=$1

python driver.py $caseName

###############################################################################
echo "the end simulation time is"
date
###############################################################################

exit 0

