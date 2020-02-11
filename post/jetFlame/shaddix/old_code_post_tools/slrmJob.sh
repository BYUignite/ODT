#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=03:00:00                    # walltime
#SBATCH --ntasks=1                       # number of processor cores (i.e. tasks)
#SBATCH --nodes=1                        # number of nodes
#SBATCH --mem-per-cpu=999M                 # memory per CPU core
#SBATCH -J "snlsootPost"                         # job name
#SBATCH --gid=fslg_crfrg
##SBATCH -p dol4                           #dol4 m7 m6beta
#SBATCH --mail-user=davidlignell@gbyu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
#export PBS_QUEUE=batch

###################################################
echo "the start time is"
date
###################################################

caseN="case9b"

module load matlab
module load python/2.7.7

./driver.sh $caseN

echo "the end time is"
date
###################################################

exit 0
