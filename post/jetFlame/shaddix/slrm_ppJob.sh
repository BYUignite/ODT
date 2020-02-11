#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=01:00:00                    # walltime
#SBATCH --ntasks=1                         # number of processor cores (i.e. tasks) 288
##SBATCH --nodes=12                         # number of nodes  12
#SBATCH --mem-per-cpu=500M                 # memory per CPU core
#SBATCH -J "PPcoldPropaneJet"                # job name
#SBATCH --gid=fslg_crfrg
#SBATCH --mail-user=davidlignell@byu.edu  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

###################################################

module load python/3/5     

###################################################


###################################################
echo "the start time is"
date
###################################################

python driver.py pj_2

####################################################
echo "the end time is"
date
###################################################

exit 0

