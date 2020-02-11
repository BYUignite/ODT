#!/bin/bash

# Submit this script with: sbatch thisFileName
# This is a parallel (ODT MPI) version that runs a parallel job

#SBATCH --time=02:00:00                    # walltime
#SBATCH --ntasks=256                       # number of processor cores (i.e. tasks) 288
#SBATCH -J "coldJet"                       # slurm job name
#SBATCH --gid=fslg_crfrg
#SBATCH --mem-per-cpu=500M                 # memory per CPU core
##SBATCH --mail-user=you@univ.edu  # email address
##SBATCH --mail-type=BEGIN
##SBATCH --mail-type=END
##SBATCH --mail-type=FAIL


###############################################################################
echo "the start time is"
date
###############################################################################

nSetsToRun=1
inputDir="../input/coldJet"

###############################################################################

runCase () {

    caseName=$1

    if [ -d "../data/$caseName" ]; then echo "../data/$caseName exists already"; exit 1; fi
    mkdir  -p "../data/$caseName/data"
    mkdir  -p "../data/$caseName/input"
    mkdir  -p "../data/$caseName/runtime"
    cp     "$inputDir/"*        "../data/$caseName/input/" > /dev/null 2>&1
    cp -r  "$inputDir/restart"* "../data/$caseName/input/" > /dev/null 2>&1
    if [ "$#" -gt 1 ]; then
        inputParamToChange=$2
        inputValue=$3
        python changeInputParam.py $caseName $inputParamToChange $inputValue
    fi

    #--------------------------------------------------------------------------

    echo "*** RUNNING ***"
    echo "Output is being written to ../$caseName/runtime/runtime_* and ../$caseName/data"

    mpiexec -np $SLURM_NTASKS ./sec.x $caseName 0       # 0 is the realization shift

    nshift=0
    it=1
    while [ $it -lt $nSetsToRun ] ; do
        nshift=$(($nshift + $SLURM_NTASKS))
        it=$(($it + 1))
        mpiexec -np $SLURM_NTASKS ./sec.x $caseName $nshift
    done

    mkdir -p "../data/$caseName/slurm"
    mv slurm-* ../data/$caseName/slurm
    mv core.* ../data/$caseName/slurm
}

###############################################################################

runCase "case_01"
#runCase "test_a" "C_param" "10"

###############################################################################
echo "the end simulation time is"
date
###############################################################################

wait

mkdir -p "../data/$caseName/slurm"
mv *$SLURM_JOB_ID* ../data/$caseName/slurm

exit 0

