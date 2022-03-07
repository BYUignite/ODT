#!/bin/bash

# Submit this script with: sbatch thefilename
# This is a non-parallel (non-MPI) version that spawns multiple jobs for each realization

#SBATCH --time=00:01:00                   # walltime
#SBATCH --ntasks=1                       # number of processor cores (i.e. tasks) 288
#SBATCH --array=0-7  
#SBATCH -J "test"                       # slurm job name
#SBATCH --mem-per-cpu=250M                 # memory per CPU core

###############################################################################
echo "the start time is"
date
###############################################################################

nSetsToRun=1  
inputDir="../input/jetFlame/shaddix"
caseName="test01"

###############################################################################

runCase () {

    caseName=$1

    mkdir -p "../data/$caseName/data"
    mkdir -p "../data/$caseName/input"
    mkdir -p "../data/$caseName/runtime"
    if [ ! -f "../data/$caseName/input/odt_input.yaml" ]; then
        cp     "$inputDir/"*        "../data/$caseName/input/" > /dev/null 2>&1
        cp -r  "$inputDir/restart"* "../data/$caseName/input/" > /dev/null 2>&1
        if [ "$#" -gt 1 ]; then
            inputParamToChange=$2
            inputValue=$3
            python changeInputParam.py $caseName $inputParamToChange $inputValue
        fi
    fi

    #--------------------------------------------------------------------------

    echo "*** RUNNING ***"
    echo "Output will be written to ../$caseName/runtime/ and ../$caseName/data"

    ../bin/odt-run $caseName $SLURM_ARRAY_TASK_ID

    nshift=0
    it=1
    while [ $it -lt $nSetsToRun ] ; do
        nshift=$(($nshift + $SLURM_ARRAY_TASK_COUNT))
        it=$(($it + 1))
        ../bin/odt-run $caseName $(($SLURM_ARRAY_TASK_ID + $nshift))
    done

}
    

###############################################################################

runCase $caseName
#runCase $caseName "C_param" "5"
#runCase $caseName "C_param" "10"

###############################################################################
echo "the end simulation time is"
date
###############################################################################

wait

mkdir -p "../data/$caseName/slurm"
mv *$SLURM_JOB_ID* ../data/$caseName/slurm

exit 0