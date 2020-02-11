#!/bin/bash

# Submit this script with: sbatch thefilename
# This is a non-parallel (non-MPI) version that spawns multiple jobs for each realization


###############################################################################
echo ""
echo "the start time is"
date
###############################################################################

nRlz=500
inputDir="../input/hips"
caseName="hips_test_03"

###############################################################################

mkdir -p "../data/$caseName/data"
mkdir -p "../data/$caseName/input"
mkdir -p "../data/$caseName/runtime"
if [ ! -f "../data/$caseName/input/odt_input.yaml" ]; then
    cp     "$inputDir/"*        "../data/$caseName/input/" > /dev/null 2>&1
    cp -r  "$inputDir/restart"* "../data/$caseName/input/" > /dev/null 2>&1
fi

###############################################################################

echo "*** RUNNING ***"
echo "Output is being written to ../$caseName/runtime/runtime_* and ../$caseName/data"

it=0
while [ $it -lt $nRlz ] ; do
    echo ""
    echo "-------------------- REALIZATION  $it -----------------------"
    echo ""
    ./sec.x $caseName $it
    it=$(($it + 1))
done

###############################################################################
echo ""
echo "the end simulation time is"
date
###############################################################################

exit 0

