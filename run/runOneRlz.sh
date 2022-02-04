#!/bin/bash

# Single realization run
# run as "./runOne.sh" or as "./runOne.sh -r"
# Change variables inputDir and caseName near the top.

###############################################################################
echo "the start time is"
date
###############################################################################

inputDir="../input/jetFlame/sooting_jet"
caseName="test_rad_soot"

###############################################################################

runCase () {

    rm -rf "../data/$caseName" > /dev/null 2>&1
    mkdir  "../data/$caseName"
    mkdir  "../data/$caseName/data"
    mkdir  "../data/$caseName/input"
    mkdir  "../data/$caseName/runtime"
    cp     "$inputDir/"*        "../data/$caseName/input/" > /dev/null 2>&1
    cp -r  "$inputDir/restart"* "../data/$caseName/input/" > /dev/null 2>&1

    #--------------------------------------------------------------------------

    echo "*** RUNNING ***"
    echo "Output is being written to ../$caseName/runtime/runtime_* and ../$caseName/data"
    ../bin/odt-run $caseName 0          # 0 is the shift (realization # here)

}

###############################################################################

rebuild () {
  echo '*** REBUILDING ***'
  cd ../build
  make -j8
  if [ $? -ne 0 ] ; then
    echo ; echo 'FATAL: error in the build' ; echo
    exit 0
  fi
  echo '*** DONE REBUILDING ***'
  cd ../run
}

###############################################################################

if [ "$1" == "-r" ]; then rebuild; fi

runCase "$caseName"

###############################################################################
echo
echo "the end simulation time is"
date
###############################################################################

exit 0
