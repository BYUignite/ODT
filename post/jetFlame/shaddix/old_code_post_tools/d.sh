#!/bin/bash

# Run as "./driver.sh myCaseDir"
# Script processes ODT data to produce means and conditional means
# FWHM and cL values
# Pictures of these as *.png files
# The eddy PDF and Pictures
# Saves the case to case directories

# This is originally written combusting, temporally-evolving cases.
# Processes data files from dumpTimes.inp (in the data dir).
# Runs Matlab, Octave, gnuplot.

# Be sure to set some parameters at the top of the Matlab files.

#---------- test input case directory from command line

echo "Post-processor normally called after parallel runs with a number of ODT realizations."
echo "  These realizations are stored in ../data/caseDir and ../run/RUNTIME/caseDir,"
echo "  where caseDir is a name provided at runtime."

if [ -z $1 ] ; then
    echo "ERROR: Need a directory argument"
    exit 0
fi
caseDir=$1

if [ ! -d "$caseDir" ] ; then
    mkdir $caseDir
fi

if [ ! -d "../data/$caseDir" ] ; then
    echo "ERROR: ../data/$caseDir not found"
    exit 0
fi
if [ ! -d "../run/RUNTIME/$caseDir" ] ; then
    echo "ERROR: ../run/RUNTIME/$caseDir not found"
    exit 0
fi


#---------- process the data --> means*.dat cmean*.dat sig*.dat csig*.dat


cd eddyMap/
./driver.sh $caseDir
cd ../
cd eddySizePdf/
./driver.sh $caseDir
cd ../
cd gridSizePdf/
./driver.sh $caseDir
cd ../

