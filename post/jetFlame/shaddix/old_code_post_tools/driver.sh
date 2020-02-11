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

#cd getConditionalPDFs
#/apps/python/2.7.7/bin/python getConditionalPDFs_on_T.py $caseDir
#/apps/python/2.7.7/bin/python getConditionalPDFs.py      $caseDir
#cd ../
#cd getRadiationProperties
#/apps/python/2.7.7/bin/python getRadiationProperties.py $caseDir
#cd ../
#cd getContourPlots
#cd ../

cd eddyMap/
./driver.sh $caseDir
cd ../
cd eddySizePdf/
./driver.sh $caseDir
cd ../
cd gridSizePdf/
./driver.sh $caseDir
cd ../

#cp -r getConditionalPDFs $casedir
#cp -r getContourPlots $caseDir
#cp -r getRadiationProperties $caseDir
cp -r eddyMap $caseDir
cp -r eddySizePdf $caseDir
cp -r gridSizePdf $caseDir

nRlz=`ls ../run/RUNTIME/$caseDir/runtime_* | wc | awk '{print $1}'`
ntimes=`head -n 1 ../input/dumpTimes.inp | awk '{print $1}'`

/fslapps/matlab/matlab_7.8/bin/matlab -nodisplay -r "caseD='$caseDir'; nrlz=$nRlz; ntimes=$ntimes; procData"
if [ $? -ne 0 ] ; then
    echo "procData.m failed"
    exit 0;
fi
#gnuplot plotmeans.gnu
#gnuplot plotsig.gnu
#gnuplot plotcmean.gnu
#gnuplot plotcsig.gnu

##---------- process the data --> fwhm and cL mixf plots
#
#/fslapps/matlab/matlab_7.8/bin/matlab -nodisplay -nojvm -r get_fwhm_profile_odt
#if [ $? -ne 0 ] ; then
#    echo "get_fwhm_profile_odt.m failed"
#    exit 0;
#fi
#gnuplot plotFWHM.gnu
#
##---------- get the eddy pdf and total number of eddies
#
#dmb="../run/RUNTIME/$caseDir/runtime_*"
#grep "^ *[0-9]" $dmb | awk '{print $6}' > eddySizes.dat
#/fslapps/octave-3.2.3/bin/octave eddyPdf.m
#if [ $? -ne 0 ] ; then
#    echo "eddyPdf.m failed"
#    exit 0;
#fi
#
#gnuplot plotEddyPdf.gnu
#
#echo "Total number of eddies for case" > nEddies.dat
#wc eddySizes.dat | awk '{print $1}' >> nEddies.dat
#echo "Total number of realizations for case" >> nEddies.dat
#ls $dmb | wc | awk '{print $1}' >> nEddies.dat
#
##---------- get the eddy space time plots
#
#mkdir "./$caseDir/eddySpaceTime"
#cd "../run/RUNTIME/$caseDir"
#
#i=0;
#while [ $i -lt "$nRlz" ] ; do
#    echo $i
#    grep "^ *[1-9]" "runtime_$i" | awk '{print $7 "  "$2"  " $6"  "}' > "eddySpaceTime_R$i.dat"
#
#cat << EOF > fixErrorBar.m
#load 'eddySpaceTime_R$i.dat';
#d=eddySpaceTime_R$i;
#d(:,3)=d(:,3)/2;
#save eddySpaceTime_R$i.dat d;
#EOF
#
#/fslapps/octave-3.2.3/bin/octave fixErrorBar.m
#
#i=$[$i+1]
#done
#
#rm -f fixErrorBar.m
#mv eddySpaceTime*.dat "../../../post/$caseDir/eddySpaceTime/"
#cd ../../../post/
#
#
##----------------- move results

cp driver.sh *.m *.gnu fwhm_dns.dat $caseDir
mv *.png eddyPDF.dat nEddies.dat eddySizes.dat fwhm_odt.dat chiStoic.dat chiStoicSig.dat means*.dat cmean*.dat csig*.dat sig*.dat $caseDir
mv uxt.dat $caseDir
cp ../input/* "../input/$caseDir"


