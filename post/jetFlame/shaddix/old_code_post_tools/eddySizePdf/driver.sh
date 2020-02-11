#!/bin/bash

module load matlab

if [ -z $1 ] ; then
    echo "ERROR: Need a directory argument"
    exit 0
fi
caseDir=$1

dmb="../../run/RUNTIME/$caseDir/runtime_*"
grep "^ *[0-9]" $dmb | awk '{print $3 "  " $8 "  " $7 "  "}' > eddyData.dat

matlab -nodisplay -r eddySizePdf
if [ $? -ne 0 ] ; then
    echo "eddySizePdf.m failed"
    exit 0;
fi

matlab -nodisplay -r eddySizePosPdfContours
if [ $? -ne 0 ] ; then
    echo "eddySizePosPdfContours.m failed"
    exit 0;
fi

