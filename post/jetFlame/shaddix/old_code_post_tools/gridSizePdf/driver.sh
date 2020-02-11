#!/bin/bash

module load matlab

if [ -z $1 ] ; then
    echo "ERROR: Need a directory argument"
    exit 0
fi
caseDir=$1

nRlz=`ls ../../run/RUNTIME/$caseDir/runtime_* | wc | awk '{print $1}'`
nDmp=`ls ../../data/$caseDir/data_0/dmp_odtl* | wc | awk '{print $1}'`

rm -f gridSizes/*

COUNTER=0
while [ $COUNTER -lt $nDmp ]; do
    let COUNTER=COUNTER+1
    echo "dump time $COUNTER"
    grep "^ *[1-9]" ../../data/$caseDir/data_[0-9]/dmp_odtl_${COUNTER}.dat | awk '{print $3}' > gridSizes/gridSizes_${COUNTER}.dat
done

matlab -nodisplay -r "nDmp=$nDmp; gridSizePdf"
if [ $? -ne 0 ] ; then
    echo "gridSizePdf.m failed"
    exit 0;
fi

