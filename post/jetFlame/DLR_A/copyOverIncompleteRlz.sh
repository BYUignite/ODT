#!/bin/bash

CASENAME=$1
RLZTOCOPYOVER=$2
RLZTOCOPYFROM=`expr $2 - 1`

RLZ=$(printf "%05d" $RLZTOCOPYOVER)
COPYFROM=$(printf "%05d" $RLZTOCOPYFROM)

echo Copying ../../../data/$CASENAME/data/data_$COPYFROM to ../../../data/$CASENAME/data/data_$RLZ

rm -rf ../../../data/$CASENAME/data/data_$RLZ

cp -r ../../../data/$CASENAME/data/data_$COPYFROM ../../../data/$CASENAME/data/data_$RLZ
