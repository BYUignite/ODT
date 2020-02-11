#!/bin/bash

if [ -z $1 ] ; then
    echo "ERROR: Need a directory argument"
    exit 0
fi
caseDir=$1

module load matlab

Ld=`../tools/getInputFileParameter.py ../../input/odtParam.inp domainLength`

echo "set bar 0;"                  > plotEddyMap.gnu
echo "set xlabel('Position (m)')" >> plotEddyMap.gnu
echo "set ylabel('Time (s)')"     >> plotEddyMap.gnu
echo "set xrange [0:$Ld"]         >> plotEddyMap.gnu

nRlz=`ls ../../run/RUNTIME/$caseDir/runtime_* | wc | awk '{print $1}'`
nRlz=40


i=0;
while [ $i -lt "$nRlz" ] ; do

echo $i
grep "^ *[1-9]" "../../run/RUNTIME/$caseDir/runtime_$i" | awk '{print $7 "  "$2"  " $6"  "}' > "eddySpaceTime/eddySpaceTime_R$i.dat"

echo "plot 'eddySpaceTime/eddySpaceTime_R$i.dat' us 1:2:3 with xerrorbars ti 'R_$i'; pause -1" >> plotEddyMap.gnu

i=$[$i+1]
done

./fixErrorBar.py

