#!/bin/bash

for OUTCAR in OUTCAR.1*
do
nsw=`grep NSW $OUTCAR | awk '{print $3}'` 
timesteps=`grep pressure $OUTCAR | wc -l`
elapsed=`grep Elapsed $OUTCAR | awk '{print $4}'`
nproc=`grep "running on" $OUTCAR | awk '{print $3}'`
per_hour=`echo $elapsed | awk '{printf "%4.2f", $1/3600}'`

echo -n "NSW=" $nsw ", "
echo -n "Actual=" $timesteps ", "
echo -n "Nproc=" $nproc ", "
echo -n "Runtime= "$elapsed 
echo " " $per_hour
 
done  
