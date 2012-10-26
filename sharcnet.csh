#!/bin/csh

foreach n (`awk '{print $1}' ~/scripts/sharcnet.txt`)
echo -n $n 
grep $n ~/scripts/sharcnet.txt | awk '{print " (" $3 " processors)"}'
set processors=`grep $n ~/scripts/sharcnet.txt | awk '{print $3}'`
ssh itamblyn@$n qstat | grep run > status.$n
awk '{r+=$1;w+=$3; print r" running (" '$processors' - r " free), " w " waiting"}' status.$n | tail -1
rm status.$n
echo 
end

