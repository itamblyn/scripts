#!/bin/tcsh

RDF < RDF.in > /dev/null &

rm -f msd.in
echo "geo_CONTCAR.xyz" > msd.in
msd.x < msd.in > /dev/null &

rm -f pressure.dat
grep Pa output.log | grep -v Pars | grep -v Path | awk '{n+=1; print n, $4/1e9}' > pressure.dat &

rm -f temperature.dat
grep "MD Temperature" output.log | awk '{n+=1; print n, $5}' > temperature.dat &

rm -f tail.xyz 
tail -516000 geo_CONTCAR.xyz > tail.xyz

RDF < RDF.tail.in > /dev/null &
