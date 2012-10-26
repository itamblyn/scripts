#!/bin/bash 

a=`cat counter`

mv -f OUTCAR      output/OUTCAR.$a
mv -f POSCAR      output/POSCAR.$a
mv -f CONTCAR     output/CONTCAR.$a
mv -f DOSCAR      output/DOSCAR.$a
mv -f vasprun.xml output/vasprun.$a.xml
mv -f XDATCAR     output/XDATCAR.$a

cp output/CONTCAR.$a POSCAR

mkdir job_files/$a
mv *.o* *.po* *.e* *.pe* job_files/$a
mv machines job_files/$a
mv log job_files/$a

a=`expr $a + 1`
rm -f counter
echo $a > counter

