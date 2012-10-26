#!/bin/bash

combined=`grep TOTEN OUTCAR | tail -1 | awk '{print $5}'`
molecule=`grep TOTEN chg/mol/OUTCAR | tail -1 | awk '{print $5}'`
surface=`grep TOTEN chg/sur/OUTCAR | tail -1 | awk '{print $5}'`

echo "combined: $combined"
echo "molecule: $molecule"
echo "surface: $surface"
echo "----------------"

binding=`echo "$combined - $molecule - $surface" | bc`

echo "binding: $binding"

grep IDIPOL INCAR chg/mol/INCAR chg/sur/INCAR
grep LDIPOL INCAR chg/mol/INCAR chg/sur/INCAR
grep NGXF OUTCAR | head -1
grep NGXF chg/mol/OUTCAR | head -1
grep NGXF chg/sur/OUTCAR | head -1
tail -3 KPOINTS | head -1
tail -3 chg/mol/KPOINTS | head -1
tail -3 chg/sur/KPOINTS | head -1
