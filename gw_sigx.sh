#!/bin/bash

touch tmp.out

neig=`awk '/bdgw/{print $3}' *.out | head -1`
nkptgw=`grep nkptgw *.out | head -1 | awk '/nkptgw/{print $3}'`

echo "# nkptgw: $nkptgw neig: $neig" > sigx.dat
echo "# E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E" >> sigx.dat

#grep -A $neig "Band     E0" *.out | grep -v "Band" | grep -v "\-\-" | awk '{print $3, $4, $5}'>> sigx.dat
grep -A $neig "Band     E0" *.out | grep -v "Band" | grep -v "\-\-" | awk '{print $3, $4, $5, $6, $7, $8, $9, $10, $11}'>> sigx.dat

head sigx.dat

rm tmp.out
