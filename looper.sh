#!/bin/bash

rm -f eqp0_k0.dat eqp1_k0.dat eqp0.dat eqp1.dat 

for directory in n?_NN

do

  cd $directory
  echo $directory 
  ~/bin/eqp_dft.py eqp0 frq0 sig.log eqp0.dat
  rm -f x??
  head -36 eqp0.dat > k0.dat
  cat eqp0.dat >> ../eqp0.dat
  rm eqp0.dat
  cat k0.dat >> ../eqp0_k0.dat

  ~/bin/eqp_dft.py eqp1 frq0 sig.log eqp1.dat
  rm -f x??
  head -36 eqp1.dat > k1.dat
  cat eqp1.dat >> ../eqp1.dat
  rm eqp1.dat
  cat k1.dat >> ../eqp1_k0.dat

  cd ..

done


echo "# nkptgw: 1  neig: 2880" > eig0.dat
echo "# E(DFT) E(GW)" >> eig0.dat
awk '{print $4, $5}' eqp0.dat >> eig0.dat


echo "# nkptgw: 1  neig: 2880" > eig1.dat
echo "# E(DFT) E(GW)" >> eig1.dat
awk '{print $4, $5}' eqp1.dat >> eig1.dat

echo "1.0" > wtk.tmp

echo "eig0.dat" > dos.in
echo "wtk.tmp" >> dos.in

~/bin/dos_kpt.x < dos.in

mv gw_dos.dat gw_eqp0_dos.dat

echo "eig1.dat" > dos.in
echo "wtk.tmp" >> dos.in

~/bin/dos_kpt.x < dos.in

mv gw_dos.dat gw_eqp1_dos.dat
