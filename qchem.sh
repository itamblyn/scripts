#!/bin/bash

rm -f machinefile

for i in {1..8}

do

  hostname >> machinefile

done
 

module load QChem/3202
export QCSCRATCH=`pwd`
export PBS_NODEFILE=./machinefile

NP=8

$QC/bin/qchem -pbs -save -np $NP input.in qchem.out scratch 
