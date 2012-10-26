#!/bin/csh

~/scripts/vasp/cat_XDATCAR

cd files

~/scripts/xdat2xyz.pl

~/bin/RDF < ../../../analysis/RDF.NN.in &
~/bin/RDF < ../../../analysis/RDF.NH.in &
~/bin/RDF < ../../../analysis/RDF.HH.in &

~/bin/fortran/unwrap_PBC.x < ../../../analysis/unwrap.in
~/bin/fortran/msd.x < ../../../analysis/msd.in

wait

