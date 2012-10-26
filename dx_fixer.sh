#!/bin/bash

grep -B 10000000 "grid connections" $1 > dx_voxels.dat &

grep -A 100000 "grid connections" $1 > dx_parameters.dat

awk '/gridpositions/{print $1, "1", $3, $4, $5, $6, $7, $8}' dx_parameters.dat > output.dx

grep origin dx_parameters.dat >> output.dx
grep delta dx_parameters.dat >> output.dx
awk '/gridconnections/{print $1, "2", $3, $4, $5, $6, $7, $8}' dx_parameters.dat  >> output.dx

wait

awk '/array/{print $1, "3", $3, $4, $5, $6, $7, $8, $9, $10, $11, $12} ! /array/{print $0} ' dx_voxels.dat | grep -v "#" | grep -v "^$" >> output.dx

