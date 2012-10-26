#!/bin/csh

grep -A 28 Summary output.log | grep "\." | awk '{print $2,$3 }' > energy.dat
