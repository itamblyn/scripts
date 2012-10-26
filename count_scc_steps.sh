#!/bin/bash

grep -B 2 "Total Energy" output.log  | grep -v Total | grep -v "\-\-" | awk '{n+=1; print n, $1; getline}' > scc_steps.dat
