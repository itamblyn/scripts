#!/bin/bash

rm -f analysis.txt

grep "dielectric constant" au*.out >> analysis.txt
echo >> analysis.txt
grep -A 9 "q-point number" au*.out >> analysis.txt
grep -A 9 "Heads and wings" au*.out >> analysis.txt
