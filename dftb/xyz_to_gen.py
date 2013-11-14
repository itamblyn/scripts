#!/usr/bin/env python2.7
# this is a bit of a hack code to make gen files from a single xyz 


import numpy
import sys

inputFile = open(sys.argv[1],'r')
alat = float(sys.argv[2])

outputFile = open('input.gen','w')

natom = int(inputFile.readline().split()[0])

atom_name = {'C':0, 'O':1, 'H':2}
atom_array = []

for i in range(len(atom_name)):

    atom_array.append([])

inputFile.readline() # skips blank line

for line in inputFile.readlines():
    name = line.split()[0]
    x,y,z = float(line.split()[1]), float(line.split()[2]), float(line.split()[3])
    atom_array[atom_name[name]].append([x,y,z])

type_counter = 1
atom_counter = 1

outputFile.write(str(natom) + ' S\n')

for name in atom_name:
    outputFile.write(name + ' ')

outputFile.write('\n')

for i in range(len(atom_name)):

    for triple in atom_array[i]:
        outputFile.write(str(atom_counter) + ' ' + str(type_counter) + ' ' + str(float(triple[0])) \
        + ' ' + str(float(triple[1])) + ' ' + str(float(triple[2])) + '\n')
        atom_counter += 1

    type_counter += 1

outputFile.write('0.000E+00 0.000E+00 0.000E+00        \n')
outputFile.write(str(alat) + ' 0.000E+00 0.000E+00     \n')
outputFile.write('0.000E+00 ' + str(alat)+ ' 0.000E+00 \n')
outputFile.write('0.000E+00 0.000E+00 ' +str(alat)+  ' \n')
