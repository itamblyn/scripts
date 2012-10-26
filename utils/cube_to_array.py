#! /usr/bin/env python

import sys
import numpy 

inputFile = open (sys.argv[1],'r')

inputFile.readline()
inputFile.readline()

number_of_atoms = int(inputFile.readline().split()[0])

print number_of_atoms

nx = int(inputFile.readline().split()[0])
ny = int(inputFile.readline().split()[0])
nz = int(inputFile.readline().split()[0])

print nx, ny, nz

for i in range(number_of_atoms):
     inputFile.readline()

INPUT_array = []

#read line into array 
for line in inputFile.readlines():

    # loop over the elemets, split by whitespace
    for i in line.split():
    
        INPUT_array.append(float(i))

inputFile.close()

OUTPUT_array = numpy.zeros((nx, ny, nz), dtype=numpy.float)

max_value = 0.0
min_value = 1E6

INPUT_counter = 0

for i in range(nx):
     for j in range(ny):
          for k in range(nz):
               OUTPUT_array[i][j][k] = INPUT_array[INPUT_counter]
               max_value = max(max_value, INPUT_array[INPUT_counter])
               min_value = min(min_value, INPUT_array[INPUT_counter])
               INPUT_counter +=1

# We now have a 3D array with all of the data

print min_value, max_value

# for i in range(nx):
#
#     outputFile = open ('plane-YZ.' + str(i + 1) + '.dat', 'w')
#
#     for j in range(ny):
#          for k in range(nz):
#               outputFile.write(repr(OUTPUT_array[i][j][k]) + ' ')
#          outputFile.write('\n')
#
#     outputFile.close()


OUTPUT_array -= min_value
OUTPUT_array /= max_value
OUTPUT_array *= 255

i = nx/2

outputFile = open('plane-YZ.' + str(i + 1) + '.pgm', 'w')
outputFile.write('P2\n')
outputFile.write('# comment line\n')
outputFile.write(str(nz) + ' ' + str(ny) + '\n')
outputFile.write('255\n')

line_counter = 0

for j in range(ny):
     for k in range(nz):
          outputFile.write(str(int(OUTPUT_array[i][j][k])) + ' ')
          line_counter += 1
          if (line_counter == 17):
               outputFile.write('\n')
               line_counter = 0
outputFile.close()

j = ny/2

outputFile = open('plane-XZ.' + str(i + 1) + '.pgm', 'w')
outputFile.write('P2\n')
outputFile.write('# comment line\n')
outputFile.write(str(nz) + ' ' + str(nx) + '\n')
outputFile.write('255\n')

line_counter = 0

for i in range(nx):
     for k in range(nz):
          outputFile.write(str(int(OUTPUT_array[i][j][k])) + ' ')
          line_counter += 1
          if (line_counter == 17):
               outputFile.write('\n')
               line_counter = 0
outputFile.close()

k = nz/2

outputFile = open('plane-XY.' + str(i + 1) + '.pgm', 'w')
outputFile.write('P2\n')
outputFile.write('# comment line\n')
outputFile.write(str(ny) + ' ' + str(nx) + '\n')
outputFile.write('255\n')

line_counter = 0

for i in range(nx):
     for j in range(ny):
          outputFile.write(str(int(OUTPUT_array[i][j][k])) + ' ')
          line_counter += 1
          if (line_counter == 17):
               outputFile.write('\n')
               line_counter = 0
outputFile.close()
