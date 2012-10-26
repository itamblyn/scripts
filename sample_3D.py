#! /usr/bin/env python

# this file now requires that hist files ions+molecules.hist are normalised to begin with. 
# note that ions.hist and molecules.hist are NOT normalised, but sum to ions+molecules.hist 

import sys

if len(sys.argv) == 1: 
     print 'usage: ' + sys.argv[0] + ' input.hist [min]'

if len(sys.argv) == 3: threshold = float(sys.argv[2])
else: threshold = 0.0

import numpy

r_min = 1E6
r_max = float(0)
phi_min = 1E6
phi_max = float(0)

HIST_sum = float(0)

number_of_r_bins = 99 
number_of_phi_bins = 30

inputFile_HIST = open (sys.argv[1],'r')
line = inputFile_HIST.readline()

r_step = float(line.split()[5])
phi_step = float(line.split()[7])*(numpy.pi/180.0)

HIST_array = []
#read line into array
for line in inputFile_HIST.readlines():
    j = 0
    for i in line.split():
        if j == 0:
           if float(i) < r_min: r_min = float(i)
           if float(i) > r_max: r_max = float(i)
        if j == 1:
           if float(i) < phi_min: phi_min = float(i)
           if float(i) > phi_max: phi_max = float(i)
        if j == 2: 
           HIST_array.append(float(i))
           HIST_sum += float(i)
        j += 1
inputFile_HIST.close()

r_max = 3.5

matrix = numpy.zeros((number_of_r_bins,number_of_phi_bins), dtype=numpy.float)

for r_count in range(number_of_r_bins):
    for phi_count in range(number_of_phi_bins):
         matrix[r_count][phi_count] = HIST_array[r_count*number_of_phi_bins + phi_count]

number_of_x_bins = 128 
number_of_y_bins = number_of_x_bins
number_of_z_bins = number_of_x_bins


x_step = (2*r_max)/number_of_x_bins  
y_step = x_step
z_step = x_step  

cart_matrix = numpy.zeros((number_of_x_bins,number_of_y_bins, number_of_z_bins), dtype=numpy.float)

for x_count in range(number_of_x_bins):

    print (float(x_count)/float(number_of_x_bins))*100

    for y_count in range(number_of_y_bins):

        for z_count in range(number_of_z_bins):

            x = float((x_count*x_step + 0.5*x_step) - r_max)
            y = float((y_count*y_step + 0.5*y_step) - r_max)
            z = float((z_count*z_step + 0.5*z_step) - r_max)

            r = (x**2 + y**2 + z**2)**0.5
            phi = numpy.arccos(z/r)

            r_bin = int(r/r_step)
            phi_bin = int(phi/phi_step)

            if r_bin < number_of_r_bins: occ = matrix[r_bin][phi_bin]
            else: occ = 0

            if occ > threshold: cart_matrix[x_count][y_count][z_count] = occ - threshold 
            else: cart_matrix[x_count][y_count][z_count] = 0


cubeFile = open ('test.cube', 'w')

cubeFile.write('COMMENT\n')
cubeFile.write('COMMENT\n')
cubeFile.write('2 0.000000 0.000000 0.000000\n')

cubeFile.write(str(number_of_x_bins) + ' ' + str(x_step) + ' 0.000000 ' + ' 0.000000 ' + '\n')
cubeFile.write(str(number_of_y_bins) + ' ' + ' 0.000000 ' + str(y_step) + ' 0.000000 ' + '\n')
cubeFile.write(str(number_of_z_bins) + ' ' + ' 0.000000 ' + ' 0.000000 ' + str(z_step) + '\n')

cubeFile.write('1 1.000000 3.500000 3.500000 4.200000\n')
cubeFile.write('1 1.000000 3.500000 3.500000 2.800000\n')

for ix in range(number_of_x_bins):
   for iy in range(number_of_y_bins):
      for iz in range(number_of_z_bins):
         voxel_value = cart_matrix[ix][iy][iz] 
         cubeFile.write('% .6E' % voxel_value + ' ')
         if (iz % 6 == 5):
            cubeFile.write('\n')
      cubeFile.write('\n');

cubeFile.close()
