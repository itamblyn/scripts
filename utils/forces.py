#! /usr/bin/env python

import sys

filename = sys.argv[1]
number_of_particles = int(sys.argv[2])

inputFile = open (filename,'r')

forces_array = []

#read line into array 
for line in inputFile.readlines():

    # add a new sublist
    forces_array.append([])

    # loop over the elemets, split by whitespace
    for i in line.split():
    
        # convert to integer and append to the last
        # element of the list
        forces_array[-1].append(float(i))

inputFile.close()

number_of_snapshots = len(forces_array)/float(number_of_particles)

outputFile = open ('forces.dat', 'w')

s = 0 # counts over snapshots

while s < number_of_snapshots:

     p = 0 # counts over particles
     
     fx = 0.0
     fy = 0.0
     fz = 0.0

     while p < number_of_particles:
          
          fx += forces_array[int(p + s*number_of_particles)][0]**2
          fy += forces_array[int(p + s*number_of_particles)][1]**2
          fz += forces_array[int(p + s*number_of_particles)][2]**2
                        
          p +=1
     
     average_net_force = (fx + fy + fz)/float(number_of_particles)

     outputFile.write(repr(s) + ' ' + repr(average_net_force) + '\n')
     
     s += 1

outputFile.close()
