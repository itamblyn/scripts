#! /usr/bin/env python

# the nearest neighbour line/ceiling is based on a bcc, which might not always make sense

import sys

filename = sys.argv[1]
number_of_particles = int(sys.argv[2])
acell = float(sys.argv[3])
print "acell: " + repr(acell)

particles_edge = (float(number_of_particles)/2.0)**(1.0/3.0)

inputFile = open (filename,'r')

TRAJECTORY_array = []

#read line into array 
for line in inputFile.readlines():

    # add a new sublist
    TRAJECTORY_array.append([])

    # loop over the elemets, split by whitespace
    for i in line.split():
    
        # convert to integer and append to the last
        # element of the list
        TRAJECTORY_array[-1].append(float(i))

inputFile.close()

number_of_snapshots = len(TRAJECTORY_array)/float(number_of_particles)

outputFile = open ('drift.dat', 'w')

s = 0 # counts over snapshots

max_displacement = 0.0

nn_distance = (1.0/particles_edge)*(3**(0.5)/2.0)

while s < number_of_snapshots:

     average_counter = 0.0
     
     p = 0 # counts over particles
     
     number_of_displaced_particles = 0

     while p < number_of_particles:
     
          
          dx = TRAJECTORY_array[p][0] - TRAJECTORY_array[s*number_of_particles + p][0]
          dy = TRAJECTORY_array[p][1] - TRAJECTORY_array[s*number_of_particles + p][1]
          dz = TRAJECTORY_array[p][2] - TRAJECTORY_array[s*number_of_particles + p][2]
                        
          distance = (1/acell)*(dx**2 + dy**2 + dz**2)**(0.5)

          if (distance >= nn_distance): number_of_displaced_particles += 1

          average_counter += distance
              
          p +=1
     

     average_displacement = float(average_counter)/float(number_of_particles)     

     fraction_displaced_particles = float(number_of_displaced_particles)/float(number_of_particles)

     outputFile.write(repr(s) + ' ' + repr(average_displacement) + ' ' + repr(nn_distance) + ' ' + repr(fraction_displaced_particles) + '\n')
     
     if (average_displacement > max_displacement): max_displacement = average_displacement

     s += 1

outputFile.close()

print "max_displacement: " + repr(int(round(max_displacement*100.0/nn_distance))) + "%"
