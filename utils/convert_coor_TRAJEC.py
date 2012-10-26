#! /usr/bin/env python

######
###
##
## this script takes a file of numbers, and converts it to a TRAJEC.xyz file. 
## It also subtracts off the center of mass motion, if there is any....
## It's pretty slow...
##
###
######


import sys

number_of_particles = int(sys.argv[1])
lattice_constant = float(sys.argv[2])    # lattice constant in Bohr

lattice_constant *= 0.529177

inputFile = open ('md.coor','r')

COOR_array = []

#read line into array 
for line in inputFile.readlines():

    # add a new sublist
    COOR_array.append([])

    # loop over the elemets, split by whitespace
    for i in line.split():
    
        # convert to integer and append to the last
        # element of the list
        COOR_array[-1].append(float(i))

inputFile.close()

number_of_snapshots = int(len(COOR_array)/number_of_particles)

element_symbol = 'H'

# determine original COM position

p = 0 

xCOMi = 0.0 
yCOMi = 0.0
zCOMi = 0.0

while p < number_of_particles:

     xCOMi += COOR_array[p][0]/number_of_particles
     yCOMi += COOR_array[p][1]/number_of_particles
     zCOMi += COOR_array[p][2]/number_of_particles

     p += 1

# translate each step back to orignal COM location

s = 0

while s < number_of_snapshots:

     p = 0

     xCOMt = 0.0
     yCOMt = 0.0
     zCOMt = 0.0

     while p < number_of_particles:
          
         xCOMt += COOR_array[p + s*number_of_particles][0]/number_of_particles
         yCOMt += COOR_array[p + s*number_of_particles][1]/number_of_particles
         zCOMt += COOR_array[p + s*number_of_particles][2]/number_of_particles

         p+= 1

     p = 0

     while p < number_of_particles:

         COOR_array[p + s*number_of_particles][0] -= (xCOMt - xCOMi)
         COOR_array[p + s*number_of_particles][1] -= (yCOMt - yCOMi)
         COOR_array[p + s*number_of_particles][2] -= (zCOMt - zCOMi)

         p+= 1

     s += 1

outputFile_TRAJEC = open('TRAJEC.xyz','w')

s = 0

while s < number_of_snapshots:
     
     outputFile_TRAJEC.write(repr(int(number_of_particles)) + '\n' + repr(int(s+1)) + '\n')
     
     p = 0
     
     while p < number_of_particles:
          
          outputFile_TRAJEC.write(element_symbol + '  ' + repr(COOR_array[s*number_of_particles + p][0]*lattice_constant) + '  ' + repr(COOR_array[s*number_of_particles + p][1]*lattice_constant) + '  ' + repr(COOR_array[s*number_of_particles + p][2]*lattice_constant) + '\n')
          p += 1
          
     s += 1
