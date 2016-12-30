#!/usr/bin/python

import math

# constants

a0 = 5.291772106712E-09 # Bohr radius in cm

pi = math.pi


# inputs

Rs = input('For what value of Rs would you like the box? ')
N = input('How many particles are there in your box? ')

volume = (1.0/3)*(4*N*pi)*(Rs*a0)**3
length = volume**(1.0/3)
density = N/volume

print '\n For ' + repr(N) + ' particles and Rs of ' + str(Rs) + ','
print 'cell must have a volume of ' + repr(volume) + ' cm^3 (' + repr(volume*(1/a0)**3) + ' A.U.)'
print 'which cooresponds to a length of ' + repr(length/1E-8) + ' A (' + repr(length/a0) + ' A.U.),' 
print 'and a density of ' + repr(density) + ' cm^-3'
