#!/usr/bin/python

import math

# constants

a0 = 0.5291772108E-10  # Bohr radius in m
pi = math.pi
kB = 1.3806503E-23


# inputs

Rs = input('For what value of Rs would you like the box? ')
T = input('Temperature in K? ')

P = (3*kB*T)/(4*pi*(Rs*a0)**3)

print 'Pressure is ' + repr(P/1E9) + ' GPa'
