#!/usr/bin/python

# backlib.py

import math

# constants

a0 = 5.291772108E-11  # Bohr radius in meters
pi = math.pi


# inputs

length = input('Volume [A.U.^3]? ')**(1.0/3.0)
N = input('Number of electrons? ')


Rs = (length)*(3/(4*N*pi))**(1.0/3.0)

print '\n For ' + repr(N) + ' particles and lattice parameter of ' + repr(length) + ',\n Rs = ' + repr(Rs)
