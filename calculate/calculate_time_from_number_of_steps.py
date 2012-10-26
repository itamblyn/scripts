#!/usr/bin/python

# backlib.py

import math


# constants

tAU = 2.4188E-17  # 1 A.U. in seconds
tPS = 1E-12       # 1 picosecond in seconds
# inputs

maxStep = input('How many steps? ')
timeStep = input('Timestep in A.U.? ')

runTimeInPS = maxStep*timeStep*tAU/tPS

print '\n For ' + repr(maxStep) + ' steps with a timestep of ' + repr(timeStep) + ',\n simluation time is ' + repr(runTimeInPS) + ' picoseconds.'
        
