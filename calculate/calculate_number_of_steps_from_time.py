#!/usr/bin/python

# backlib.py

import math


# constants

tAU = 2.4188E-17  # 1 A.U. in seconds
tPS = 1E-12       # 1 picosecond in seconds
# inputs

runTimeInPS = input('Desired time (ps)? ')
timeStep = input('Timestep in A.U.? ')

maxStep= (tPS/tAU)*runTimeInPS/timeStep

print '\n For ' + repr(runTimeInPS) + ' ps with a timestep of ' + repr(timeStep) + ',\n you will need ' + repr(int(maxStep)) + ' steps.'
        
