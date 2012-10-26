#!/usr/bin/python

#QP = -7.430 - -5.839# BDA, HOMO-1, 829 bands

#QP = -26.005 - -22.514 # BDA, LOMO, 829 bands

#QP = -2.7789 # BDA, HOMO, experiment
#QP = -1.956 # BDA, HOMO, 3330 bands
QP = -1.433 # BDA, HOMO, 829 bands  

print 'REMEMBER TO SEE IF ITS IN DIRECT COORINDATES OR NOT'

image_plane = 1.47 

surface_z = 12.094969999999968

import numpy

inputFile = open('CONTCAR', 'r') # this makes sure it's direct coordinates

for i in range(2):

    inputFile.readline()

a = []

for i in range(3):

    a.append(float(inputFile.readline().split()[i]))

for i in range(3):

    inputFile.readline()


z_coor = []

for i in range(6):

    z_coor.append(a[2]*float(inputFile.readline().split()[2]))

carbon_z = numpy.average(z_coor)
print "carbon_z: ", carbon_z
print "surface_z: ", surface_z
d = ((carbon_z - surface_z) - image_plane)/0.529177
W = -13.605692/(2*d)
print "d: ", d, "W: ", W, "QP: ", QP 
sigma = QP - W 

print "c6_z - surface_z: ", carbon_z - surface_z, "sigma: ", sigma

