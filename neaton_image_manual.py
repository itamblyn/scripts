#!/global/software/centos-5.x86_64/modules/Python/2.7.1/bin/python

import sys
import numpy

QP = -3.1  

image_plane = 1.47
#image_plane = 1.

surface_z = 0.0
carbon_z = surface_z + 3.7

#carbon_z = surface_z + 5.4

print "carbon_z: ", carbon_z
print "surface_z: ", surface_z
print "image_plane: ", image_plane
d = ((carbon_z - surface_z) - image_plane)/0.529177
W = -13.605692/(2*d)
print "d: ", d, "W: ", W, "QP: ", QP 
sigma = QP - W 

print "c6_z - surface_z: ", carbon_z - surface_z, "sigma: ", sigma

