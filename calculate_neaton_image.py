#!/global/software/centos-5.x86_64/modules/Python/2.6.5/bin/python

#QP = -7.430 - -5.839# BDA, HOMO-1, 829 bands
#QP = -26.315 - -22.514 # BDA LOMO 750, SR
#QP = -26.035 - -22.514 #, BDA, LOMO, 1024 bands
QP = -26.079 - -22.514 # BDA, LOMO, 1424
#QP = -26.005 - -22.514 # BDA, LOMO, 829 bands



#QP = -2.7789 # BDA, HOMO, experiment
#QP = -2.73 # BDA, NL
#QP = -1.956 # BDA, HOMO, 3330 bands
#QP = -1.433 # BDA, HOMO, 829 bands  
#QP = -5.568872 - -4.050642# BDA, HOMO, 1024 bands
#QP = -5.699 - -4.050642 # BDA, HOMO, 1424
#QP = -5.752 - -4.05064 # BDA, HOMO, 1630 bands


image_plane = 1.#57 

surface_z = 22.8562
carbon_z = 29.47029 

print "carbon_z: ", carbon_z
print "surface_z: ", surface_z
d = ((carbon_z - surface_z) - image_plane)
W = -13.606/(2*d)
print "d: ", d, "W: ", W, "QP: ", QP 
sigma = QP - W 

print "c6_z - surface_z: ", carbon_z - surface_z, "sigma: ", sigma
print "image plane: ", image_plane
