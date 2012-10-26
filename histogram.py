#! /usr/bin/env python

import numpy, matplotlib, pylab, sys

inputFile = open ('data.dat','r')

index = 0 

data_array = []
for line in inputFile.readlines():
#   for i in line.split():
#       BAND_array.append(float(i))
    data_array.append(float(line.split()[index]))

inputFile.close()

print 'input read complete'

#pylab.hist(BAND_array, 250, 1)
pylab.hist(data_array)
#pylab.axis([-4,2.0,0,1.6])

#pylab.savefig(sys.argv[1])
pylab.show()
