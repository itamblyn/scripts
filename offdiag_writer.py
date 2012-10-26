#!/usr/bin/env python

import numpy

band_index_min = 16 
band_index_max = 22

#j = 21

sigma_matrix = 21


noffdiag = int((((band_index_max+1) - band_index_min)**2 - ((band_index_max+1) - band_index_min))/2.)

print "number_offdiag", noffdiag
print "begin offdiag"

for i in numpy.arange(band_index_min,band_index_max+1):

    for j in numpy.arange(i, band_index_max+1):
        if (i != j):

#            sigma_matrix = j
            print " ",i,j,sigma_matrix


print "end"
#number_offdiag   noffdiag
#begin offdiag
#   n_1   m_1   l_1
#   n_2   m_2   l_2
#   ...
#   n_noffdiag   m_noffdiag   l_noffdiag
#end

