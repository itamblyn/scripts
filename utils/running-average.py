#!/usr/bin/python

import sys

# usage: running-average.py filename column average_over

filename = sys.argv[1]
column = int(sys.argv[2]) - 1
n = int(sys.argv[3])

inputFile = open (sys.argv[1],'r')

data_array = []

#read line into array 
for line in inputFile.readlines():

    # add a new sublist
    data_array.append([])

    # loop over the elemets, split by whitespace
    for i in line.split():
    
        # convert to integer and append to the last
        # element of the list
        data_array[-1].append(float(i))

inputFile.close()

outputFile = open(filename + '.running', 'w')

t = n - 1

while t < len(data_array):

   i = t - (n - 1)

   average_count = float(0)

   while i <= t:

     average_count += data_array[i][column]
     i += 1

   outputFile.write(repr(t) + ' ' + repr(average_count/float(n)) + '\n')

   t += 1
 
outputFile.close()
