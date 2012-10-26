#!/usr/bin/python

import os, sys
import numpy as np

def example_function(a,b,c=17.5,d=True):
    """
          variables that do not
    """
    # Now just code in the function
    sum = a+b+c
    if d:
        print 'd = True'
    else:
        print 'd != True'

    return sum


def main():
    """
    read an array from file, compute the running average, write to disk    
    """
    # A clean way to ask for user input
    try:
        # Attempt to retrieve required input from user
        prog = sys.argv[0]
        inputFilename = sys.argv[1]
        stride = int(sys.argv[2])
    except IndexError:
        print '\nusage: '+prog+' inputFilename\n'
        sys.exit(0)

    # Execute the function defined above

    inputFile = open(inputFilename,'r')
    
    x = []
    y = []

    for line in inputFile.readlines():

       if line.split() != '#':   # looks for header

          x.append(float(line.split()[0]))
          y.append(float(line.split()[1]))


    outputFile = open('raverage.dat','w')

#    for i in np.arange(stride, len(y) - stride):
    for i in np.arange(0, len(y)):
        if i > stride or i > len(y) - stride:
            outputFile.write(str(x[i]) +' '+str(np.average(y[(i-stride):(i+stride+1)]))+'\n')
        else:
            outputFile.write(str(x[i]) +' '+str(y[i])+'\n')
    outputFile.close()

# This executes main() only if executed from shell
if __name__ == '__main__':
    main()
