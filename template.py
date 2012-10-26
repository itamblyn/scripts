#!/usr/bin/env python

import os, sys

def example_function(a,b,c=17.5,d=True):
    """
    This is an example function
    a: A number that must be handed to the function
    b: A number that must be handed to the function
    c: A value that can be handed to the function but if not
       has a default value of 17.5
    d: A value that can be handed to the function but if not
       has a default value of True
    NOTE: Variables with default values must be given after 
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
    This is the main function.
    """
    # A clean way to ask for user input
    try:
        # Attempt to retrieve required input from user
        prog = sys.argv[0]
        a = float(sys.argv[1])
        b = float(sys.argv[2])
    except IndexError:
        # Tell the user what they need to give
        print '\nusage: '+prog+' a b    (where a & b are numbers)\n'
        # Exit the program cleanly
        sys.exit(0)

    # Execute the function defined above
    sum = example_function(a,b)

    print 'sum =',sum


# This executes main() only if executed from shell
if __name__ == '__main__':
    main()
