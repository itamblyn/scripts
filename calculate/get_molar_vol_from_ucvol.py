#! /usr/bin/python

""" Calculate molar volume from abinit ucvol in cm^3 """

import sys

# Number of atoms in cell
global NATOM
NATOM = 2

def main():
    try:
        ucvol = float(sys.argv[1])
    except IndexError:
        print '\nPlease provide ucvol in bohr^3\n'
        sys.exit(0)
    
    # Convert ucvol from bohr^3 to cm^3
    conversion = 1.4850979e-25
    ucvol *= conversion
    
    avagadro = 6.0221415E+23
    NMOL = NATOM / avagadro
    
    molar_volume = ucvol / NMOL
    
    print '\nFrom a value of ucvol = '+str(ucvol)+' cm^3'
    print '\n --->  Molar Volume = '+str(molar_volume)+' cm^3/mol\n'

if __name__ == '__main__':
    main()
