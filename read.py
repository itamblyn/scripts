#!/usr/local/bin/python


import sys, os
sys.path += ["/g/g91/ogitsu/bin"]
from phys import *

confs = []
atoms = []
for file in sys.arg[1:]:
    for line in open(file):
        tmp = line.split()
        if len(tmp) == 1:
            numatom = int(tmp[0])
            if len(atoms) > 0:
                confs += [cell,atoms]
            atoms = []
        elif len(tmp) == 9:
            (ax,ay,az,bx,by,bz,cx,cy,cz) = tmp
            aa = D3v(float(ax),float(ay),float(az))
            bb = D3v(float(bx),float(by),float(bz))
            cc = D3v(float(cy),float(cy),float(cz))
            cell = PeriodicCell(aa,bb,cc)
        elif len(tmp) == 4:
            (name,x,y,z) = tmp
            r = D3v(float(x),float(y),float(z))
            atoms += [{'element':name,'pos':r}]


for conf in confs:
    cell = conf.pop(0)
    for atoms in conf:
        for atom in atoms:
            print atom['element'],atom['pos']
