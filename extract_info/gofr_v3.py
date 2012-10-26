#!/usr/local/bin/python

import sys,math,copy,string
import Structure
from phys import * 
from Parser import *

#read input file. use the sys format of qbox.
file = sys.argv.pop()
datatype = sys.argv.pop()
if datatype == "-sys":
	conf = SysParser(file)
	atoms = conf.atoms
	cell = conf.cell
elif datatype =="-pwscf":
	pwout = PWParser(file)
	atoms = pwout.conf.pop()
	cell = pwout.cell
elif datatype =="-pwscfin":
	pwin = PWParserIn(file)
	atoms = pwin.conf.pop()
	cell = pwin.cell
elif datatype == "-xyz":
	xyzout = XYZParser(file)
	atoms = xyzout.conf.pop()
	cell = xyzout.cell
else:
	print >> sys.stderr, "Wronge data type"
	sys.exit()


#read the labels of the atoms for which the distance will be calculated
#the idea is to call this program as many times as there are types of bonds
#for example, if we want the bonds Si-Si, Si-N, N-N, we will call the
#program 3 times.
atomb = sys.argv.pop()
atoma = sys.argv.pop()

#sys file syntax
#atom uniquelabel species x y z fx fy fz

#declare the atoms and the distances lists

suspicious_distances = []
pair_correlation_func = []
natom=0
natoma=0
natomb=0

#fill the list with the atom's labels and their positions under the
#condition that they match the type a or b.

#shave the numbers of the name to get the type.
#we are interested in atoms being of the type atoma or atomb.
#we can use this to find the atoma or atomb substring inside the name
#string to identify the type of atom.
for atom in atoms:
	name = atom['name']
	if name.count(atoma) == 1: 
		type = atoma
		natoma=natoma+1
	if name.count(atomb) == 1:
		type = atomb
		natomb=natomb+1
	if (name.count(atomb) == 0 and name.count(atoma) == 0):
		type = "other"
	atom['type'] = type



distances = Structure.Distances(cell,atoms,atoma,atomb)

#we can now evaluate the volume	
vol=cell.volume

#loop over u for the pair-correlation function

bc = cell.b.crossprod(cell.c)
ca = cell.c.crossprod(cell.a)
ab = cell.a.crossprod(cell.b)

radius=min(abs(cell.a.dot(bc)/abs(bc)),abs(cell.b.dot(ca)/abs(ca)),abs(cell.c.dot(ab)/abs(ab)))/2



#step is on the order of 0.1 bohr to differentiate between types of bonds
#for example the C-C bond is approximatively:
#1.53 A : 2.89 bohr single
#1.40 A : 2.65 bohr double
#1.20 A : 2.27 bohr triple

#nbin = 100
#step = radius/float(nbin)
step=0.02

# initialize array, and weight

# for the case atoma == atomb
c1 =2.0/((natoma-1)*float(natoma)*4.0*math.pi*step/vol)

# for the case atoma != aomb
c2 = 1.0/(float(natoma*natomb)*4.0*math.pi*step/vol)	
u = 0.0
weight = []
while u <= radius:
	u2 = u+step/2.0
	pair_correlation_func += [[u2,0.0]]
	uu = 1.0/u2**2
	if atoma==atomb: 
		# division by natoma-1 instead of natoma and extra division by 2 in this case. see note.
		wgt=c1*uu
	else: 
		wgt=c2*uu
	weight += [wgt]
	u += step

natom=len(atoms)
for dst in distances:
	rv = dst['dist']
	if rv > radius:
		continue
        if rv  < 1.0: print >>sys.stderr, dst
	# count the number of atom2 in the shell around atom1
	nn = int(rv/step)
	
	pair_correlation_func[nn][1] += weight[nn]

#note that if atoma and atomb are of the same type, the distances are accounted for twice. 
#this is taken care of by an extra division by 2 in the pair-correlation function calculation (and histograms).

#the formula for the pair-correlation is g(r)=[counts for shell r,r+dr]/vol of said shell * [total volume of box/(natom1*natom2)]
#hence, if the system has no correlation, i.e. its bond density is the average bond density (total number of bonds)/volume of box
#then g(r) -> 1. the total number of bonds (in the loose sense) is natom1*natom2.

if atoma==atomb:
	print '# coordination number for '+atoma
	sum = 0.0
	for pc in pair_correlation_func:
		sum += pc[1]*4*math.pi*((pc[0]+step*0.5)**2*step*float(natoma)/vol)
		print pc[0],pc[1],sum
else:
	sum = 0.0
	print '# coordination number for '+atoma+' as seen from '+atomb
	for pc in pair_correlation_func:
		sum += pc[1]*4*math.pi*(pc[0]+step*0.5)**2*step*float(natoma)/vol
		print pc[0],pc[1],sum
	sum = 0.0
	print '# coordination number for '+atomb+' as seen from '+atoma
	for pc in pair_correlation_func:
		sum += pc[1]*4*math.pi*(pc[0]+step*0.5)**2*step*float(natomb)/vol
		print pc[0],pc[1],sum
