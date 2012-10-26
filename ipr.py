#!/usr/local/bin/python
import os,sys, math
from phys import *
from base64 import decodestring
from array import array
from cStringIO import StringIO
#useage: qbox2plot.py format xmlfile 1st_orbital Last_orbital plot(output)file
#
#Modify Tadashi's script to handle hexagonal cells
#We need to parse the XML output from Qbox and export:  
#-Lattice vectors
#-Atomic coordinates
#-wavefunction file(s)
#-density file(s)
#
#It would be nice to be able to export using:
#-OpenDX format
#-molmol
#-map3Dv
#-other
#
#for now, I'll export to dx (and try not to break the map3Dv related part
#norm=sum(phi**2)*(lx*ly*lz)/(nx*ny*nz)
#inverse_participation_ratio=ntot*(sum(phi**4)*(lx*ly*lz)/(nx*ny*nz))/norm
#

if len(sys.argv) != 4:
    print >>sys.stderr,"Usage: ipr.py xmlfile 1st_orbital Last_orbital"
    sys.exit()

ust = int(sys.argv.pop())
lst = int(sys.argv.pop())
file = sys.argv.pop()

totalcharge = []
f = open(file,"r",-1)

waveaddress = []

data = []
nst = 0
atoms = []
locpara = []

f.seek(0,0)

while 1:
    line = f.readline()
    if line.find("<grid_function ") != -1:
        tmp = line.split('\"')
        nx = int(tmp[3])
        ny = int(tmp[5])
        nz = int(tmp[7])
        if len(totalcharge) == 0:
            for i in range(0,nx*ny*nz):
                totalcharge += [0.0]

            index = []
            for iz in range(0,nz):
                if iz < nz/2:
                    jz = iz + nz/2
                else:
                    jz = iz - nz/2
                for iy in range(0,ny):
                    if iy < ny/2:
                        jy = iy + ny/2
                    else:
                        jy = iy - ny/2
                    for ix in range(0,nx):
                        if ix < nx/2:
                            jx = ix + nx/2
                        else:
                            jx = ix - nx/2
                        jj = jx*ny*nz+jy*nz+jz
                        index += [jj]

        sizeofwavefun = nx*ny*nz*8
        waveaddress += [f.tell()]
        break
    elif line.find("<slater_determinant ") != -1:
        print >>sys.stderr,"Start reading wavefunction"
    elif line.find("<density_matrix ") != -1:
        numberofstates = int(line.split('"')[-2])
    elif line.find("<domain a=") != -1:
        tmp = line.split('\"')
        tmp2 = tmp[1].split() 
        D3v.X = float(tmp2[0])  #this just works for rectangular lattices but is needed for PlotAtomsAndCharge?
        a_latt_vec=D3v(float(tmp2[0]),float(tmp2[1]),float(tmp2[2]))
    elif line.find("b=") != -1:
        tmp = line.split('\"')
        tmp2 = tmp[1].split()
        D3v.Y = float(tmp2[1])  #this just works for rectangular lattices but is needed for PlotAtomsAndCharge?
        b_latt_vec=D3v(float(tmp2[0]),float(tmp2[1]),float(tmp2[2]))
    elif line.find("c=") != -1:
        tmp = line.split('\"')
        tmp2 = tmp[1].split()
        D3v.Z = float(tmp2[2])  #this just works for rectangular lattices but is needed for PlotAtomsAndCharge?
        c_latt_vec=D3v(float(tmp2[0]),float(tmp2[1]),float(tmp2[2]))
        cell = PeriodicCell(a_latt_vec,b_latt_vec,c_latt_vec)
        vol = cell.volume
        r= PeriodicD3v(0.0,0.0,0.0)
        r.SetCell(a_latt_vec,b_latt_vec,c_latt_vec) 

print >> sys.stderr, "Count the number of states. Total is supposed to be",numberofstates

fact = 1.35
f.seek(waveaddress[0],0)
print >> sys.stderr,1,waveaddress[0],f.tell()
nst = 1
attheend = False
while 1:
    f.seek(sizeofwavefun*fact,1)
    nl = 0
    while 1:
        line = f.readline()
        if line.find("<grid_function ") != -1:
            waveaddress += [f.tell()]
            print >>sys.stderr,nst+1,f.tell(),nl
            nst += 1
            break
        elif line.find("</wavefunction>") != -1:
            attheend = True
            break
        elif f.tell() > sizeofwavefun*numberofstates*2 or line.find("</slater_derminant>") != -1:
            print >>sys.stderr,"Please chose smaller fact than",fact
            sys.exit()
        nl += 1
    if attheend: break


# plotting

if ust > 0:
    mxst = ust
else:
    mxst = nst+ust+1

if lst > 0:
    mnst = lst-1
else:
    mnst = nst+lst

for j in range(mnst,mxst):
    f.seek(waveaddress[j],0)
    wavefun = StringIO()
    while 1:
        line = f.readline()
        if line.find("</grid_function>") != -1: break
        wavefun.write(line)
 
    data = decodestring(wavefun.getvalue()) #base64 decoding to string
    a = array('d',[]) #assign array object with double elements
    a.fromstring(data) #string to double array
    b = []
    for aa in a:
        b += [aa**2]
 
 
    dn = 1.0/float(len(a))

    sum = 0.0
    for dd in b:
        sum += dd*dn
    print "Normalization = ",sum 

    ipr = 0.0
    for dd in b:
        ipr += dd**2*dn
    print "Inverse participation ratio = ",ipr*(nx*ny*nz/sum)

    vs = 1.0/math.sqrt(vol)
    ii = 0
    for aa in a:
        b[index[ii]] = aa*vs
        ii += 1

print >>sys.stderr,"end reading",nst
print >>sys.stderr,"outputing the total charge of states: ",mnst+1," to ",mxst

f.close()
