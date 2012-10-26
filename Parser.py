from phys import *
import sys, os
from types import *

class PsParser:
    def __init__(self,file):
        ps = open(file,"r")
        lines = ps.readlines()
        ps.close()

        line = lines.pop(0)
        while line.find('#') != -1:
            line = lines.pop(0)
        head = line
        [nmesh,zval,gamma,mass,lmax,psname,color,radii,trans,lloc,nquad] = head.split()
        self.nmesh = int(nmesh)
        self.zval = int(zval)
        self.gamma = float(gamma)
        self.mass = float(mass)
        self.lmax = int(lmax)
        self.psname = psname
        self.color = color
        self.radii = float(radii)
        self.trans = float(trans)
        self.lloc = int(lloc)
        self.nquad = int(nquad)
        self.pot = []
        self.rad = []
        self.wf = []
        l = 0
        while l <= self.lmax:
            data = lines[0:self.nmesh]
            del lines[0:self.nmesh]
            if len(lines) > 0:
                tmp = lines.pop(0)
                tmp.strip()
                nmesh = int(tmp)
                if nmesh != self.nmesh:
                    sys.stderr.write("mesh size inconsistent\n")
                    sys.exit()
            rtmp = []
            pstmp = []
            for elem in data:
                elem.strip()
                tmp = elem.split()
                rtmp += [float(tmp[0])]
                pstmp += [float(tmp[1])]
            self.rad = rtmp
            self.pot += [pstmp]

            data = lines[0:self.nmesh]
            del lines[0:self.nmesh]
            if len(lines) > 0:
                tmp = lines.pop(0)
                tmp.strip()
                nmesh = int(tmp)
                if nmesh != self.nmesh:
                    sys.stderr.write("mesh size inconsistent\n")
                    sys.exit()
            wftmp = []
            for elem in data:
                elem.strip()
                tmp = elem.split()
                pstmp += [float(tmp[1])]
            self.wf += [wftmp]
            l += 1


class GPParser:
    def __init__(self,filename):
        if type(filename) == StringType:
            sys.stderr.write("Reading a GP output, " + filename + ": ")
            file = open(filename,"r")
            out = file.readlines()
            file.close()
            sys.stderr.write("Done!\n")
        elif type(filename) == ListType:
            out = filename
        else:
            sys.stderr.write("Unsupported input for GPParser\n")
            sys.exit()

        energy = [] 
        atoms = []
        cell = []
        mlwf = []
        mlwfe = []
        syslines = []
        for line in out:
	    if line.find("0:") != -1:
                tmp0 = line.split("0:")
                tmp = tmp0[1].split()
            else:
                tmp = line.split()
            if line.find("%%") != -1:
                del tmp[0:1]
                energy += [tmp]
            elif line.find("##") != -1:
                pos = tmp.index("##")
                atoms += [{"name": tmp[pos+1],"position": float(D3v(tmp[pos+2],tmp[pos+3],tmp[pos+4])), "force": float(D3v(tmp[pos+5],tmp[pos+6],tmp[pos+7]))}]
            elif line.find("set_cell") != -1:
                pos = tmp.index("set_cell")
                print >>sys.stderr, "PBC for GPParser is not implemented correctly"
                sys.exit()
                D3v.X, D3v.Y, D3v.Z = float(tmp[pos+1]), float(tmp[pos+2]), float(tmp[pos+3])
            elif line.find("dimension") != -1:
                cell += [tmp[-3:]]
            elif line.find("atoms,") != -1:
                self.num_atoms = int(tmp[tmp.index("atoms,")-1])
                self.num_states = int(tmp[-2])
            elif line.find("iprint") != -1:
                if tmp[tmp.index("iprint")+1] == "=":
                    self.iprint = int(tmp[tmp.index("iprint")+2])
            elif line.find("&&") != -1:
                pos = tmp.index("&&")
                mlwf += [{"#":int(tmp[pos+1]),"position": float(D3v(tmp[pos+2],tmp[pos+3],tmp[pos+4])),"spread":float(tmp[pos+5])}]
            elif line.find("&%") != -1:
                pos = tmp.index("&%")
                mlwfe += [{"eigenval":float(tmp[pos+2]),"espread":float(tmp[pos+3])}]

            elif (len(tmp) == 10) and (tmp[-9] == "atom"):
                syslines += [{"name":tmp[-8],"specy":tmp[-7],"position":float(D3v(tmp[-6],tmp[-5],tmp[-4])),"velocity":float(D3v(tmp[-3],tmp[-2],tmp[-1]))}]
            elif (len(tmp) == 7) and (tmp[-6] == "atom"):
                syslines += [{"name":tmp[-5],"specy":tmp[-4],"position":float(D3v(tmp[-3],tmp[-2],tmp[-1])),"velocity":D3v(0.0,0.0,0.0)}]
            elif line.find("nst") != -1:
                if tmp[0] == "nst":
                    self.num_states = int(tmp[2])
        psnamelist = {}
        self.psname = {}
        
        for line in syslines:
            psnamelist[line["specy"]] = 1
            self.psname[line["name"]] = line["specy"]
        self.ps = {}
        for ps in psnamelist.keys():
            self.ps[ps] = PsParser(ps)
        self.mass = {}
        for name in self.psname.keys():
            ps = self.psname[name]
            self.mass[name] = self.ps[ps].mass

        num_atoms = len(syslines)
        #        if not self.num_atoms:
        self.num_atoms = int(num_atoms)
#        print >>sys.stderr, "no status command. number of atoms is estimated to be ",num_atoms

        self.data = []
        if len(atoms) > 0:
            while len(atoms) > 0:
                eachconf = {}
                if len(cell) > 0:
                    tmp = cell.pop(0)
                    eachconf["cell"] = float(D3v(tmp[0],tmp[1],tmp[2]))
                if len(mlwf) > 0:
                    eachconf["mlwf"] = mlwf[0:self.num_states]
                    del mlwf[0:self.num_states]
                etmp = energy[0:self.iprint]
                del energy[0:self.iprint]
                eline = etmp.pop(0)
                eline.pop(0)
                eachconf["atoms"] = atoms[0:self.num_atoms]
                del atoms[0:self.num_atoms]
                eachconf["temp"] = float(eline.pop(0))
                eachconf["etotal"] = float(eline.pop(0))
                eachconf["econst"] = float(eline.pop(0))
                eachconf["ekin_el"] = float(eline.pop(0))
                eachconf["ekin_ion"] = float(eline.pop(0))

                self.data += [eachconf]

        else:
            eachconf = {}
            # in this case you have only one configuration
            # there is no force here, of course.
            eachconf["cell"] = D3v(D3v.X,D3v.Y,D3v.Z)
            if len(mlwf) > 0:
                eachconf["mlwf"] = mlwf[-self.num_states:]

            eline = energy.pop()
            eline.pop(0)
            eachconf["atoms"] = syslines

            eachconf["temp"] = float(eline.pop(0))
            eachconf["etotal"] = float(eline.pop(0))
            eachconf["econst"] = float(eline.pop(0))
            eachconf["ekin_el"] = float(eline.pop(0))
            eachconf["ekin_ion"] = float(eline.pop(0))

            self.data += [eachconf]
            
                

class SysParser:
    def __init__(self,sysfile):
        f = open(sysfile,"r")
        self.ps = {}
        self.atoms = []
        for line in f:
            tmp = line.split()
            if line.find("set_cell") != -1:
                a = D3v(float(tmp[1]),0.0,0.0)
                b = D3v(0.0,float(tmp[2]),0.0)
                c = D3v(0.0,0.0,float(tmp[3]))
                self.cell = PeriodicCell(a,b,c)
            elif line.find("set ") != -1 and line.find("cell") != -1 and line.find("ref_cell") == -1:
                tmp.pop(0)
                a = float(D3v(tmp[1],tmp[2],tmp[3]))
                b = float(D3v(tmp[4],tmp[5],tmp[6]))
                c = float(D3v(tmp[7],tmp[8],tmp[9]))
                self.cell = PeriodicCell(a,b,c)
            elif line.find("species") != -1:
                self.ps[tmp[-2]] = tmp[-1]
            elif line.find("atom") != -1:
                name = tmp[1]
                specy = tmp[2]
                x = float(tmp[3])
                y = float(tmp[4])
                z = float(tmp[5])
                self.atoms += [{'name':name,'specy':specy,'position':D3v(x,y,z)}]
        f.close()
                
            
                
class PWParser:
    def __init__(self,pwout):
        self.converged = False
        lines = open(pwout,"r").readlines()
        nn = 0
        for line in lines:
            if line.find("PWSCF") != -1:
                nn += 1
            elif line.find('convergence') != -1 and line.find('achieved') != -1:
                self.converged = True
            elif line.find('converged') != -1 and line.find('bfgs') != -1:
                self.converged = True
        if not self.converged:
            return
 #        if nn != 2:
#            return
            
        na = 0
        ic = -1
        self.conf = []
        atoms = []
        for line in lines:
            if line.find("number") != -1 and line.find("atoms/cell") != -1:
                na = int(line.split().pop())
            elif ic > 0:
                (type, X, Y, Z) = line.split()
                atoms += [{'name':type+str(ic),'type':type,'position':float(D3v(X,Y,Z))}]
                ic += 1
                if ic > na:
                    self.conf += [atoms]
                    ic = -1
                    atoms = []
            elif line.find("ATOMIC_POSITIONS") != -1:
                ic = 1
            elif line.find("tau") != -1:
                tmp = line.split()
                type = tmp[1]
                (X,Y,Z) = tmp[-4:-1]
                num = int(line.split("=").pop(0).split(')').pop(0).split('(').pop())
                atoms += [{'name':type+str(num),'type':type,'position':a_0*float(D3v(X,Y,Z))}]
                if num == na:
                    self.conf += [atoms]
                    atoms = []
            elif line.find("energy") != -1:
                if line.find("!") != -1:
                    self.energy = float(line.split()[-2])
            elif line.find("bravais-lattice") != -1:
                if int(line.split().pop()) != 0:
                    print >>sys.stderr,"Only ibrav = 0 is implemented"
                    sys.exit()
            elif line.find('celldm(1)=') != -1:
                a_0 = float(line.split().pop(1))
            elif line.find('a(1)') != -1:
                (ax,ay,az) = line.split()[-4:-1]
            elif line.find('a(2)') != -1:
                (bx,by,bz) = line.split()[-4:-1]
            elif line.find('a(3)') != -1:
                (cx,cy,cz) = line.split()[-4:-1]
                a = float(D3v(ax,ay,az))*a_0
                b = float(D3v(bx,by,bz))*a_0
                c = float(D3v(cx,cy,cz))*a_0
                self.cell = PeriodicCell(a,b,c)
            elif line.find('a_0') != -1:
                if line.find('lattice') != -1:
                    a_0 = float(line.split()[-2])
            elif line.find('convergence') != -1 and line.find('achieved') != -1:
                self.converged = True
                
class PWParserIn:
    def __init__(self,pwout):
        lines = open(pwout,"r").readlines()
        na = 0
        ic = -1
        self.conf = []
        atoms = []
        icell = 0
        for line in lines:
            if line.find("nat") != -1:
                if line.find(','):
                    tmp = line.split(',').pop(0)
                    tmp1 = tmp.split('=').pop()
                else:
                    tmp1 = line.split('=').pop()
                na = int(tmp1)
            elif ic > 0:
                (type, X, Y, Z) = line.split()
#                print ic,type,X,Y,Z
                atoms += [{'name':type+str(ic),'type':type,'position':float(PeriodicD3v(X,Y,Z))}]
                ic += 1
                if ic > na:
                    self.conf += [copy.deepcopy(atoms)]
                    ic = -1
#                    atoms = []
            elif line.find("ATOMIC_POSITIONS") != -1:
                ic = 1

            elif line.find("ibrav") != -1:
                ibrav = int(line.split("=").pop())
                if(not (ibrav == 0 or ibrav == 5)):
                    sys.exit("Only ibrav = 0 or 5 is implemented: "+str(ibrav))
            elif icell == 3:
                (cx,cy,cz) = line.split()
                a = float(D3v(ax,ay,az))
                b = float(D3v(bx,by,bz))
                c = float(D3v(cx,cy,cz))
#                dummy = PeriodicD3v(0.0,0.0,0.0)
#                dummy.SetCell(a,b,c)
#                print str(a),str(b),str(c)
                self.cell = PeriodicCell(a,b,c)
                icell += 1
            elif icell == 2:
                (bx,by,bz) = line.split()
                icell += 1
            elif icell == 1:
                (ax,ay,az) = line.split()
                icell += 1
            elif line.find('CELL_PARAMETERS') != -1:
                icell = 1
            elif line.find('celldm(1)') != -1:
                tmp = line.split("=")
#                print 'tmp',tmp
                celldm1 = float(tmp.pop())
            elif line.find('celldm(4)') != -1:
                tmp = line.split("=")
                celldm4 = float(tmp.pop())

        if ibrav == 5:
            aa = celldm1
            alpha = math.acos(celldm4)
            xx = aa*math.sin(alpha*0.5)
            yy = 2.0/3.0*aa*math.sin(alpha*0.5)*math.cos(math.pi/6.0)
            zz = math.sqrt(aa**2-(2.0*yy)**2)
            a = D3v(xx,-yy,zz)
            b = D3v(0.0,2.0*yy,zz)
            c = D3v(-xx,-yy,zz)
            self.cell = PeriodicCell(a,b,c)
                
#    sys.exit()
class XYZParser:
    def __init__(self,file):
        lines =  open(file).readlines()
        totallines = len(lines)
        self.conf = []
        while lines:
            na = int(lines.pop(0))
            (ax,ay,az,bx,by,bz,cx,cy,cz) = lines.pop(0).split()

            a = float(D3v(ax,ay,az))
            b = float(D3v(bx,by,bz))
            c = float(D3v(cx,cy,cz))
            self.cell = PeriodicCell(a,b,c)
            atoms = []
            for ii in range(0,na):
                line = lines.pop(0)
                (na, x, y, z) = line.split()
                type = na
                name = na+str(ii)
                atoms += [{'position':float(D3v(x,y,z)),'name':name,'type':type}]
                ii += 1
            self.conf += [atoms]
