import math, sys, copy

ao = 0.529177249

# following constants are adpoted from Burkhard's Physics.h

eCharge    = 1.602176462e-19
mu0       = 4.0*math.pi*1e-7
c         = 299792458.0
eps0      = 1.0/(mu0*c*c)
h         = 6.62606876e-34
hBar      = h/(2.0*math.pi)
kb        = 1.3806503e-23
avogadro  = 6.02214199e23
fConst    = eCharge*eCharge/(4.0*math.pi*eps0)
me        = 9.10938188e-31
mp        = 1.67262158e-27
md        = 1.99900750083*mp
ab        = hBar*hBar/(fConst*me)
Ha        = fConst*fConst*me/hBar/hBar
Ry        = Ha/2.0
	     
e0        = -0.5*1.1676     # Ha per atom
rho0D2GCC =0.171            # rho0 D2 g/cc
rho0D2KGM3=rho0D2GCC*1000.0 # rho0 D2 kg/m^3
nn0       =rho0D2KGM3/md*ab*ab*ab # number of D atoms in ab^3
p0        =0.0
	     
GPaToAU  = 1e9/(Ha/ab/ab/ab)
AUToGPa  = 1.0/GPaToAU
MbarToAU = 100.0*GPaToAU
AUToMbar = 1.0/MbarToAU
AUToK    = Ha/kb
KToAU    = kb/Ha
K4ToAU   = 10000.0*kb/Ha
AUToeV   = Ha/eCharge
eVToAU   = 1.0/AUToeV
timeAU   = hBar/Ha
sToAU    = Ha/hBar               # == 1/tAU(sec)
e6cmPerSToAU = 1e6*0.01/ab/sToAU # 10^4*tAU/lAU
kmPerSToAU = 1e3/ab/sToAU        # 10^3*tAU/lAU
AUTokmPerS = 1.0/kmPerSToAU
AUToe6cmPerS = 1.0/e6cmPerSToAU
AUTogcc  = md/(ab*ab*ab)/1000.0
gccToAU  = 1.0/AUTogcc
eVToK    = eVToAU*AUToK
	     
mdAU     = md/me
rho0AlKGM3 = 2710.0  # Nellis paper Tab. III
rho0AlAU = rho0AlKGM3/me*ab*ab*ab

class D3v:
    X = 0
    Y = 0
    Z = 0
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def pbc(self):
        if D3v.X == 0 or D3v.Y == 0 or D3v.Z == 0:
            sys.stderr.write("Cell parameter must be set\n")
            sys.exit()

        self.x /= D3v.X
        while self.x < -0.5:
            self.x += 1
        while self.x >= 0.5:
            self.x -= 1
        self.x *= D3v.X
            
        self.y /= D3v.Y
        while self.y < -0.5:
            self.y += 1
        while self.y >= 0.5:
            self.y -= 1
        self.y *= D3v.Y

        self.z /= D3v.Z
        while self.z < -0.5:
            self.z += 1
        while self.z >= 0.5:
            self.z -= 1
        self.z *= D3v.Z


    def __add__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            return D3v(self.x + other.x, self.y + other.y, self.z + other.z)
        else:
            return D3v(self.x + other, self.y + other, self.z + other)

    __radd__ = __add__


    def __sub__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            return D3v(self.x - other.x, self.y - other.y, self.z - other.z)
        else:
            return D3v(self.x - other, self.y - other, self.z - other)

    __rsub__ = __sub__

    def __neg__(self):
        return D3v(-self.x, -self.y, -self.z)

    def __pos__(self):
        return D3v(self.x, self.y, self.z)

    def __mul__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            return D3v(self.x * other.x, self.y * other.y, self.z * other.z)
        else:
            return D3v(self.x * other, self.y * other, self.z * other)

    __rmul__ = __mul__

    def __div__(self,other):
        return D3v(self.x/other,self.y/other,self.z/other)

    def dot(self,other):
        return self.x*other.x + self.y*other.y + self.z*other.z

    def norm(self):
        return self.dot(self)

    def __abs__(self):
        return math.sqrt(self.norm())

    def __iadd__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            self.x += other.x
            self.y += other.y
            self.z += other.z
        else:
            self.x += other
            self.y += other
            self.z += other
        return self
        
    def __imul__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            self.x *= other.x
            self.y *= other.y
            self.z *= other.z
        else:
            self.x *= other
            self.y *= other
            self.z *= other
        return self

    def __isub__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            self.x -= other.x
            self.y -= other.y
            self.z -= other.z
        else:
            self.x -= other
            self.y -= other
            self.z -= other
        return self

    def __idiv__(self,other):
        if hasattr(other,"__class__") and other.__class__ is D3v:
            self.x /= other.x
            self.y /= other.y
            self.z /= other.z
        else:
            self.x /= other
            self.y /= other
            self.z /= other
        return self
        

    def __str__(self):
        return "%g %g %g" % (self.x, self.y, self.z)

    def __float__(self):
        return D3v(float(self.x), float(self.y), float(self.z))
        
    def crossprod(self,other):
        return D3v( self.y * other.z - self.z * other.y ,
                         self.z * other.x - self.x * other.z ,
                         self.x * other.y - self.y * other.x )

    def volume(self,other1,other2):
        return self.dot(other1.crossprod(other2))

    def normalized(self):
        return self/abs(self)

    def normHalfCell(self):
        return (D3v.X*D3v.X + D3v.Y*D3v.Y + D3v.Z*D3v.Z)/4

    def normQuatCell(self):
        return (D3v.X*D3v.X + D3v.Y*D3v.Y + D3v.Z*D3v.Z)/64

    def pbc_v2(self,a,b,c):
        bc = b.crossprod(c)
	v = a.dot(bc)

        x = self.dot(bc)/v
        while x < -0.5:
            x += 1.0
        while x >= 0.5:
            x -= 1.0

        ca = c.crossprod(a)
	v = b.dot(ca)
        y = self.dot(ca)/v
        while y < -0.5:
            y += 1.0
        while y >= 0.5:
            y -= 1.0

        ab = a.crossprod(b)
	v = c.dot(ab)
        z = self.dot(ab)/v
        while z < -0.5:
            z += 1.0
        while z >= 0.5:
            z -= 1.0

        self.x = x*a.x+y*b.x+z*c.x
        self.y = x*a.y+y*b.y+z*c.y
        self.z = x*a.z+y*b.z+z*c.z
#        self = copy.deepcopy(x*a+y*b+z*c)
#        sys.stderr.write(str(self)+": ")

    def pbc_v3(self,a,b,c):
        bc = b.crossprod(c)
	v = a.dot(bc)

        x = self.dot(bc)/v
        while x < 0.0:
            x += 1.0
        while x >= 1.0:
            x -= 1.0

        ca = c.crossprod(a)
	v = b.dot(ca)
        y = self.dot(ca)/v
        while y < 0.0:
            y += 1.0
        while y >= 1.0:
            y -= 1.0

        ab = a.crossprod(b)
	v = c.dot(ab)
        z = self.dot(ab)/v
        while z < 0.0:
            z += 1.0
        while z >= 1.0:
            z -= 1.0

        self.x = x*a.x+y*b.x+z*c.x
        self.y = x*a.y+y*b.y+z*c.y
        self.z = x*a.z+y*b.z+z*c.z

    def rotz(self,angle):
        cs = math.cos(angle)
        sn = math.sin(angle)
        x = self.x
        y = self.y
        self.x = cs*x-sn*y
        self.y = sn*x+cs*y

    def rotation(self,axis,angle):
        cs = math.cos(angle)
        sn = math.sin(angle)
        r=self
        axis.normalized()
        #r' = r cos(a)+n(n.r)(1-cos(a))+(r x n)sin(a) 
        #angle is clockwise
        return r*cs + axis*r.dot(axis)*(1-cs) + (r.crossprod(axis))*sn


class PeriodicCell:
    def __init__(self,a,b,c):
        self.a = a
        self.b = b
        self.c = c
        self.__ab = a.crossprod(b)
        self.__bc = b.crossprod(c)
        self.__ca = c.crossprod(a)
        self.__va = a.dot(self.__bc)
        self.__vb = b.dot(self.__ca)
        self.__vc = c.dot(self.__ab)
        self.volume = abs(self.__va)
        self.__large = 100000.5
    def Reset(self,a,b,c):
        self.a = a
        self.b = b
        self.c = c
        self.__ab = a.crossprod(b)
        self.__bc = b.crossprod(c)
        self.__ca = c.crossprod(a)
        self.__va = a.dot(self.__bc)
        self.__vb = b.dot(self.__ca)
        self.__vc = c.dot(self.__ab)
        self.volume = abs(self.__va)
        self.__large = 100000.5
    def PBC(self,r):
        x = r.dot(self.__bc)/self.__va
        x += self.__large
        x -= int(x)+0.5
        
        y = r.dot(self.__ca)/self.__vb
        y += self.__large
        y -= int(y)+0.5

        z = r.dot(self.__ab)/self.__vc
        z += self.__large
        z -= int(z)+0.5

        return D3v(x*self.a.x+y*self.b.x+z*self.c.x, x*self.a.y+y*self.b.y+z*self.c.y, x*self.a.z+y*self.b.z+z*self.c.z)
    
class PeriodicD3v(D3v):
    cell  = PeriodicCell(D3v(0.0,0.0,0.0),D3v(0.0,0.0,0.0),D3v(0.0,0.0,0.0))
    def SetCell(self,a,b,c):
        PeriodicD3v.cell.Reset(a,b,c)
    def pbc(self):
        rr = PeriodicD3v.cell.PBC(self)
        self.x = rr.x
        self.y = rr.y
        self.z = rr.z

class PeriodicPoints:
    def __init__(self,a,b,c):
        self.cell = PeriodicCell(a,b,c)
        self.r = []
    def AddPoint(self,r):
        self.r += [r]
    def AddPoints(self,r):
        self.r += r
            
class DiAtomicMolecule:
    def __init__(self,atoms,rcut):
        self.rcut = rcut
        pairs = {}
        rcut2 = rcut**2
        self.rcut2 = rcut
        for atom2 in atoms:
            for atom1 in atoms:
                if atom2 == atom1:
                    continue
                rr = atom2["position"] - atom1["position"]
                rr.pbc()
                d2 = rr.norm()
                if d2 < rcut2:
                    d = math.sqrt(d2)
                    if pairs.has_key(atom2["name"]):
                        pair = pairs[atom2["name"]]
                        if d < pair["dist"]:
                            pairs[atom2["name"]] = {"name":atom1["name"], "dist":d}
                    else:
                        pairs[atom2["name"]] = {"name":atom1["name"], "dist":d}

        self.pairlist = {}
        for atom1 in pairs.keys():
            if pairs[atom1]["dist"] > rcut:
                sys.stderr.write("too long bond: "+str(pairs[atom1]["dist"]))
                sys.exit()

            self.pairlist[atom1] = pairs[atom1]["name"]

    def getCenterOfMass(self,atoms):
        self.cmlist = []
        namelist = []
        for atom in atoms:
            if namelist.count(atom["name"]) >= 1:
                continue
            name1 = atom["name"]
            namelist += [name1]
            r1 = atom["position"]
            f1 = atom["force"]
            for atom2 in atoms:
                name2 = atom2["name"]
                if name2 == self.pairlist[name1]:
                    namelist += [name2]
                    r2 = atom2["position"]
	            f2 = atom["force"]
		    cmf = (f1+f2)*0.5
                    rr = r2 - r1
                    rr.pbc()
                    d = abs(rr)
                    if d > self.rcut:
                        sys.stderr.write("bond length too large: " + str(d))
                        sys.exit()
                    r2 = r1 + rr
                    cm = (r1 + r2)/2
                    self.cmlist += [{"name1":name1, "name2":name2, "position":cm, "length":d, "v": rr, "cmf": cmf}]
        return self.cmlist
        
