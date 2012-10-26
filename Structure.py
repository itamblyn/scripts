from phys import *
def Distances_orig(cell,atoms,atoma,atomb):
    distances = []
    for atom1 in atoms:
	for atom2 in atoms:
            if atom1['name'] == atom2['name']:
                continue
            elif atom1['type'] == atoma and atom2['type'] == atomb:
                rr = atom2['position'] - atom1['position']
                rr2 = cell.PBC(rr)
                rv= abs(rr2)
                #if the distance is different from zero, add it to the distances list
                dname=atom1['name']+'-'+atom2['name']
                dtype=atoma+'-'+atomb
                distances += [{'dname':dname,'dtype':dtype,'dist':rv,'atoma':atoma,'atomb':atomb}]
    return distances

def Distances(cell,atoms,atoma,atomb):
    distances = []
    for atom1 in atoms:
	for atom2 in atoms:
            if atom1['name'] == atom2['name']:
                break
            elif atom1['type'] == atoma and atom2['type'] == atomb:
                rr = atom2['position'] - atom1['position']
                rr2 = cell.PBC(rr)
                rv= abs(rr2)
                #if the distance is different from zero, add it to the distances list
                dname1=atom1['name']+'-'+atom2['name']
                dtype1=atoma+'-'+atomb
                distances += [{'dname':dname1,'dtype':dtype1,'dist':rv,'atoma':atoma,'atomb':atomb}]
            elif atom1['type'] == atomb and atom2['type'] == atoma:
                rr = atom2['position'] - atom1['position']
                rr2 = cell.PBC(rr)
                rv= abs(rr2)
                dname2=atom2['name']+'-'+atom1['name']
                dtype2=atomb+'-'+atoma
                distances += [{'dname':dname2,'dtype':dtype2,'dist':rv,'atoma':atomb,'atomb':atoma}]
    return distances
