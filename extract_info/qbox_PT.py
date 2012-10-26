#!/usr/local/bin/python

import xml.sax
import qboxhandler
import sys
import os

parser = xml.sax.make_parser()

argv = []
argv = sys.argv
del argv[0]

map_s = {}
map_s["deuterium"] = "H"
map_s["hydrogen_solid"] = "H_S"
map_s["hydrogen_liquid"] = "H_L"
map_s["lithium_solid"] = "Li_S"
map_s["lithium_liquid"] = "Li_L"
map_s["boron_84a"] = "B"
map_s["boron_84b"] = "B"
map_s["boron_84i"] = "B"
map_s["boron_10"] = "B"
map_s["boron_10p"] = "B"
map_s["boron_10i"] = "B"
map_s["boron_12"] = "B"
map_s["boron_12p"] = "B"
map_s["boron_28"] = "B"
map_s["boron_28p"] = "B"
map_s["boron_28i"] = "B"
map_s["boron"] = "B"
map_s["aluminum"] = "Al"
map_s["silicon"] = "Si"

for file in argv:
  if os.path.isfile(file):
    if os.path.isfile(file+".dat"):
      continue
    cmd = "tail -1 "+file
    lasts = os.popen(cmd,"r",-1)
    lst = lasts.readline()
    if lst.find("</qbox:simulation>") == -1:
      print >> sys.stderr,file,"is not complete!"
      continue
    ofd = open(file+".dat","w")
    ofxyz = open(file+".xyz","w")
    ofpos = open(file+".pos","w")
    offxyz = open(file+".fxyz","w")
    ofvxyz = open(file+".vxyz","w")
    print >> ofd, "# 1 temp, 2 s_xx, 3 s_yy, 4 s_zz, 5 s_xy, 6 s_yz, 7 s_xz, 8-16 a:b:c, 17 volume [au], 18 etotal, 19 ekin_ion, 20 econst, 21 pv, 22 enthalpy "
    sys.stderr.write(file + " is processed.\n")
    handler = qboxhandler.QboxHandler()
    parser.setContentHandler(handler)
    parser.parse(file)
    print >>sys.stderr, "Parsing finished"

    for data in handler.data:
#      print data["atoms"]
      print >> ofpos, len(data["atoms"])
      print >> ofxyz, len(data["atoms"])
      print >> offxyz, len(data["atoms"])
      print >> ofvxyz, len(data["atoms"])
      a = data["unit_cell"]["a"]
      b = data["unit_cell"]["b"]
      c = data["unit_cell"]["c"]
      volume = a.volume(b,c)/float(len(data["atoms"]))
      cell = str(a.x) + " " + str(a.y) + " " + str(a.z) + " " + str(b.x) + " " + str(b.y) + " " + str(b.z) + " " + str(c.x) + " " + str(c.y) + " " + str(c.z)
      print >> ofxyz, cell
      print >> ofpos, cell
      print >> offxyz, cell
      print >> ofvxyz, cell
      stress = data["stress"]
      print >> ofd, data["temp_ion"], stress["sigma_xx"], stress["sigma_yy"], stress["sigma_zz"], stress["sigma_xy"], stress["sigma_yz"], stress["sigma_xz"], cell, volume, data["etotal"], data["ekin_ion"], data["econst"], data["pv"], data["enthalpy"]
      for atom in data["atoms"]:
        print >> ofpos, atom["name"],atom["specy"], atom["position"].x, atom["position"].y, atom["position"].z
        print >> ofxyz, map_s[atom["specy"]], atom["position"].x, atom["position"].y, atom["position"].z
        print >> offxyz, map_s[atom["specy"]], atom["force"].x, atom["force"].y, atom["force"].z
        print >> ofvxyz, map_s[atom["specy"]], atom["velocity"].x, atom["velocity"].y, atom["velocity"].z

  else:
    sys.stderr.write(file + " has not been found.\n")
    sys.exit()
  ofd.close()
  ofpos.close()
  ofxyz.close()
  offxyz.close()
  ofvxyz.close()

