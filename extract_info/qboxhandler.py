import xml.sax.handler
import sys, copy
from phys import D3v
from cStringIO import StringIO
class QboxHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.inData = 0
    self.inCell = 0
    self.iteration = {}
    self.data = []
    self.temp = {}
    self.atomsbuff = []
    self.buffer = StringIO()

  def startElement(self, name, attributes):
    if name == "iteration":
      self.iteration = {}
      self.iteration["count"]  = int(attributes["count"])
      self.atomsbuff = []
    elif name == "etotal":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "ekin_ion":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "ekin":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "econf":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "eps":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "enl":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "ecoul":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "exc":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "esr":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "eself":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "pv":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "enthalpy":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "econst":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "temp_ion":
      self.buffer = StringIO()
      self.inData = 1
    elif name == "atom":
      self.name = attributes["name"]
      self.specy = attributes["species"]
      self.temp = {}
    elif name == "position":	
      self.inData = 1
      self.buffer = StringIO()
    elif name == "force":	
      self.inData = 1
      self.buffer = StringIO()
    elif name == "velocity":	
      self.inData = 1
      self.buffer = StringIO()
#    elif name == "cell":
#      self.tmp = attributes["a"].split()
#      self.tmpa = float(D3v(self.tmp[0],self.tmp[1],self.tmp[2]))
#      self.tmp = attributes["b"].split()
#      self.tmpb = float(D3v(self.tmp[0],self.tmp[1],self.tmp[2]))
#      self.tmp = attributes["c"].split()
#      self.tmpc = float(D3v(self.tmp[0],self.tmp[1],self.tmp[2]))
#      self.cell = [self.tmpa,self.tmpb,self.tmpc]
#    elif name == "unit_cell" or "cell":
    elif name == "unit_cell":
#      print name, attributes
      self.inCell = 1
      self.cellbuff = {}
    elif self.inCell == 1 and (name == "a" or name == "b" or name == "c"):
#      print name,attribute
      self.inData = 1
      self.buffer = StringIO()
    elif name == "refcell":
      (x,y,z) = attributes["a"].split()
      self.tmpa = D3v(float(x),float(y),float(z))
      (x,y,z) = attributes["b"].split()
      self.tmpb = D3v(float(x),float(y),float(z))
      (x,y,z) = attributes["c"].split()
      self.tmpc = D3v(float(x),float(y),float(z))
      self.refcell = [self.tmpa,self.tmpb,self.tmpc]
    elif name == "stress_tensor":
      self.inStress = 1
      self.strbuff = {}
    elif name.find("sigma_") != -1 and self.inStress == 1:
      self.inData = 1
      self.buffer = StringIO()

  def characters(self, data):
    if self.inData:
#        print data
        self.buffer.write(data)

  def endElement(self, name):
    if name == "etotal":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "ekin_ion":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "ekin":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "econf":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "eps":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "enl":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "ecoul":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "exc":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "esr":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "eself":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "pv":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "enthalpy":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "econst":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "temp_ion":
      self.inData = 0
      self.fbuff = float(self.buffer.getvalue())
      self.iteration[name] = self.fbuff
    elif name == "atom":
      self.inData = 0
      self.temp["name"] = self.name
      self.temp["specy"] = self.specy
      self.atomsbuff += [self.temp]
#      print self.temp
    elif name == "position":
      self.inData = 0
      (x,y,z) = self.buffer.getvalue().split()
      self.temp["position"] = D3v(float(x),float(y),float(z))
    elif name == "force":
      self.inData = 0
      (x,y,z) = self.buffer.getvalue().split()
      self.temp["force"] = D3v(float(x),float(y),float(z))
    elif name == "velocity":
      self.inData = 0
      (x,y,z) = self.buffer.getvalue().split()
      self.temp["velocify"] = D3v(float(x),float(y),float(z))
    elif self.inCell == 1 and (name == "a" or name == "b" or name == "c"):
      self.inData = 0
      (x,y,z) = self.buffer.getvalue().split()
      self.cellbuff[name] = D3v(float(x),float(y),float(z))
    elif name == "unit_cell":
      self.inCell = 0
      self.iteration["unit_cell"] = copy.deepcopy(self.cellbuff)
    elif name.find("sigma_") != -1 and self.inStress == 1:
      self.inData = 0
#      print name,self.buffer.getvalue()
      self.strbuff[name] = float(self.buffer.getvalue())
    elif name == "stress_tensor":
      self.inStress = 0
      self.iteration["stress"] = copy.deepcopy(self.strbuff)
      
    elif name == "iteration":
      self.iteration["atoms"] = copy.deepcopy(self.atomsbuff)
      self.data += [self.iteration]

class QboxPsHandler(xml.sax.handler.ContentHandler):
  def __init__(self):
    self.inData = 0
    self.buffer = StringIO()

  def startElement(self, name, attributes):
    if name == "mass":
      self.inData = 1
    elif name == "valence_charge":
      self.inData = 1

  def characters(self, data):
    if self.inData:
        self.buffer.write(data)

  def endElement(self, name):
    if name == "mass":
      self.inData = 0
      self.mass = float(self.buffer.getvalue())
      self.buffer = StringIO()
    elif name == "valence_charge":
      self.inData = 0
      self.valence = float(self.buffer.getvalue())
      self.buffer = StringIO()

