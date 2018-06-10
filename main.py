from __future__ import print_function, division
from six.moves import range

import sys
from optparse import OptionParser
if sys.version_info[0] == 3:
    import subprocess as commands
else:
    import commands
import json
from twoeq import *

class Pygodunov:

  def __init__(self, jsonfile):
      self.jsonfile = jsonfile
      self.jsonfilewithextension = jsonfile + ".json"
      parser = OptionParser(usage="usage: %prog [options] filename",version="%prog 1.0")
      parser.add_option('--json',dest="jsonfilename",default=self.jsonfilewithextension,type="string", action="store",help="Input json filename")
      opt, remainder = parser.parse_args()
      with open(opt.jsonfilename) as json_file:
          self.input = json.load(json_file)


  def run_all_simulations(self):
    self.prefix = self.input["prefix"]
    self.peraff = self.input["peraff"]
    self.version = self.input["version"]
    self.Lx = self.input["Lx"]
    self.N = self.input["N"]
    self.t0 = self.input["t0"]
    self.tf = self.input["tf"]
    self.CFL = self.input["CFL"]
    self.gamma = self.input["gamma"]
    self.diffD = self.input["diffD"]
    self.diffL = self.input["diffL"]
    self.onde = self.input["onde"]
    self.rhoeL = self.input["rhoeL"]
    self.rhohL = self.input["rhohL"]
    self.peL = self.input["peL"]
    self.pL = self.input["pL"]
    self.vhL = self.input["vhL"]
    self.rhoeR = self.input["rhoeR"]
    self.rhohR = self.input["rhohR"]
    self.peR = self.input["peR"]
    self.pR = self.input["pR"]
    self.vhR = self.input["vhR"]
    if self.input["splitting"] == "True":
        self.splitting = True
        if self.input["diff_substeps"] == "True":
            self.diff_substeps = True
        else:
            self.diff_substeps = False
    else:
        self.splitting = False
        self.diff_substeps = False
    self.scheme_conv = self.input['scheme_conv']
    self.scheme_nc = self.input['scheme_nc']

    self.liste_plot = self.input.get("liste_plot", None)
    self.dir = self.prefix+"sol_" + self.jsonfile + "/"
    menage=commands.getoutput('rm -rf '+self.dir)
    createdir=commands.getoutput('mkdir -p '+self.dir)
    system = twoeq(self)
    print(system)
    system.run()

if __name__ == "__main__":
    import fnmatch
    import os

    liste  = []
    for file in sorted(os.listdir('.')):
        if fnmatch.fnmatch(file, '*.json'):
            fileName,fileExtension = os.path.splitext(file)
            liste.append((fileName))

    print("Please, choose a json test file: ")
    for k, f in enumerate(liste):
        print("{0:2d} -> {1:s}".format(k, f))
    numFile = int(input())
    if numFile<0 or numFile>len(liste):
        print("Error in the choice")
    else:
        c = Pygodunov(jsonfile = liste[numFile])
        c.run_all_simulations()