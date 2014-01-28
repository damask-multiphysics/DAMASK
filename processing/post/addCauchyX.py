#!/usr/bin/env python

import os,re,sys,math,string,h5py
import numpy as np
import damask
from optparse import OptionParser, Option

# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Cauchy stress based on given column(s) of
deformation gradient and first Piola--Kirchhoff stress.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-f','--defgrad',     dest='defgrad', type='string', \
                                        help='heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',      dest='stress', type='string', \
                                        help='heading of columns containing first Piola--Kirchhoff stress [%default]')
parser.add_option('-o','--output',      dest='output', type='string', \
                                        help='group containing requested data [%default]')
parser.set_defaults(defgrad = 'f')
parser.set_defaults(stress  = 'p')
parser.set_defaults(output  = 'crystallite')

(options,filenames) = parser.parse_args()

if options.defgrad == None or options.stress == None or options.output == None:
  parser.error('missing data column...')


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'file':h5py.File(name,"a")})

# ------------------------------------------ loop over input files ------------------------------------  

for myFile in files:
  print(myFile['name'])

# ------------------------------------------ loop over increments -------------------------------------
  for inc in myFile['file']['increments'].keys():
    print("Current Increment: "+inc)
    for instance in myFile['file']['increments/'+inc+'/'+options.output].keys():
      dsets = myFile['file']['increments/'+inc+'/'+options.output+'/'+instance].keys()
      if (options.defgrad in dsets and options.stress in dsets):
        defgrad = myFile['file']['increments/'+inc+'/'+options.output+'/'+instance+'/'+options.defgrad]
        stress = myFile['file']['increments/'+inc+'/'+options.output+'/'+instance+'/'+options.stress]
        cauchy=np.zeros(np.shape(stress),'f')
        for p in range(stress.shape[0]):
          cauchy[p,...] = 1.0/np.linalg.det(defgrad[p,...])*np.dot(stress[p,...],defgrad[p,...].T)  # [Cauchy] = (1/det(F)) * [P].[F_transpose]
        cauchyFile = myFile['file']['increments/'+inc+'/'+options.output+'/'+instance].create_dataset('cauchy', data=cauchy)
        cauchyFile.attrs['units'] = 'Pa'

  

