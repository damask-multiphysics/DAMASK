#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,string,h5py
import numpy as np
from optparse import OptionParser
import damask

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Cauchy stress based on given column(s) of
deformation gradient and first Piola--Kirchhoff stress.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-f','--defgrad',     dest='defgrad', \
                                        help='heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',      dest='stress', \
                                        help='heading of columns containing first Piola--Kirchhoff stress [%default]')
parser.add_option('-o','--output',      dest='output', \
                                        help='group containing requested data [%default]')
parser.set_defaults(defgrad = 'f')
parser.set_defaults(stress  = 'p')
parser.set_defaults(output  = 'crystallite')

(options,filenames) = parser.parse_args()

if options.defgrad is None or options.stress is None or options.output is None:
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
