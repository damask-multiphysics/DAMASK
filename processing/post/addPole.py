#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add x,y coordinates of stereographic projection of given direction (pole) in crystal frame.

""", version = scriptID)

parser.add_option('-p', '--pole',
                  dest = 'pole',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'crystal frame direction for pole figure [%default]')
parser.add_option('--polar',
                  dest = 'polar',
                  action = 'store_true',
                  help = 'output polar coordinates r,phi [%default]')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'Euler angles label')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest = 'matrix',
                  type = 'string', metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option('-a',
                  dest = 'a',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame a vector label')
parser.add_option('-b',
                  dest = 'b',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame b vector label')
parser.add_option('-c',
                  dest = 'c',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar = 'string',
                  help = 'quaternion label')

parser.set_defaults(pole = (1.0,0.0,0.0),
                    degrees = False,
                    polar   = False,
                   )

(options, filenames) = parser.parse_args()

input = [options.eulers     != None,
         options.a          != None and \
         options.b          != None and \
         options.c          != None,
         options.matrix     != None,
         options.quaternion != None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format.')

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested
toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
pole = np.array(options.pole)
pole /= np.linalg.norm(pole)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  
  if not np.all(table.label_dimension(label) == dim):  errors.append('input {} has wrong dimension {}.'.format(label,dim))
  else:  column = table.label_index(label)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['{}_pole_{}{}{}'.format(i+1,*options.pole) for i in xrange(2)])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if inputtype == 'eulers':
      o = damask.Orientation(Eulers = np.array(map(float,table.data[column:column+3]))*toRadians)
    elif inputtype == 'matrix':
      o = damask.Orientation(matrix = np.array(map(float,table.data[column:column+9])).reshape(3,3).transpose())
    elif inputtype == 'frame':
      o = damask.Orientation(matrix = np.array(map(float,table.data[column[0]:column[0]+3] + \
                                                         table.data[column[1]:column[1]+3] + \
                                                         table.data[column[2]:column[2]+3])).reshape(3,3))
    elif inputtype == 'quaternion':
      o = damask.Orientation(quaternion = np.array(map(float,table.data[column:column+4])))

    rotatedPole = o.quaternion*pole                                                                 # rotate pole according to crystal orientation
    (x,y) = rotatedPole[0:2]/(1.+abs(pole[2]))                                                      # stereographic projection

    table.data_append([np.sqrt(x*x+y*y),np.arctan2(y,x)] if options.polar else [x,y])                                                                      # cartesian coordinates

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
