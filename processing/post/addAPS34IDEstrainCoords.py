#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Transform X,Y,Z,F APS BeamLine 34 coordinates to x,y,z APS strain coordinates.

""", version = scriptID)

parser.add_option('-f','--frame',dest='frame', nargs=3, metavar='string string string',
                                 help='APS X,Y,Z coords')
parser.add_option('--depth',     dest='depth',          metavar='string',
                                 help='depth')

(options,filenames) = parser.parse_args()

if options.frame is None:
  parser.error('frame not specified')
if options.depth is None:
  parser.error('depth not specified')

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks -----------------------------------------
  errors  = []
  if table.label_dimension(options.frame) != 3:
    errors.append('input {} does not have dimension 3.'.format(options.frame))
  if table.label_dimension(options.depth) != 1:
    errors.append('input {} does not have dimension 1.'.format(options.depth))
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# ------------------------------------------ assemble header ---------------------------------------
  table.labels_append(['%i_coord'%(i+1) for i in range(3)])                                         # extend ASCII header with new labels
  table.head_write()
  
# ------------------------------------------ process data ------------------------------------------
  theta=-0.75*np.pi
  RotMat2TSL=np.array([[1.,  0.,            0.],
                       [0.,  np.cos(theta), np.sin(theta)],                                         # Orientation to account for -135 deg
                       [0., -np.sin(theta), np.cos(theta)]])                                        # rotation for TSL convention
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    coord = list(map(float,table.data[table.label_index(options.frame):table.label_index(options.frame)+3]))
    depth = float(table.data[table.label_index(options.depth)])

    table.data_append(np.dot(RotMat2TSL,np.array([-coord[0],-coord[1],-coord[2]+depth])))

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
