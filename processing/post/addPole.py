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
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar = 'string',
                  help = 'quaternion label')

parser.set_defaults(pole = (1.0,0.0,0.0),
                    quaternion = 'orientation',
                    polar   = False,
                   )

(options, filenames) = parser.parse_args()

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

  if not table.label_dimension(options.quaternion) == 4:
    damask.util.croak('input {} does not have dimension 4.'.format(options.quaternion))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue

  column = table.label_index(options.quaternion)

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['{}_pole_{}{}{}'.format(i+1,*options.pole) for i in range(2)])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    o = damask.Orientation(quaternion = np.array(list(map(float,table.data[column:column+4]))))

    rotatedPole = o.quaternion*pole                                                                 # rotate pole according to crystal orientation
    (x,y) = rotatedPole[0:2]/(1.+abs(pole[2]))                                                      # stereographic projection

    table.data_append([np.sqrt(x*x+y*y),np.arctan2(y,x)] if options.polar else [x,y])               # cartesian coordinates

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
