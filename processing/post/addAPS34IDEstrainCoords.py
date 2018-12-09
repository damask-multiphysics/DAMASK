#!/usr/bin/env python2.7
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
Transform X,Y,Z,F APS BeamLine 34 coordinates to x,y,z APS strain coordinates.

""", version = scriptID)

parser.add_option('-f','--frame',       dest='frame', nargs=4, type='string', metavar='string string string string',
                                        help='APS X,Y,Z coords, and depth F')
parser.set_defaults(frame = None)

(options,filenames) = parser.parse_args()

if options.frame is None:
  parser.error('no data column specified...')


datainfo = {'len':3,
            'label':[]
           }

datainfo['label']  += options.frame

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

  if not table.label_dimension(options.quaternion) == 4:
    damask.util.croak('input {} does not have dimension 4.'.format(options.quaternion))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue



  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}
  columnMissing = False

  for label in datainfo['label']:    
    key = label
    if key in table.labels(raw = True):
        active.append(label)
        column[label] = table.labels.index(key)                                                     # remember columns of requested data
    else:
      file['croak'].write('column %s not found...\n'%label)
      columnMissing = True
      
  if columnMissing: continue

# ------------------------------------------ assemble header ---------------------------------------
  table.labels_append(['%i_coord'%(i+1) for i in range(3)])                                         # extend ASCII header with new labels
  table.head_write()
  
# ------------------------------------------ process data ------------------------------------------
  theta=-0.75*np.pi
  RotMat2TSL=np.array([[1.,  0.,            0.],
                       [0.,  np.cos(theta), np.sin(theta)],                                         # Orientation to account for -135 deg
                       [0., -np.sin(theta), np.cos(theta)]])                                        # rotation for TSL convention
  vec         = np.zeros(4)
  
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for i,label in enumerate(active):
      vec[i] = table.data[column[label]]

    table.data_append(np.dot(RotMat2TSL,np.array([-vec[0], -vec[1],-vec[2]+vec[3]])))

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
