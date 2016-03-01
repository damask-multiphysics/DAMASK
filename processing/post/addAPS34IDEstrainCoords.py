#!/usr/bin/env python
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

parser.add_option('-f','--frame',       dest='frame', nargs=4, type='string', metavar='<string string string string>',
                                        help='APS X,Y,Z coords, and depth F')
parser.set_defaults(frame = None)

(options,filenames) = parser.parse_args()

if options.frame == None:
  parser.error('no data column specified...')


datainfo = {'len':3,
            'label':[]
           }

if options.frame  != None:    datainfo['label']  += options.frame

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}
  columnMissing = False

  for label in datainfo['label']:    
    key = label
    if key in table.labels:
        active.append(label)
        column[label] = table.labels.index(key)                                                     # remember columns of requested data
    else:
      file['croak'].write('column %s not found...\n'%label)
      columnMissing = True
      
  if columnMissing: continue

# ------------------------------------------ assemble header ---------------------------------------
  table.labels_append(['%i_coord'%(i+1) for i in xrange(3)])                                       # extend ASCII header with new labels
  table.head_write()
  
# ------------------------------------------ process data ------------------------------------------
  theta=-0.75*np.pi
  RotMat2TSL=np.array([[1., 0., 0.],
                       [0.,  np.cos(theta), np.sin(theta)],
                       [0., -np.sin(theta), np.cos(theta)]])                                        # Orientation Matrix to account for -135 degree rotation for TSL Convention[Adapted from Chen Zhang's code]
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
