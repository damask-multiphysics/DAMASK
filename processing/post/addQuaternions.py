#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add Quaternions based on Crystal Frame Coordinates.

""", version = scriptID)

parser.add_option('-f','--frame',       dest='frame', nargs=4, type='string', metavar='<string string string string>',
                                        help='heading of columns containing b* vector components and three frame vectors in that order')
parser.add_option('-s','--symmetry', dest='crysym', nargs=1,type='string',metavar='<string>',
                                     help='crystal symmetry definition')                                        
parser.set_defaults(frame = None)

(options,filenames) = parser.parse_args()

if options.frame == None:
  parser.error('no data column specified...')

datainfo = {'len':4,
            'label':[]
           }

if options.frame  != None:    datainfo['label']  += options.frame

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

  table.head_read()                                                                                 # read ASCII header info

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}

  for label in datainfo['label']:
    key = '1_'+label if datainfo['len'] > 1 else label                                              # non-special labels have to start with '1_'
    if key in table.labels:
        active.append(label)
        column[label] = table.labels.index(key)                                                     # remember columns of requested data
    else:
      damask.util.croak('column %s not found...'%label)

# ------------------------------------------ assemble header ---------------------------------------
  
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['Q_%i'%(i+1) for i in xrange(4)])                              # extend ASCII header with new labels [1 real, 3 imaginary components]
  table.head_write()
  
# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    vec    = np.zeros([4,3])
    for i,label in enumerate(active):
      vec[i,:] = np.array(table.data[column[label]:
                                     column[label]+3])
      
    if sys.argv[1:][6]=='hexagonal':                                                      # Ensure Input matrix is orthogonal
      M=np.dot(vec[0,:],vec[2,:])
      vec[1,:]=vec[1,:]/np.linalg.norm(vec[1,:])
      vec[2,:]=M*(vec[0,:]/np.linalg.norm(vec[0,:]))
      vec[3,:]=vec[3,:]/np.linalg.norm(vec[3,:])
    else:
      vec[1,:]=vec[1,:]/np.linalg.norm(vec[1,:])
      vec[2,:]=vec[2,:]/np.linalg.norm(vec[2,:])
      vec[3,:]=vec[3,:]/np.linalg.norm(vec[3,:])
     
    
    Ori=damask.Orientation(matrix=vec[1:,:],symmetry=sys.argv[1:][6])
   
    table.data_append(np.asarray(Ori.asQuaternion()))
    
    
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.close()                                                                                     # close ASCII tables
