#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,h5py,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile[s]]', description = """
Convert DREAM3D file to ASCIItable

""", version = scriptID)

(options, filenames) = parser.parse_args()

rootDir ='DataContainers/ImageDataContainer'

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: parser.error('no input file specified.')

for name in filenames:
  try:
    table = damask.ASCIItable(outname = os.path.splitext(name)[0]+'.txt',
                              buffered = False
                              )
  except: continue
  damask.util.report(scriptName,name)

  inFile = h5py.File(name, 'r')
      
  grid = inFile[rootDir+'/_SIMPL_GEOMETRY/DIMENSIONS'][...]
  size = grid * inFile[rootDir+'/_SIMPL_GEOMETRY/SPACING'][...]
  origin = inFile[rootDir+'/_SIMPL_GEOMETRY/ORIGIN'][...]

# --- read comments --------------------------------------------------------------------------------
  dat = np.hstack( (inFile[rootDir+'/CellData/EulerAngles'][...].reshape(grid.prod(),3),
                    inFile[rootDir+'/CellData/Phases'][...].reshape(grid.prod(),1),
                    inFile[rootDir+'/CellData/Confidence Index'][...].reshape(grid.prod(),1),
                    inFile[rootDir+'/CellData/Fit'][...].reshape(grid.prod(),1),
                    inFile[rootDir+'/CellData/Image Quality'][...].reshape(grid.prod(),1)))

  print dat.shape
  sys.exit()
  table.labels_clear()
  table.labels_append(['1_Euler','2_Euler','3_Euler',
                       '1_pos','2_pos',
                       'IQ','CI','PhaseID','Intensity','Fit',
                      ],                                                                            # OIM Analysis 7.2 Manual, p 403 (of 517)
                      reset = True)

# ------------------------------------------ assemble header ---------------------------------------

  table.head_write()

#--- write remainder of data file ------------------------------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():
    outputAlive = table.data_write()

# ------------------------------------------ finalize output ---------------------------------------

  table.close()
