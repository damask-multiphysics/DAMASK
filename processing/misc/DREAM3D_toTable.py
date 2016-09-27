#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,h5py
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

  
# --- read comments --------------------------------------------------------------------------------
                    
  coords = (np.mgrid[0:grid[2], 0:grid[1], 0: grid[0]]).reshape(3, -1).T
  coords = (np.fliplr(coords)*inFile[rootDir+'/_SIMPL_GEOMETRY/SPACING'][...] \
            + inFile[rootDir+'/_SIMPL_GEOMETRY/ORIGIN'][...] \
            + inFile[rootDir+'/_SIMPL_GEOMETRY/SPACING'][...]*0.5)

  table.data = np.hstack( (coords,
                           inFile[rootDir+'/CellData/EulerAngles'][...].reshape(grid.prod(),3),
                           inFile[rootDir+'/CellData/Phases'][...].reshape(grid.prod(),1),
                           inFile[rootDir+'/CellData/Confidence Index'][...].reshape(grid.prod(),1),
                           inFile[rootDir+'/CellData/Fit'][...].reshape(grid.prod(),1),
                           inFile[rootDir+'/CellData/Image Quality'][...].reshape(grid.prod(),1)))
                    
  
  labels = ['1_pos','2_pos','3_pos',
            '1_Euler','2_Euler','3_Euler',
            'PhaseID','CI','Fit','IQ']
  try:
    table.data = np.hstack((table.data, inFile[rootDir+'/CellData/FeatureIds'][...].reshape(grid.prod(),1)))
    labels.append(['FeatureID'])
  except Exception:
    pass
    
# ------------------------------------------ assemble header ---------------------------------------
  table.labels_clear()
  table.labels_append(labels,reset = True)
  table.head_write()
  
# ------------------------------------------ finalize output ---------------------------------------
  table.data_writeArray() #(fmt='%e2.2')
  table.close()
