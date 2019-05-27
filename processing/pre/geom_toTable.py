#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser
from io import StringIO

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile(s)]', description = """
Translate geom description into ASCIItable containing position and microstructure.

""", version = scriptID)


(options, filenames) = parser.parse_args()

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  damask.util.croak(geom)  
  microstructure = geom.get_microstructure().flatten('F')
  grid           = geom.get_grid()
  size           = geom.get_size()
  
  for i,line in enumerate(geom.get_comments()):
    if line.lower().strip().startswith('origin'):
      origin= np.array([float(line.split()[j]) for j in [2,4,6]])                                    # assume correct order (x,y,z)

#--- generate grid --------------------------------------------------------------------------------

  x = (0.5 + np.arange(grid[0],dtype=float))/grid[0]*size[0]+origin[0]
  y = (0.5 + np.arange(grid[1],dtype=float))/grid[1]*size[1]+origin[1]
  z = (0.5 + np.arange(grid[2],dtype=float))/grid[2]*size[2]+origin[2]

  xx = np.tile(          x,        grid[1]* grid[2])
  yy = np.tile(np.repeat(y,grid[0]        ),grid[2])
  zz =         np.repeat(z,grid[0]*grid[1])

# ------------------------------------------ finalize output ---------------------------------------

  table = damask.ASCIItable(outname = os.path.splitext(name)[0]+'.txt' if name else name)
  table.info_append([scriptID + '\t' + ' '.join(sys.argv[1:])] + geom.get_comments())
  table.labels_append(['{}_{}'.format(1+i,'pos') for i in range(3)]+['microstructure'])
  table.head_write()
  table.output_flush()
  table.data = np.squeeze(np.dstack((xx,yy,zz,microstructure)),axis=0)
  table.data_writeArray()
  table.close()
