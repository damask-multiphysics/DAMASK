#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]
  
#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [file[s]]', description = """
Produce VTK rectilinear mesh of structure data from geom description

""", version = scriptID)

parser.add_option('-m','--nodata',
                  dest   = 'data',
                  action = 'store_false',
                  help   = 'generate mesh without microstructure index data')

parser.set_defaults(data = True,
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = ['STDIN']

for name in filenames:
  if not (name == 'STDIN' or os.path.exists(name)): continue
  table = damask.ASCIItable(name = name, outname = None,
                            buffered = False, labeled = False, readonly = True)
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name != 'STDIN' else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    file['croak'](errors)
    table.close(dismiss = True)
    continue

# --- generate VTK rectilinear grid --------------------------------------------------------------------------------

  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions([x+1 for x in info['grid']])
  for i in xrange(3):
    temp = vtk.vtkDoubleArray()
    temp.SetNumberOfTuples(info['grid'][i]+1)
    for j in xrange(info['grid'][i]+1):
      temp.InsertTuple1(j,j*info['size'][i]/info['grid'][i]+info['origin'][i])
    if   i == 0: grid.SetXCoordinates(temp)
    elif i == 1: grid.SetYCoordinates(temp)
    elif i == 2: grid.SetZCoordinates(temp)

#--- read microstructure information --------------------------------------------------------------

  if options.data:
    microstructure = table.microstructure_read(info['grid'])                                        # read microstructure

    structure = vtk.vtkIntArray()
    structure.SetName('Microstructures')

    for idx in microstructure:
      structure.InsertNextValue(idx)

    grid.GetCellData().AddArray(structure)

# --- write data -----------------------------------------------------------------------------------

  if name == 'STDIN':
    writer = vtk.vtkRectilinearGridWriter()
    writer.WriteToOutputStringOn()
    writer.SetFileTypeToASCII()
    writer.SetHeader('# powered by '+scriptID)
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(grid)
    else:                          writer.SetInputData(grid)
    writer.Write()
    sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])
  else:
    (dir,filename) = os.path.split(name)
    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(dir,'mesh_'+os.path.splitext(filename)[0]
                                                   +'.'+writer.GetDefaultFileExtension()))
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(grid)
    else:                          writer.SetInputData(grid)
    writer.Write()

  table.close()
