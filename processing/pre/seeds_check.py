#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,vtk
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])
     
#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog [seedsfile[s]]', description = """
Produce VTK point mesh from seeds file

""", version = scriptID)

parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'x,y,z size of hexahedral box [1.0 along largest grid point number]')
parser.add_option('-p','--position',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'column label for coordinates [%default]')

parser.set_defaults(size = [0.0,0.0,0.0],
                    position = 'pos',
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, readonly = True)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  remarks = []
  errors = []

  if     np.any(info['grid'] < 1):   remarks.append('invalid grid a b c.')
  if     np.any(info['size'] <= 0.0) \
     and np.all(info['grid'] < 1):   errors.append('invalid size x y z.')
  else:
    for i in xrange(3):
      if info['size'][i] <= 0.0:                                                                    # any invalid size?
        info['size'][i] = float(info['grid'][i])/max(info['grid'])                                  # normalize to grid
        remarks.append('rescaling size {} to {}...'.format({0:'x',1:'y',2:'z'}[i],info['size'][i]))
  if table.label_dimension(options.position) != 3: errors.append('columns "{}" have dimension {}'.format(options.position,
                                                                                                         table.label_dimension(options.position)))
  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

  labels = ['{dim}_{label}'.format(dim = 1+i,label = options.position) for i in xrange(3)]
  hasGrains = table.label_index('microstructure') != -1
  labels += ['microstructure'] if hasGrains else []

  table.data_readArray(labels)                                                                      # read ASCIItable columns

  coords = table.data[:,:3]*info['size']                                                            # assign coordinates (rescaled to box size)
  grain = table.data[:,3].astype('i') if hasGrains else 1+np.arange(len(coords),dtype='i')          # assign grains
 
# --- generate grid --------------------------------------------------------------------------------

  grid = vtk.vtkUnstructuredGrid()
  pts = vtk.vtkPoints()

# --- process microstructure information -----------------------------------------------------------

  IDs = vtk.vtkIntArray()
  IDs.SetNumberOfComponents(1)
  IDs.SetName("GrainID")

  for i,item in enumerate(coords):
    IDs.InsertNextValue(grain[i])
    pid = pts.InsertNextPoint(item[0:3])
    pointIds = vtk.vtkIdList()
    pointIds.InsertId(0, pid)
    grid.InsertNextCell(1, pointIds)

  grid.SetPoints(pts)
  grid.GetCellData().AddArray(IDs)

#--- write data -----------------------------------------------------------------------------------
  if name:
    writer = vtk.vtkXMLRectilinearGridWriter()
    (directory,filename) = os.path.split(name)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(directory,os.path.splitext(filename)[0]
                                              +'.'+writer.GetDefaultFileExtension()))
  else:
    writer = vtk.vtkDataSetWriter()
    writer.WriteToOutputStringOn()
    writer.SetHeader('# powered by '+scriptID)
  
  if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(grid)
  else:                          writer.SetInputData(grid)
  writer.Write()
  if name == None:  sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()

