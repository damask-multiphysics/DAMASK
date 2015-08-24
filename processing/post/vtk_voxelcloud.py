#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create hexahedral voxels around points in an ASCIItable.

""", version = scriptID)

parser.add_option('-p', '--positions',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'coordinate label [%default]')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'x,y,z size of voxel')
parser.add_option('-o', '--origin',
                  dest = 'origin',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'x,y,z origin of coordinate system')
parser.add_option('-g', '--geom',
                  dest = 'geom', action='store_true',
                  help = 'derive geometry from geom-file header information [%default]')
parser.set_defaults(position = 'pos',
                    origin = (0.0,0.0,0.0),
                    geom = False,
                   )

(options, filenames) = parser.parse_args()

if options.size == None and not options.geom:
  parser.error('no size specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, readonly = True)
  except: continue
  table.croak(damask.util.emph(scriptName)+(': '+name if name else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  
  if options.geom:
    info,extra_header = table.head_getGeom()
  
    table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
                 'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
                 'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
                 'homogenization:  %i'%info['homogenization'],
                 'microstructures: %i'%info['microstructures'],
                ])
  else:
    info = {}
    info['size'] = np.ones(3)
    info['grid'] = info['size'] / options.size
    info['origin'] = options.origin
  
  errors = []
  if table.label_dimension(options.position) != 3: errors.append('columns "{}" have dimension {}'.format(options.position,
                                                                                                         table.label_dimension(options.position)))
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ process data ---------------------------------------  

  table.data_readArray(options.position)
  table.data[:,0:3] *= info['size']
  table.data[:,0:3] += info['origin']

  hexPoints = np.array([[-1,-1,-1],
                        [ 1,-1,-1],
                        [ 1, 1,-1],
                        [-1, 1,-1],
                        [-1,-1, 1],
                        [ 1,-1, 1],
                        [ 1, 1, 1],
                        [-1, 1, 1],
                       ])

  halfDelta = 0.5*info['size']/info['grid']

  Points = vtk.vtkPoints()
  Hex    = vtk.vtkHexahedron()
  uGrid  = vtk.vtkUnstructuredGrid()
    
  for p in table.data:
    for i,h in enumerate(hexPoints):
      id = Points.InsertNextPoint(p+h*halfDelta)
      Hex.GetPointIds().SetId(i,id)

    uGrid.InsertNextCell(Hex.GetCellType(), Hex.GetPointIds())

  uGrid.SetPoints(Points)

# ------------------------------------------ output result ---------------------------------------  

  if name:
    (dir,filename) = os.path.split(name)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(dir,os.path.splitext(filename)[0]+'_{}'.format(options.position) \
                                        +'.'+writer.GetDefaultFileExtension()))
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(uGrid)
    else:                          writer.SetInputData(uGrid)
    writer.Write()

  else:
    writer = vtk.vtkUnstructuredGridWriter()
    writer.WriteToOutputStringOn()
    writer.SetFileTypeToASCII()
    writer.SetHeader('# powered by '+scriptID)
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(uGrid)
    else:                          writer.SetInputData(uGrid)
    writer.Write()
    sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()                                                                                     # close input ASCII table

