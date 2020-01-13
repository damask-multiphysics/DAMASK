#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import vtk
import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Create regular voxel grid from points in an ASCIItable.

""", version = scriptID)

parser.add_option('-m',
                  '--mode',
                  dest = 'mode',
                  metavar='string',
                  type = 'choice', choices = ['cell','point'],
                  help = 'cell-centered or point-centered coordinates')
parser.add_option('-p',
                  '--pos', '--position',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')

parser.set_defaults(mode   = 'cell',
                    pos    = 'pos',
                   )

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

  coords = [np.unique(table.get(options.pos)[:,i]) for i in range(3)]
  if options.mode == 'cell':
    coords = [0.5 * np.array([3.0 * coords[i][0] - coords[i][0 + int(len(coords[i]) > 1)]] + \
                             [coords[i][j-1] + coords[i][j] for j in range(1,len(coords[i]))] + \
                             [3.0 * coords[i][-1] - coords[i][-1 - int(len(coords[i]) > 1)]]) for i in range(3)]

  grid = np.array(list(map(len,coords)),'i')
  N = grid.prod() if options.mode == 'point' else (grid-1).prod()

# ------------------------------------------ process data ---------------------------------------

  rGrid = vtk.vtkRectilinearGrid()
  coordArray = [vtk.vtkDoubleArray(),
                vtk.vtkDoubleArray(),
                vtk.vtkDoubleArray(),
               ]

  rGrid.SetDimensions(*grid)
  for i,points in enumerate(coords):
    for point in points:
      coordArray[i].InsertNextValue(point)

  rGrid.SetXCoordinates(coordArray[0])
  rGrid.SetYCoordinates(coordArray[1])
  rGrid.SetZCoordinates(coordArray[2])


# ------------------------------------------ output result ---------------------------------------  

  if name:
    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(os.path.split(name)[0],
                                    os.path.splitext(os.path.split(name)[1])[0] +
                                    '_{}({})'.format(options.pos, options.mode) +
                                    '.' + writer.GetDefaultFileExtension()))
  else:
    writer = vtk.vtkDataSetWriter()
    writer.SetHeader('# powered by '+scriptID)
    writer.WriteToOutputStringOn()
  
  writer.SetInputData(rGrid)

  writer.Write()

  if name is None: sys.stdout.write(writer.GetOutputString())
