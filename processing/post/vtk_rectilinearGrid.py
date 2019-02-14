#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,vtk
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
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

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False,
                                    labeled  = True,
                                    readonly = True,
                                   )
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  remarks = []
  errors  = []
  coordDim = table.label_dimension(options.pos)
  if not 3 >= coordDim >= 1: errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.pos))
  elif coordDim < 3:        remarks.append('appending {} dimension{} to coordinates "{}"...'.format(3-coordDim,
                                                                                                    's' if coordDim < 2 else '',
                                                                                                    options.pos))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray(options.pos)
  if table.data.shape[1] < 3:
    table.data = np.hstack((table.data,
                            np.zeros((table.data.shape[0],
                                      3-table.data.shape[1]),dtype='f')))                           # fill coords up to 3D with zeros

  coords = [np.unique(table.data[:,i]) for i in range(3)]

  if options.mode == 'cell':
    coords = [0.5 * np.array([3.0 * coords[i][0] - coords[i][0 + int(len(coords[i]) > 1)]] + \
                             [coords[i][j-1] + coords[i][j] for j in range(1,len(coords[i]))] + \
                             [3.0 * coords[i][-1] - coords[i][-1 - int(len(coords[i]) > 1)]]) for i in range(3)]

  grid = np.array(list(map(len,coords)),'i')
  N = grid.prod() if options.mode == 'point' else (grid-1).prod()

  if N != len(table.data):
    errors.append('data count {} does not match grid {}x{}x{}.'.format(N,*(grid - (options.mode == 'cell')) ))
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

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
  
  if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(rGrid)
  else:                          writer.SetInputData(rGrid)

  writer.Write()

  if name is None: sys.stdout.write(writer.GetOutputString()[:writer.GetOutputStringLength()])      # limiting of outputString is fix for vtk <7.0

  table.close()
