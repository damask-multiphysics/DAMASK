#!/usr/bin/env python
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
Create regular voxel grid from points in an ASCIItable (or geom file).

""", version = scriptID)

parser.add_option('-m', '--mode',
                  dest = 'mode',
                  type = 'choice', choices = ['cell','point'],
                  help = 'cell-centered or point-centered coordinates ')
parser.add_option('-c', '--coordinates',
                  dest = 'coords',
                  type = 'string', metavar = 'string',
                  help = 'coordinate label [%default]')
parser.set_defaults(coords = 'pos',
                    mode   = 'cell'
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  isGeom = name.endswith('.geom')
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False,
                                    labeled  = not isGeom,
                                    readonly = True,
                                   )
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  remarks = []
  errors  = []
  coordDim = 3 if isGeom else table.label_dimension(options.coords)
  if not 3 >= coordDim >= 1: errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.coords))
  elif coordDim < 3:         remarks.append('appending {} dimensions to coordinates "{}"...'.format(3-coordDim,options.coords))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  if isGeom:
    info,extra_header = table.head_getGeom()
    coords = [np.linspace(info['origin'][i],
                          info['origin'][i]+info['size'][i],
                          num = info['grid'][i]+1,
                          endpoint = True,
                         ) for i in xrange(3)]
  else:
    table.data_readArray(options.coords)
    if len(table.data.shape) < 2: table.data.shape += (1,)                                            # expand to 2D shape
    if table.data.shape[1] < 3:
      table.data = np.hstack((table.data,
                              np.zeros((table.data.shape[0],
                                        3-table.data.shape[1]),dtype='f')))                           # fill coords up to 3D with zeros

    coords = [np.unique(table.data[:,i]) for i in xrange(3)]

    if options.mode == 'cell':
      coords = [0.5 * np.array([3.0 * coords[i][0] - coords[i][0 + len(coords[i]) > 1]] + \
                               [coords[i][j-1] + coords[i][j] for j in xrange(1,len(coords[i]))] + \
                               [3.0 * coords[i][-1] - coords[i][-1 - (len(coords[i]) > 1)]]) for i in xrange(3)]

  grid = np.array(map(len,coords),'i')
  N = grid.prod() if options.mode == 'point' or isGeom else (grid-1).prod()

  if not isGeom and N != len(table.data):
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
    (directory,filename) = os.path.split(name)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(directory,os.path.splitext(filename)[0] +
                                              ('' if isGeom else '_{}({})'.format(options.coords, options.mode)) +
                                              '.' + writer.GetDefaultFileExtension()))
  else:
    writer = vtk.vtkDataSetWriter()
    writer.WriteToOutputStringOn()
    writer.SetHeader('# powered by '+scriptID)

  if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(rGrid)
  else:                          writer.SetInputData(rGrid)
  writer.Write()
  if name is None:  sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()
