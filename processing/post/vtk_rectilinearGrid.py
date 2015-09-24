#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create regular voxel grid from points in an ASCIItable.

""", version = scriptID)

parser.add_option('-m', '--mode',
                  dest = 'mode',
                  type = 'choice', choices = ['cell','point'],
                  help = 'cell-centered or point-centered coordinates ')
parser.add_option('-p', '--positions',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'coordinate label [%default]')
parser.set_defaults(position = 'pos',
                   )

(options, filenames) = parser.parse_args()

if not options.mode: parser.error("No coordinate type specified.")
# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, readonly = True)
  except: continue
  damask.util.croak(damask.util.emph(scriptName)+(': '+name if name else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  
# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  
  if table.label_dimension(options.position) != 3:  errors.append('coordinates {} are not a vector.'.format(options.position))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray(options.position)

  coords = [np.unique(table.data[:,i]) for i in xrange(3)]
  if options.mode == 'cell':
    coords = [0.5 * np.array([3.0 * coords[i][0] - coords[i][0 + len(coords[i]) > 1]] + \
                             [coords[i][j-1] + coords[i][j] for j in xrange(1,len(coords[i]))] + \
                             [3.0 * coords[i][-1] - coords[i][-1 - (len(coords[i]) > 1)]]) for i in xrange(3)]
  grid   = np.array(map(len,coords),'i')
  N = grid.prod() if options.mode == 'point' else (grid-1).prod()
  
  if N != len(table.data): errors.append('data count {} does not match grid {}x{}x{}.'.format(N,*(grid - options.mode == 'cell') ))
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

  damask.util.croak('{} points and {} cells...'.format(rGrid.GetNumberOfPoints(),rGrid.GetNumberOfCells(),))
 
# ------------------------------------------ output result ---------------------------------------  

  if name:
    (dir,filename) = os.path.split(name)
    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(dir,os.path.splitext(filename)[0]+'_{}'.format(options.position) \
                                        +'.'+writer.GetDefaultFileExtension()))
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(rGrid)
    else:                          writer.SetInputData(rGrid)
    writer.Write()

  else:
    writer = vtk.vtkRectilinearGridWriter()
    writer.WriteToOutputStringOn()
    writer.SetFileTypeToASCII()
    writer.SetHeader('# powered by '+scriptID)
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(rGrid)
    else:                          writer.SetInputData(rGrid)
    writer.Write()
    sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()                                                                                     # close input ASCII table

