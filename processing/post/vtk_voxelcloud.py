#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create hexahedral voxels around points in an ASCIItable.

""", version = scriptID)

parser.add_option('-d', '--deformed',
                  dest = 'deformed',
                  type = 'string', metavar = 'string',
                  help = 'deformed coordinate label [%default]')
parser.add_option('-c','--coordinates',
                  dest = 'coords',
                  type = 'string', metavar='string',
                  help = 'undeformed coordinates label [%default]')
parser.set_defaults(deformed = 'ipdeformedcoord',
                    coords   = 'ipinitialcoord',
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

# --------------- interprete header -----------------------------------------------------------------
  table.head_read()
  errors=[]
  if table.label_dimension(options.deformed) != 3: errors.append('columns "{}" have dimension {}'.format(options.deformed,
                                                                                                   table.label_dimension(options.deformed)))
  if table.label_dimension(options.coords)   != 3: errors.append('coordinates {} are not a vector.'.format(options.coords))

  table.data_readArray([options.coords,options.deformed])

# --------------- figure out size and grid ---------------------------------------------------------
  coords = [{},{},{}]
  for i in xrange(len(table.data)):  
    for j in xrange(3):
      coords[j][str(table.data[i,table.label_index(options.coords)+j])] = True
  grid = np.array(map(len,coords),'i')
  size = grid/np.maximum(np.ones(3,'d'),grid-1.0)* \
            np.array([max(map(float,coords[0].keys()))-min(map(float,coords[0].keys())),\
                      max(map(float,coords[1].keys()))-min(map(float,coords[1].keys())),\
                      max(map(float,coords[2].keys()))-min(map(float,coords[2].keys())),\
                      ],'d')                                                                        # size from bounding box, corrected for cell-centeredness

  size = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                               # spacing for grid==1 equal to smallest among other spacings

# ------------------------------------------ process data ---------------------------------------  
  hexPoints = np.array([[-1,-1,-1],
                        [ 1,-1,-1],
                        [ 1, 1,-1],
                        [-1, 1,-1],
                        [-1,-1, 1],
                        [ 1,-1, 1],
                        [ 1, 1, 1],
                        [-1, 1, 1],
                       ])

  halfDelta = 0.5*size/grid

  Points = vtk.vtkPoints()
  Hex    = vtk.vtkHexahedron()
  uGrid  = vtk.vtkUnstructuredGrid()
  
  for p in table.data:
    for i,h in enumerate(hexPoints):
      pointID = Points.InsertNextPoint(p[table.label_index(options.deformed):table.label_index(options.deformed)+3]+h*halfDelta)
      Hex.GetPointIds().SetId(i,pointID)

  uGrid.InsertNextCell(Hex.GetCellType(), Hex.GetPointIds())

  uGrid.SetPoints(Points)

# ------------------------------------------ output result ---------------------------------------  

  if name:
    writer = vtk.vtkXMLUnstructuredGridWriter()
    (directory,filename) = os.path.split(name)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(directory,os.path.splitext(filename)[0]
                                              +'.'+writer.GetDefaultFileExtension()))
  else:
    writer = vtk.vtkDataSetWriter()
    writer.WriteToOutputStringOn()
    writer.SetHeader('# powered by '+scriptID)
  
  if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(uGrid)
  else:                          writer.SetInputData(uGrid)
  writer.Write()
  if name == None:  sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()                                                                             # close input ASCII table

