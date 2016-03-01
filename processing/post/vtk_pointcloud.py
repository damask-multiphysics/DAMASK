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
Produce a VTK point cloud dataset based on coordinates given in an ASCIItable.

""", version = scriptID)

parser.add_option('-d', '--deformed',
                  dest = 'deformed',
                  type = 'string', metavar = 'string',
                  help = 'deformed coordinate label [%default]')

parser.set_defaults(deformed = 'ipdeformedcoord'
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False,
                              readonly = True)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  errors =  []
  remarks = []
  coordDim = table.label_dimension(options.deformed)
  if not 3 >= coordDim >= 1: errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.deformed))
  elif coordDim < 3:         remarks.append('appending {} dimensions to coordinates "{}"...'.format(3-coordDim,options.deformed))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# ------------------------------------------ process data ---------------------------------------  

  table.data_readArray(options.deformed)
  if len(table.data.shape) < 2: table.data.shape += (1,)                                            # expand to 2D shape
  if table.data.shape[1] < 3:
    table.data = np.hstack((table.data,
                            np.zeros((table.data.shape[0],
                                      3-table.data.shape[1]),dtype='f')))                           # fill coords up to 3D with zeros

  Polydata = vtk.vtkPolyData()
  Points   = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()

  for p in table.data:
    pointID = Points.InsertNextPoint(p)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(pointID)

  Polydata.SetPoints(Points)
  Polydata.SetVerts(Vertices)
  Polydata.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5: Polydata.Update()
 
# ------------------------------------------ output result ---------------------------------------  

  if name:
    writer = vtk.vtkXMLPolyDataWriter()
    (directory,filename) = os.path.split(name)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(directory,os.path.splitext(filename)[0]
                                              +'.'+writer.GetDefaultFileExtension()))
  else:
    writer = vtk.vtkDataSetWriter()
    writer.WriteToOutputStringOn()
    writer.SetHeader('# powered by '+scriptID)
  
  if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(Polydata)
  else:                          writer.SetInputData(Polydata)
  writer.Write()
  if name == None:  sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()
