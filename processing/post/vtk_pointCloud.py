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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Produce a VTK point cloud dataset based on coordinates given in an ASCIItable.

""", version = scriptID)

parser.add_option('-p',
                  '--pos', '--position',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')

parser.set_defaults(pos = 'pos',
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False,
                                    readonly = True)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  errors =  []
  remarks = []
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

# ------------------------------------------ process data ---------------------------------------  

  table.data_readArray(options.pos)
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
 
# ------------------------------------------ output result ---------------------------------------  

  if name:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(os.path.split(name)[0],
                                    os.path.splitext(os.path.split(name)[1])[0] +
                                    '.' + writer.GetDefaultFileExtension()))
  else:
    writer = vtk.vtkDataSetWriter()
    writer.SetHeader('# powered by '+scriptID)
    writer.WriteToOutputStringOn()
  
  
  writer.SetInputData(Polydata)

  writer.Write()

  if name is None: sys.stdout.write(writer.GetOutputString())

  table.close()
