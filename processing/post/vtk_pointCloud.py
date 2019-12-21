#!/usr/bin/env python3

import os
import sys
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
if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

# ------------------------------------------ process data ---------------------------------------  

  Polydata = vtk.vtkPolyData()
  Points   = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()

  for p in table.get(options.pos):
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
