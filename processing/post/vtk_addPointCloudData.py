#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import vtk
from vtk.util import numpy_support

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [ASCIItable(s)]',
                      description = "Add scalars, vectors, tensors, and/or an RGB tuple from ASCIItable "
                                  + "VTK point cloud (.vtp).",
                      version = scriptID)

parser.add_option(      '--vtk',
                  dest = 'vtk',
                  type = 'string', metavar = 'string',
                  help = 'VTK file name')
parser.add_option('-r', '--render',
                  dest = 'render',
                  action = 'store_true',
                  help = 'open output in VTK render window')
parser.add_option('-d', '--data',
                  dest = 'data',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'scalar/vector value(s) label(s)')
parser.add_option('-t', '--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'tensor (3x3) value label(s)')
parser.add_option('-c', '--color',
                  dest = 'color',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'RGB color tuple label')

parser.set_defaults(data = [],
                    tensor = [],
                    color = [],
)

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if not options.vtk:                 parser.error('No VTK file specified.')
if not os.path.exists(options.vtk): parser.error('VTK file does not exist.')

vtk_file,vtk_ext = os.path.splitext(options.vtk)

if vtk_ext == '.vtp':
  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  Polydata = reader.GetOutput()
elif vtk_ext == '.vtk':
  reader = vtk.vtkGenericDataObjectReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  Polydata = reader.GetPolyDataOutput()
else:
  parser.error('unsupported VTK file type extension.')

Npoints   = Polydata.GetNumberOfPoints()
Ncells    = Polydata.GetNumberOfCells()
Nvertices = Polydata.GetNumberOfVerts()

if Npoints != Ncells or Npoints != Nvertices:
  parser.error('number of points, cells, and vertices in VTK differ from each other.')

damask.util.croak('{}: {} points/vertices/cells...'.format(options.vtk,Npoints))

for name in filenames:
  damask.util.report(scriptName,name)

  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
  
  VTKarray = {}
  for data in options.data:
    VTKarray[data] = numpy_support.numpy_to_vtk(table.get(data).copy(),
                                                deep=True,array_type=vtk.VTK_DOUBLE)
    VTKarray[data].SetName(data)
  
  for color in options.color:
    VTKarray[color] = numpy_support.numpy_to_vtk((table.get(color)*255).astype(int).copy(),
                                                deep=True,array_type=vtk.VTK_UNSIGNED_CHAR)
    VTKarray[color].SetName(color)

  for tensor in options.tensor:
    data = damask.mechanics.symmetric(table.get(tensor).reshape((-1,3,3))).reshape((-1,9))
    VTKarray[tensor] = numpy_support.numpy_to_vtk(data.copy(),
                                                deep=True,array_type=vtk.VTK_DOUBLE)
    VTKarray[tensor].SetName(tensor)


  for data in VTKarray:
     Polydata.GetPointData().AddArray(VTKarray[data])
  Polydata.Modified()

# ------------------------------------------ output result ---------------------------------------

  writer = vtk.vtkXMLPolyDataWriter()
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.SetFileName(vtk_file+'.'+writer.GetDefaultFileExtension())
  writer.SetInputData(Polydata)
  writer.Write()

# ------------------------------------------ render result ---------------------------------------

if options.render:
  mapper = vtk.vtkDataSetMapper()
  mapper.SetInputData(Polydata)
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)

# Create the graphics structure. The renderer renders into the
# render window. The render window interactively captures mouse events
# and will perform appropriate camera or actor manipulation
# depending on the nature of the events.

  ren = vtk.vtkRenderer()

  renWin = vtk.vtkRenderWindow()
  renWin.AddRenderer(ren)

  ren.AddActor(actor)
  ren.SetBackground(1, 1, 1)
  renWin.SetSize(200, 200)

  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)

  iren.Initialize()
  renWin.Render()
  iren.Start()
