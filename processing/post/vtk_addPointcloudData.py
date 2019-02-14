#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,vtk
import damask
from collections import defaultdict
from optparse import OptionParser
from vtk.util import numpy_support

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description = """Add scalar and RGB tuples from ASCIItable to existing VTK point cloud (.vtp).""",
                      version = scriptID)

parser.add_option(      '--vtk',
                  dest = 'vtk',
                  type = 'string', metavar = 'string',
                  help = 'VTK file name')
parser.add_option(      '--inplace',
                  dest = 'inplace',
                  action = 'store_true',
                  help = 'modify VTK file in-place')
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
parser.add_option('-c', '--color',   dest='color', action='extend',
                  metavar ='<string LIST>',
                  help = 'RGB color tuples')

parser.set_defaults(data = [],
                    tensor = [],
                    color = [],
                    inplace = False,
                    render = False,
)

(options, filenames) = parser.parse_args()

if not options.vtk:                 parser.error('no VTK file specified.')
if not os.path.exists(options.vtk): parser.error('VTK file does not exist.')

if os.path.splitext(options.vtk)[1] == '.vtp':
  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  Polydata = reader.GetOutput()
elif os.path.splitext(options.vtk)[1] == '.vtk':
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

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False,
                                    readonly = True)
  except: continue
  damask.util.report(scriptName, name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  remarks = []
  errors  = []
  VTKarray = {}
  active = defaultdict(list)

  for datatype,dimension,label in [['data',0,options.data],
                                   ['tensor',9,options.tensor],
                                   ['color' ,3,options.color],
                                   ]:
    for i,dim in enumerate(table.label_dimension(label)):
      me = label[i]
      if dim == -1:          remarks.append('{} "{}" not found...'.format(datatype,me))
      elif dimension > 0 \
       and dim != dimension: remarks.append('"{}" not of dimension {}...'.format(me,dimension))
      else:
        remarks.append('adding {}{} "{}"...'.format(datatype if dim > 1 else 'scalar',
                                                    '' if dimension > 0 or dim == 1 else '[{}]'.format(dim),
                                                    me))
        active[datatype].append(me)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --------------------------------------- process and add data -----------------------------------

  table.data_readArray([item for sublist in active.values() for item in sublist])                 # read all requested data

  for datatype,labels in active.items():                                                          # loop over scalar,color
    for me in labels:                                                                             # loop over all requested items
      VTKtype = vtk.VTK_DOUBLE
      VTKdata = table.data[:, table.label_indexrange(me)].copy()                                  # copy to force contiguous layout

      if datatype == 'color':
        VTKtype = vtk.VTK_UNSIGNED_CHAR
        VTKdata = (VTKdata*255).astype(int)                                                       # translate to 0..255 UCHAR
      elif datatype == 'tensor':
        VTKdata[:,1] = VTKdata[:,3] = 0.5*(VTKdata[:,1]+VTKdata[:,3])
        VTKdata[:,2] = VTKdata[:,6] = 0.5*(VTKdata[:,2]+VTKdata[:,6])
        VTKdata[:,5] = VTKdata[:,7] = 0.5*(VTKdata[:,5]+VTKdata[:,7])

      VTKarray[me] = numpy_support.numpy_to_vtk(num_array=VTKdata,deep=True,array_type=VTKtype)
      VTKarray[me].SetName(me)

      if datatype == 'color':
        Polydata.GetPointData().SetScalars(VTKarray[me])
        Polydata.GetCellData().SetScalars(VTKarray[me])
      else:
        Polydata.GetPointData().AddArray(VTKarray[me])
        Polydata.GetCellData().AddArray(VTKarray[me])


  table.input_close()                                                                            # close input ASCII table

# ------------------------------------------ output result ---------------------------------------

  Polydata.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5: Polydata.Update()

  writer = vtk.vtkXMLPolyDataWriter()
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.SetFileName(os.path.splitext(options.vtk)[0]+('.vtp' if options.inplace else '_added.vtp'))
  if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(Polydata)
  else:                          writer.SetInputData(Polydata)
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
