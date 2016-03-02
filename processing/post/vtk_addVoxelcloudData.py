#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,vtk
import damask
from collections import defaultdict
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add scalar and RGB tuples from ASCIItable to existing VTK voxel cloud (.vtu/.vtk).

""", version = scriptID)

parser.add_option('-v', '--vtk',   dest='vtk', \
                  help = 'VTK file name')
parser.add_option('-s', '--scalar',   dest='scalar', action='extend', \
                  help = 'scalar values')
parser.add_option('-c', '--color',   dest='color', action='extend', \
                  help = 'RGB color tuples')

parser.set_defaults(scalar = [],
                    color = [],
                    render = False,
)

(options, filenames) = parser.parse_args()

if options.vtk is None or not os.path.exists(options.vtk):
  parser.error('VTK file does not exist')

if os.path.splitext(options.vtk)[1] == '.vtu':
  reader = vtk.vtkXMLUnstructuredGridReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  uGrid = reader.GetOutput()
elif os.path.splitext(options.vtk)[1] == '.vtk':
  reader = vtk.vtkGenericDataObjectReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  uGrid = reader.GetUnstructuredGridOutput()
else:
  parser.error('unsupported VTK file type extension')

Npoints = uGrid.GetNumberOfPoints()
Ncells  = uGrid.GetNumberOfCells()

sys.stderr.write('{}: {} points and {} cells...\n'.format(damask.util.emph(options.vtk),Npoints,Ncells))

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

  remarks = []
  errors  = []
  VTKarray = {}
  active = defaultdict(list)
  
  for datatype,dimension,label in [['scalar',1,options.scalar],
                                   ['color',3,options.color],
                                   ]:
    for i,dim in enumerate(table.label_dimension(label)):
      me = label[i]
      if dim == -1: remarks.append('{} "{}" not found...'.format(datatype,me))
      elif dim > dimension: remarks.append('"{}" not of dimension{}...'.format(me,dimension))
      else:
        damask.util.croak('adding {} {}'.format(datatype,me))
        active[datatype].append(me)

        if datatype in ['scalar']:
          VTKarray[me] = vtk.vtkDoubleArray()
        elif datatype == 'color':
          VTKarray[me] = vtk.vtkUnsignedCharArray()

        VTKarray[me].SetNumberOfComponents(dimension)
        VTKarray[me].SetName(label[i])

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue
       
# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                                          # read next data line of ASCII table
    
    for datatype,labels in active.items():                                                          # loop over scalar,color
      for me in labels:                                                                             # loop over all requested items
        theData = [table.data[i] for i in table.label_indexrange(me)]                               # read strings
        if datatype == 'color':
          VTKarray[me].InsertNextTuple3(*map(lambda x: int(255.*float(x)),theData))
        elif datatype == 'scalar':
          VTKarray[me].InsertNextValue(float(theData[0]))

# ------------------------------------------ add data ---------------------------------------  

  for datatype,labels in active.items():                                   # loop over scalar,color
    if datatype == 'color':
      uGrid.GetCellData().SetScalars(VTKarray[active['color'][0]])
    for label in labels:                                                   # loop over all requested items
      uGrid.GetCellData().AddArray(VTKarray[me])

  uGrid.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5:
    uGrid.Update()

# ------------------------------------------ output result ---------------------------------------  

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetDataModeToBinary()
writer.SetCompressorTypeToZLib()
writer.SetFileName(os.path.splitext(options.vtk)[0]+'_added.vtu')
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(uGrid)
else:
    writer.SetInputData(uGrid)
writer.Write()
