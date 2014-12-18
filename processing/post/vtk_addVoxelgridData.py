#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import damask
from collections import defaultdict
from optparse import OptionParser

scriptID   = string.replace('$Id: addCalculation.py 3465 2014-09-12 14:14:55Z MPIE\m.diehl $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add scalar and RGB tuples from ASCIItable to existing VTK voxel grid (.vtr).

""", version = scriptID)

parser.add_option('-v', '--vtk',     dest='vtk',    type='string',
                  help = 'VTK file name')
parser.add_option('-s', '--scalar',  dest='scalar', action='extend',
                  help = 'scalar values')
parser.add_option('-c', '--color',   dest='color',  action='extend',
                  help = 'RGB color tuples')
parser.add_option('-r', '--render',  dest='render', action='store_true',
                  help = 'open output in VTK render window')

parser.set_defaults(scalar = [])
parser.set_defaults(color = [])
parser.set_defaults(render = False)

(options, filenames) = parser.parse_args()

datainfo = {                                                               # list of requested labels per datatype
             'scalar':     {'len':1,
                            'label':[]},
             'color':      {'len':3,
                            'label':[]},
           }

if options.vtk == None or not os.path.exists(options.vtk):
  parser.error('VTK file does not exist')

if os.path.splitext(options.vtk)[1] == '.vtr':
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  rGrid = reader.GetOutput()
elif os.path.splitext(options.vtk)[1] == '.vtk':
  reader = vtk.vtkGenericDataObjectReader()
  reader.SetFileName(options.vtk)
  reader.Update()
  rGrid = reader.GetRectilinearGridOutput()
else:
  parser.error('unsupported VTK file type extension')


Npoints = rGrid.GetNumberOfPoints()
Ncells  = rGrid.GetNumberOfCells()

#if Npoints != Ncells or Npoints != Nvertices:
#  parser.error('Number of points, cells, and vertices in VTK differ from each other'); sys.exit()
if options.scalar != None:  datainfo['scalar']['label'] += options.scalar
if options.color  != None:  datainfo['color']['label']  += options.color

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':sys.stderr, 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info

# --------------- figure out columns to process
  active = defaultdict(list)
  column = defaultdict(dict)
  array = defaultdict(dict)
  
  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
          if datatype == 'scalar':
            array[datatype][label] = vtk.vtkDoubleArray()
            array[datatype][label].SetNumberOfComponents(1)
            array[datatype][label].SetName(label)
          elif datatype == 'color':
            array[datatype][label] = vtk.vtkUnsignedCharArray()
            array[datatype][label].SetNumberOfComponents(3)
            array[datatype][label].SetName(label)
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table
    
    for datatype,labels in active.items():                                  # loop over scalar,color
      for label in labels:                                                  # loop over all requested items
        theData = table.data[column[datatype][label]:\
                             column[datatype][label]+datainfo[datatype]['len']]  # read strings
        if datatype == 'color':
          theData = map(lambda x: int(255.*float(x)),theData)
          array[datatype][label].InsertNextTuple3(theData[0],theData[1],theData[2],)
        elif datatype == 'scalar':
          array[datatype][label].InsertNextValue(float(theData[0]))

  file['input'].close()                                                   # close input ASCII table

# ------------------------------------------ add data ---------------------------------------  

  for datatype,labels in active.items():                                   # loop over scalar,color
    if datatype == 'color':
#      Polydata.GetPointData().SetScalars(array[datatype][labels[0]])
      rGrid.GetCellData().SetScalars(array[datatype][labels[0]])
    for label in labels:                                                   # loop over all requested items
#      Polydata.GetPointData().AddArray(array[datatype][label])
      rGrid.GetCellData().AddArray(array[datatype][label])

  rGrid.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5:
    rGrid.Update()

# ------------------------------------------ output result ---------------------------------------  

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetDataModeToBinary()
writer.SetCompressorTypeToZLib()
writer.SetFileName(os.path.splitext(options.vtk)[0]+'_added.vtr')
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(rGrid)
else:
    writer.SetInputData(rGrid)
writer.Write()

# ------------------------------------------ render result ---------------------------------------  

if options.render:
  mapper = vtk.vtkDataSetMapper()
  mapper.SetInputData(rGrid)
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)

  # Create the graphics structure. The renderer renders into the
  # render window. The render window interactor captures mouse events
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
 
  #ren.ResetCamera()
  #ren.GetActiveCamera().Zoom(1.5)
 
  iren.Initialize()
  renWin.Render()
  iren.Start()
