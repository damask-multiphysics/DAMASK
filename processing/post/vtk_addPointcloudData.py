#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,vtk
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add scalar and RGB tuples from ASCIItable to existing VTK point cloud (.vtp).

""", version = scriptID)

parser.add_option('-v', '--vtk',   dest='vtk', \
                  help = 'VTK file name')
parser.add_option('-s', '--scalar',   dest='scalar', action='extend', \
                  help = 'scalar values')
parser.add_option('-c', '--color',   dest='color', action='extend', \
                  help = 'RGB color tuples')

parser.set_defaults(scalar = [])
parser.set_defaults(color = [])

(options, filenames) = parser.parse_args()

datainfo = {                                                               # list of requested labels per datatype
             'scalar':     {'len':1,
                            'label':[]},
             'color':      {'len':3,
                            'label':[]},
           }

if not os.path.exists(options.vtk):
  parser.error('VTK file does not exist'); sys.exit()

reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(options.vtk)
reader.Update()
Npoints = reader.GetNumberOfPoints()
Ncells = reader.GetNumberOfCells()
Nvertices = reader.GetNumberOfVerts()
Polydata = reader.GetOutput()

if Npoints != Ncells or Npoints != Nvertices:
  parser.error('Number of points, cells, and vertices in VTK differ from each other'); sys.exit()
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
  active = {}
  column = {}

  array = {}
  
  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          if datatype not in active: active[datatype] = []
          if datatype not in column: column[datatype] = {}
          if datatype not in array:   array[datatype] = {}
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

  table.input_close()                                                     # close input ASCII table

# ------------------------------------------ add data ---------------------------------------  

  for datatype,labels in active.items():                                   # loop over scalar,color
    if datatype == 'color':
      Polydata.GetPointData().SetScalars(array[datatype][labels[0]])
      Polydata.GetCellData().SetScalars(array[datatype][labels[0]])
    for label in labels:                                                   # loop over all requested items
      Polydata.GetPointData().AddArray(array[datatype][label])
      Polydata.GetCellData().AddArray(array[datatype][label])

  Polydata.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5:
    Polydata.Update()

# ------------------------------------------ output result ---------------------------------------  

writer = vtk.vtkXMLPolyDataWriter()
writer.SetDataModeToBinary()
writer.SetCompressorTypeToZLib()
writer.SetFileName(os.path.splitext(options.vtk)[0]+'_added.vtp')
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(Polydata)
else:
    writer.SetInputData(Polydata)
writer.Write()
