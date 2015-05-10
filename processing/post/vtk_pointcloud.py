#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add grain index based on similitude of crystal lattice orientation.

""", version = scriptID)

parser.add_option('-p', '--positions', dest='pos', metavar='string',
                  help = 'coordinate label [%default]')
parser.set_defaults(pos = 'pos')

(options, filenames) = parser.parse_args()

datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
           }

datainfo['vector']['label'] += [options.pos]

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':os.path.splitext(name)[0]+'.vtp', 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['croak'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info

# --------------- figure out columns to process
  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          if datatype not in active: active[datatype] = []
          if datatype not in column: column[datatype] = {}
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                 # remember columns of requested data
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
        break


# ------------------------------------------ process data ---------------------------------------  

  Polydata = vtk.vtkPolyData()
  Points = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()
  table.data_readArray(range(column['vector'][options.pos],\
                             column['vector'][options.pos]+datainfo['vector']['len']))

  for p in table.data:
    id = Points.InsertNextPoint(p)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)

  Polydata.SetPoints(Points)
  Polydata.SetVerts(Vertices)
  Polydata.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5:
    Polydata.Update()
 
# ------------------------------------------ output result ---------------------------------------  

  writer = vtk.vtkXMLPolyDataWriter()
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.SetFileName(file['output'])
  if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(Polydata)
  else:
      writer.SetInputData(Polydata)
  writer.Write()

  table.input_close()                                                     # close input ASCII table
