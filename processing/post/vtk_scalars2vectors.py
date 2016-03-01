#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,shutil
import damask
from optparse import OptionParser
from vtk import *

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """

""", version = scriptID)

parser.add_option('-v','--vector', nargs=3, dest='vector', \
                  help='suffices indicating vector components [%default]')
parser.add_option('-s','--separator', dest='separator', \
                  help='separator between label and suffix [%default]')

parser.set_defaults(vector = ['x','y','z'])
parser.set_defaults(separator = '.')

(options, filenames) = parser.parse_args()


# --- sanity checks

if filenames == []:
  parser.print_help()
  parser.error('no file specified...')

for filename in filenames:
  if not os.path.isfile(filename):
    parser.print_help()
    parser.error('invalid file "%s" specified...'%filename)


# --- ITERATE OVER FILES AND PROCESS THEM

for filename in filenames:

  # Read the source file
  
  sys.stdout.write('read file "%s" ...'%filename)
  sys.stdout.flush()
  suffix = os.path.splitext(filename)[1]
  if suffix == '.vtk':
    reader = vtkUnstructuredGridReader()
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
  elif suffix == '.vtu':
    reader = vtkXMLUnstructuredGridReader()
  else:
    parser.error('filetype "%s" not supported'%suffix)
  reader.SetFileName(filename)
  reader.Update() 
  uGrid = reader.GetOutput()
  sys.stdout.write(' done\n')
  sys.stdout.flush()

 
  # Read the scalar data 

  scalarData = {}
  scalarsToBeRemoved = []
  Nscalars = uGrid.GetCellData().GetNumberOfArrays()
  for i in range(Nscalars):
    sys.stdout.write("\rread scalar data %d%%" %(100*i/Nscalars))
    sys.stdout.flush()
    scalarName = uGrid.GetCellData().GetArrayName(i)
    if scalarName.split(options.separator)[-1] in options.vector:
      label,suffix = scalarName.split(options.separator)
      if label not in scalarData:
        scalarData[label] = [[],[],[]]
      uGrid.GetCellData().SetActiveScalars(scalarName)
      scalarData[label][options.vector.index(suffix)] = uGrid.GetCellData().GetScalars(scalarName)
      scalarsToBeRemoved.append(scalarName)
  for scalarName in scalarsToBeRemoved:
    uGrid.GetCellData().RemoveArray(scalarName)
  # uGrid.UpdateData()
  sys.stdout.write('\rread scalar data done\n')
  sys.stdout.flush()


  # Convert the scalar data to vector data

  NscalarData = len(scalarData)
  for n,label in enumerate(scalarData):
    sys.stdout.write("\rconvert to vector data %d%%" %(100*n/NscalarData))
    sys.stdout.flush()
    Nvalues = scalarData[label][0].GetNumberOfTuples()
    vectorData = vtkDoubleArray()
    vectorData.SetName(label)
    vectorData.SetNumberOfComponents(3) # set this before NumberOfTuples !!!
    vectorData.SetNumberOfTuples(Nvalues)
    for i in range(Nvalues):
      for j in range(3):
        vectorData.SetComponent(i,j,scalarData[label][j].GetValue(i))
    uGrid.GetCellData().AddArray(vectorData)
    # uGrid.GetCellData().SetActiveVectors(label)
  sys.stdout.write('\rconvert to vector data done\n')
  
  
  # Write to new vtk file

  outfilename = os.path.splitext(filename)[0]+'.vtu'
  sys.stdout.write('write to file "%s" ...'%outfilename)
  sys.stdout.flush()
  writer = vtkXMLUnstructuredGridWriter()
  writer.SetFileName(outfilename+'_tmp')
  writer.SetDataModeToAscii()
  writer.SetInput(uGrid)
  writer.Write()
  sys.stdout.write(' done\n')
  sys.stdout.flush()
  shutil.move(outfilename+'_tmp',outfilename)
  


# ---------------------------       DONE     --------------------------------
