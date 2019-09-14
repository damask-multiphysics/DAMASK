#!/usr/bin/env python3

import os
import argparse

import numpy as np
import vtk
from vtk.util import numpy_support

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = argparse.ArgumentParser()

#ToDo:  We need to decide on a way of handling arguments of variable lentght
#https://stackoverflow.com/questions/15459997/passing-integer-lists-to-python

#parser.add_argument('--version', action='version', version='%(prog)s {}'.format(scriptID))
parser.add_argument('filenames', nargs='+',
                    help='DADF5 files')
parser.add_argument('-d','--dir', dest='dir',default='postProc',metavar='string',
                    help='name of subdirectory to hold output')
parser.add_argument('--mat', nargs='+',
                    help='labels for materialpoint/homogenization',dest='mat')
parser.add_argument('--con', nargs='+',
                    help='labels for constituent/crystallite/constitutive',dest='con')

options = parser.parse_args()

if options.mat is None: options.mat=[]
if options.con is None: options.con=[]

# --- loop over input files ------------------------------------------------------------------------

for filename in options.filenames:
  results = damask.DADF5(filename)
  
  if results.structured:                                                                               # for grid solvers use rectilinear grid
    rGrid = vtk.vtkRectilinearGrid()
    coordArray = [vtk.vtkDoubleArray(),
                  vtk.vtkDoubleArray(),
                  vtk.vtkDoubleArray(),
                 ]

    rGrid.SetDimensions(*(results.grid+1))
    for dim in [0,1,2]:
      for c in np.linspace(0,results.size[dim],1+results.grid[dim]):
        coordArray[dim].InsertNextValue(c)

    rGrid.SetXCoordinates(coordArray[0])
    rGrid.SetYCoordinates(coordArray[1])
    rGrid.SetZCoordinates(coordArray[2])


  for i,inc in enumerate(results.increments):
    print('Output step {}/{}'.format(i+1,len(results.increments)))
    vtk_data = []
    results.active['increments'] = [inc]
    
    results.materialpoint_set([])
    results.constituent_set(results.constituents)
    for label in options.con:
      
      for o in results.constituent_output_iter():
        if o != 'generic':
          for c in results.constituent_iter():
            x = results.get_dataset_location(label)
            if len(x) == 0:
              continue
            array = results.read_dataset(x,0)
            shape = [array.shape[0],np.product(array.shape[1:])]
            vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
            vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
            rGrid.GetCellData().AddArray(vtk_data[-1])
        else:
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          array = results.read_dataset(x,0)
          shape = [array.shape[0],np.product(array.shape[1:])]
          vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
          vtk_data[-1].SetName('1_'+x[0].split('/')[1]+'/generic/'+label)
          rGrid.GetCellData().AddArray(vtk_data[-1])
    
    results.constituent_set([])
    results.materialpoint_set(results.materialpoints)
    for label in options.mat:       
      for o in results.materialpoint_output_iter():
        if o != 'generic':
          for m in results.materialpoint_iter():
            x = results.get_dataset_location(label)
            if len(x) == 0:
              continue
            array = results.read_dataset(x,0)
            shape = [array.shape[0],np.product(array.shape[1:])]
            vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
            vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
            rGrid.GetCellData().AddArray(vtk_data[-1])
        else:
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          array = results.read_dataset(x,0)
          shape = [array.shape[0],np.product(array.shape[1:])]
          vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
          vtk_data[-1].SetName('1_'+x[0].split('/')[1]+'/generic/'+label)
          rGrid.GetCellData().AddArray(vtk_data[-1])
          
    if results.structured:
      writer = vtk.vtkXMLRectilinearGridWriter()


    dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
    try:
      os.mkdir(dirname)
    except FileExistsError:
      pass
    file_out = '{}_inc{:04d}.{}'.format(filename.split('.')[0],inc['inc'],writer.GetDefaultFileExtension())
    
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(dirname,file_out))
    if results.structured:
      writer.SetInputData(rGrid)
    
    writer.Write()
