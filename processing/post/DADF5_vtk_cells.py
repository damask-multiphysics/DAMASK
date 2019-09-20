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
                    help='name of subdirectory relative to the location of the DADF5 file to hold output')
parser.add_argument('--mat', nargs='+',
                    help='labels for materialpoint',dest='mat')
parser.add_argument('--con', nargs='+',
                    help='labels for constituent',dest='con')

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


  for i,inc in enumerate(results.iter_visible('increments')):
    print('Output step {}/{}'.format(i+1,len(results.increments)))
    vtk_data = []
    
    results.set_visible('materialpoints',False)
    results.set_visible('constituents',  True)
    for label in options.con:
      
      for p in results.iter_visible('con_physics'):
        if p != 'generic':
          for c in results.iter_visible('constituents'):
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
          vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
          rGrid.GetCellData().AddArray(vtk_data[-1])
    
    results.set_visible('constituents',  False)
    results.set_visible('materialpoints',True)
    for label in options.mat:       
      for p in results.iter_visible('mat_physics'):
        if p != 'generic':
          for m in results.iter_visible('materialpoints'):
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
          vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
          rGrid.GetCellData().AddArray(vtk_data[-1])
          
    if results.structured:
      writer = vtk.vtkXMLRectilinearGridWriter()


    dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
    if not os.path.isdir(dirname):
      os.mkdir(dirname,0o755)
    file_out = '{}_{}.{}'.format(os.path.splitext(os.path.split(filename)[-1])[0],inc,writer.GetDefaultFileExtension())
    
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(dirname,file_out))
    if results.structured:
      writer.SetInputData(rGrid)
    
    writer.Write()