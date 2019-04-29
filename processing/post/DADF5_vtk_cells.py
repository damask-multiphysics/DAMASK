#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,vtk
import numpy as np
import argparse
import damask
from vtk.util import numpy_support

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

options = parser.parse_args()

options.labels = ['Fe','Fp','xi_sl']

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
    for label in options.labels:
      for o in results.c_output_types:
        results.active['c_output_types'] = [o]
        if o != 'generic':
          for c in results.constituents:
            results.active['constituents'] = [c]
            x = results.get_dataset_location(label)
            if len(x) == 0:
              continue
            array = results.read_dataset(x,0)
            shape = [array.shape[0],np.product(array.shape[1:])]
            vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
            vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
            rGrid.GetCellData().AddArray(vtk_data[-1])
        else:
          results.active['constituents'] = results.constituents
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

    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(os.path.split(filename)[0],
                                    os.path.splitext(os.path.split(filename)[1])[0] +
                                    '_inc{:04d}'.format(i) +                                        # ToDo: adjust to length of increments
                                    '.' + writer.GetDefaultFileExtension()))
    if results.structured:
      writer.SetInputData(rGrid)
    
    writer.Write()
