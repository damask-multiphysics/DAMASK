#!/usr/bin/env python3

import os
import argparse
import re

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
  
  Points   = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()  
  for c in results.cell_coordinates():
    pointID = Points.InsertNextPoint(c)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(pointID)

  Polydata = vtk.vtkPolyData()
  Polydata.SetPoints(Points)
  Polydata.SetVerts(Vertices)
  Polydata.Modified()
  
  N_digits = int(np.floor(np.log10(int(results.increments[-1][3:]))))+1
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
            Polydata.GetCellData().AddArray(vtk_data[-1])
        else:
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          ph_name = re.compile(r'(\/[1-9])_([A-Z][a-z]*)_(([a-z]*)|([A-Z]*))')  #looking for phase name in dataset name
          array = results.read_dataset(x,0)
          shape = [array.shape[0],np.product(array.shape[1:])]
          vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
          dset_name = '1_' + re.sub(ph_name,r'',x[0].split('/',1)[1])           #removing phase name from generic dataset
          vtk_data[-1].SetName(dset_name)
          Polydata.GetCellData().AddArray(vtk_data[-1])
    
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
            Polydata.GetCellData().AddArray(vtk_data[-1])
        else:
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          array = results.read_dataset(x,0)
          shape = [array.shape[0],np.product(array.shape[1:])]
          vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),deep=True,array_type= vtk.VTK_DOUBLE))
          vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
          Polydata.GetCellData().AddArray(vtk_data[-1])
          
    writer = vtk.vtkXMLPolyDataWriter()


    dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
    if not os.path.isdir(dirname):
      os.mkdir(dirname,0o755)
    file_out = '{}_inc{}.{}'.format(os.path.splitext(os.path.split(filename)[-1])[0],
                                   inc[3:].zfill(N_digits),
                                   writer.GetDefaultFileExtension())
    
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(dirname,file_out))
    writer.SetInputData(Polydata)
    
    writer.Write()
