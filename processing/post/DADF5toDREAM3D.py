#!/usr/bin/env python3

import argparse
import os

import h5py
import numpy as np

import damask

class AttributeManagerNullterm(h5py.AttributeManager): 
  """
  Attribute management for DREAM.3D hdf5 files.
  
  String attribute values are stored as fixed-length string with NULLTERM
  
  References
  ----------
    https://stackoverflow.com/questions/38267076
    https://stackoverflow.com/questions/52750232

  """ 

  def create(self, name, data, shape=None, dtype=None):
    if isinstance(data,str):
      tid = h5py.h5t.C_S1.copy()
      tid.set_size(len(data + ' '))
      super().create(name=name,data=data+' ',dtype = h5py.Datatype(tid))
    else:
      super().create(name=name,data=data,shape=shape,dtype=dtype)
     

h5py._hl.attrs.AttributeManager = AttributeManagerNullterm # 'Monkey patch'


# --------------------------------------------------------------------
# Crystal structure specifications
# --------------------------------------------------------------------
Crystal_structures = {'fcc': 1,
                      'bcc': 1,
                      'hcp': 0,
                      'bct': 7,
                      'ort': 6} #TODO: is bct Tetragonal low/Tetragonal high?
Phase_types = {'Primary': 0} #further additions to these can be done by looking at 'Create Ensemble Info' filter


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Creating a file for DREAM3D from DAMASK data')
parser.add_argument('filenames',nargs='+',help='HDF5 based output file')
parser.add_argument('--inc',nargs='+',help='Increment for which DREAM3D to be used, eg. 00025',type=int)
parser.add_argument('-d','--dir', dest='dir',default='postProc',metavar='string',
                    help='name of subdirectory to hold output')

options = parser.parse_args()

# --------------------------------------------------------------------
# loop over input files
for filename in options.filenames:
  f = damask.DADF5(filename)  #DAMASK output file
  count = 0
  for increment in f.increments:
    if int(increment[3:]) not in options.inc: 
      count = count + 1
      continue

    #-------output file creation-------------------------------------
    dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
    print(dirname)
    try:
      os.mkdir(dirname)
    except FileExistsError:
      pass

    o = h5py.File(dirname + '/' + os.path.splitext(filename)[0] + '_{}.dream3D'.format(increment),'w')
    #-----------------------------------------------------------------
    o.attrs['DADF5toDREAM3D'] = '1.0'
    o.attrs['FileVersion']    = '7.0' 
    #-----------------------------------------------------------------
    

    for g in ['DataContainerBundles','Pipeline']: # empty groups (needed)
      o.create_group(g)

    data_container_label = 'DataContainers/ImageDataContainer'        
    cell_data_label      = data_container_label + '/CellData'


    # Phase information of DREAM.3D is constituent ID in DAMASK
    o[cell_data_label + '/Phases'] = f.get_constituent_ID().reshape(tuple(f.grid)+(1,))  
    # Data quaternions
    DAMASK_quaternion = f.read_dataset(f.get_dataset_location('orientation'),0)
    DREAM_3D_quaternion = np.empty((np.prod(f.grid),4),dtype=np.float32)
    # Convert: DAMASK uses P = -1, DREAM.3D uses P = +1. Also change position of imagninary part
    DREAM_3D_quaternion = np.hstack((-DAMASK_quaternion['x'],-DAMASK_quaternion['y'],-DAMASK_quaternion['z'],
                                      DAMASK_quaternion['w']))
    o[cell_data_label + '/Quats'] = DREAM_3D_quaternion.reshape(tuple(f.grid)+(4,))
    
    # Attributes to CellData group
    o[cell_data_label].attrs['AttributeMatrixType'] = np.array([3],np.uint32)
    o[cell_data_label].attrs['TupleDimensions']     = f.grid.astype(np.uint64)
  
    # Common Attributes for groups in CellData
    for group in ['/Phases','/Quats']:
      o[cell_data_label + group].attrs['DataArrayVersion']      = np.array([2],np.int32)
      o[cell_data_label + group].attrs['Tuple Axis Dimensions'] = 'x={},y={},z={}'.format(*f.grid)
      
    # phase attributes
    o[cell_data_label + '/Phases'].attrs['ComponentDimensions'] = np.array([1],np.uint64)
    o[cell_data_label + '/Phases'].attrs['ObjectType']          = 'DataArray<int32_t>'
    
    # Quats attributes
    o[cell_data_label + '/Quats'].attrs['ComponentDimensions'] = np.array([4],np.uint64)
    o[cell_data_label + '/Quats'].attrs['ObjectType']          = 'DataArray<float>'        

    # Create EnsembleAttributeMatrix
    ensemble_label = data_container_label + '/EnsembleAttributeMatrix' 
    
    # Data CrystalStructures
    o[ensemble_label + '/CrystalStructures'] = np.uint32(np.array([999,\
                                                   Crystal_structures[f.get_crystal_structure()]])).reshape((2,1))
    o[ensemble_label + '/PhaseTypes']        = np.uint32(np.array([999,Phase_types['Primary']])).reshape((2,1))    # ToDo
   
    # Attributes Ensemble Matrix
    o[ensemble_label].attrs['AttributeMatrixType'] = np.array([11],np.uint32)
    o[ensemble_label].attrs['TupleDimensions']     = np.array([2], np.uint64)
    
    # Attributes for data in Ensemble matrix
    for group in ['CrystalStructures','PhaseTypes']: # 'PhaseName' not required MD: But would be nice to take the phase name mapping
      o[ensemble_label+'/'+group].attrs['ComponentDimensions']   = np.array([1],np.uint64)
      o[ensemble_label+'/'+group].attrs['Tuple Axis Dimensions'] = 'x=2'
      o[ensemble_label+'/'+group].attrs['DataArrayVersion']      = np.array([2],np.int32)
      o[ensemble_label+'/'+group].attrs['ObjectType']            = 'DataArray<uint32_t>'
      o[ensemble_label+'/'+group].attrs['TupleDimensions']       = np.array([2],np.uint64)
      
  
    # Create geometry info
    geom_label = data_container_label + '/_SIMPL_GEOMETRY'
    
    o[geom_label + '/DIMENSIONS'] = np.int64(f.grid)
    o[geom_label + '/ORIGIN']     = np.float32(np.zeros(3))
    o[geom_label + '/SPACING']    = np.float32(f.size)
  
    o[geom_label].attrs['GeometryName']     = 'ImageGeometry'
    o[geom_label].attrs['GeometryTypeName'] = 'ImageGeometry'
    o[geom_label].attrs['GeometryType']          = np.array([0],np.uint32) 
    o[geom_label].attrs['SpatialDimensionality'] = np.array([3],np.uint32) 
    o[geom_label].attrs['UnitDimensionality']    = np.array([3],np.uint32) 
