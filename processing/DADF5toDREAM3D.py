#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import h5py
import numpy as np

import damask
from damask import Rotation
from damask import Orientation

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

class DAMASKtoDREAM3D():
    """
    This class can convert the DAMASK data to DREAM3D compatible data. 
    There can be various different types of ways DAMASK data can be represented. 
    Therefore, there are multiple functions available for different purposes.
    """
    def __init__(self,job_file,geom_file,load_file):
        """
        Defining the common quantities for all the functions in this class.
       
        Parameters
        ----------
        job_file: str or pathlib.Path
            Full path of the DAMASK results file.
        geom_file : str
          name of the geom file.
        load_file : 
          name of the load file.
        """
        self.job_file = Path(job_file).expanduser().absolute()
        self.geom_file = geom_file   
        self.load_file = load_file   

    def DAMASKtoDREAM3D(self,dx,inc):
        """
        Creates a dream3D file from DAMASK output.
        Without any regridding. 
        Considers the original grid from DAMASK. 
        
        Parameters:
        -----------
        dx : float
          The grid spacing.
        inc: int
          increment of interest for DREAM3D processing.
        """
        os.chdir(self.job_file.parents[0])
        #--------------------------------------------------------------------------
        #Build array of euler angles for each cell
        #--------------------------------------------------------------------------
        d = damask.Result(self.job_file)
        inc_data = d.view(increments=inc)  # selecting only relevant data to reduce overload
    
        f = h5py.File(self.job_file,'r')
        cells = f['geometry'].attrs['cells']
        size = f['geometry'].attrs['size']
        dx = size/cells
    
        O_dict = inc_data.get('O') 
        
        cell_orientation_array = np.zeros((np.prod(cells),3))
    
        phase_ID_array = np.zeros((np.prod(cells)),dtype=np.int32) #need to reshape it later
    
        for count,p in enumerate(d.phases):
            phase_index = np.where(f['cell_to/phase']['label'] == f'{p}'.encode())[0]
            if len(d.phases) > 1:
                cell_orientation_array[phase_index,:] = Rotation(O_dict[p]).as_Euler_angles()
            else:
                cell_orientation_array[phase_index,:] = Rotation(O_dict).as_Euler_angles()
            phase_ID_array[phase_index] = count + 1
        
        #--------------------------------------------------------------------------
        job_file_no_ext = os.path.splitext(self.job_file)[0]
        o = h5py.File(f'{job_file_no_ext}_increment{inc}.dream3D','w')
        o.attrs['DADF5toDREAM3D'] = '1.0'
        o.attrs['FileVersion']    = '7.0'
    
        for g in ['DataContainerBundles','Pipeline']: # empty groups (needed)
          o.create_group(g)
    
        data_container_label = 'DataContainers/SyntheticVolumeDataContainer'        
        cell_data_label      = data_container_label + '/CellData'
    
        # Data phases
        o[cell_data_label + '/Phases'] = np.reshape(phase_ID_array, \
                                                    tuple(np.flip(cells))+(1,))
    
        # Data eulers
        orientation_data = cell_orientation_array.astype(np.float32)
        o[cell_data_label + '/Eulers'] = orientation_data.reshape(tuple(np.flip(cells))+(3,))
    
        # Attributes to CellData group
        o[cell_data_label].attrs['AttributeMatrixType'] = np.array([3],np.uint32)
        o[cell_data_label].attrs['TupleDimensions']     = np.array(cells,np.uint64)
    
        # Common Attributes for groups in CellData
        for group in ['/Phases','/Eulers']:
          o[cell_data_label + group].attrs['DataArrayVersion']      = np.array([2],np.int32)
          o[cell_data_label + group].attrs['Tuple Axis Dimensions'] = 'x={},y={},z={}'.format(*np.array(cells))
        
        # phase attributes
        o[cell_data_label + '/Phases'].attrs['ComponentDimensions'] = np.array([1],np.uint64)
        o[cell_data_label + '/Phases'].attrs['ObjectType']          = 'DataArray<int32_t>'
        o[cell_data_label + '/Phases'].attrs['TupleDimensions']     = np.array(cells,np.uint64)
        
        # Eulers attributes
        o[cell_data_label + '/Eulers'].attrs['ComponentDimensions'] = np.array([3],np.uint64)
        o[cell_data_label + '/Eulers'].attrs['ObjectType']          = 'DataArray<float>'        
        o[cell_data_label + '/Eulers'].attrs['TupleDimensions']     = np.array(cells,np.uint64)
    
        # Create EnsembleAttributeMatrix
        ensemble_label = data_container_label + '/CellEnsembleData'
    
        # Data CrystalStructures
        #o[ensemble_label + '/CrystalStructures'] = np.uint32(np.array([999,1]))
        o[ensemble_label + '/CrystalStructures'] = np.uint32(np.array([999] + [1]*len(d.phases)))
        # assuming only cubic crystal structures
        # Damask can give the crystal structure info but need to look into dream3d which crystal structure corresponds to which number
        o[ensemble_label + '/PhaseTypes']        = np.uint32(np.array([999] + [Phase_types['Primary']]*len(d.phases))).reshape((len(d.phases)+1,1))
        # also assuming Primary phases
        # there can be precipitates etc as well
    
        # Attributes Ensemble Matrix
        o[ensemble_label].attrs['AttributeMatrixType'] = np.array([11],np.uint32)
        o[ensemble_label].attrs['TupleDimensions']     = np.array([len(d.phases) + 1], np.uint64)
    
        # Attributes for data in Ensemble matrix
        for group in ['CrystalStructures','PhaseTypes']: # 'PhaseName' not required MD: But would be nice to take the phase name mapping
          o[ensemble_label+'/'+group].attrs['ComponentDimensions']   = np.array([1],np.uint64)
          o[ensemble_label+'/'+group].attrs['Tuple Axis Dimensions'] = f'x={len(d.phases)+1}'
          o[ensemble_label+'/'+group].attrs['DataArrayVersion']      = np.array([2],np.int32)
          o[ensemble_label+'/'+group].attrs['ObjectType']            = 'DataArray<uint32_t>'
          o[ensemble_label+'/'+group].attrs['TupleDimensions']       = np.array([len(d.phases) + 1],np.uint64)
    
        # Create geometry info
        geom_label = data_container_label + '/_SIMPL_GEOMETRY'
        
        o[geom_label + '/DIMENSIONS'] = np.int64(np.array(cells))
        o[geom_label + '/ORIGIN']     = np.float32(np.zeros(3))
        o[geom_label + '/SPACING']    = np.float32(dx)
            
        o[geom_label].attrs['GeometryName']     = 'ImageGeometry'
        o[geom_label].attrs['GeometryTypeName'] = 'ImageGeometry'
        o[geom_label].attrs['GeometryType']          = np.array([0],np.uint32) 
        o[geom_label].attrs['SpatialDimensionality'] = np.array([3],np.uint32) 
        o[geom_label].attrs['UnitDimensionality']    = np.array([3],np.uint32) 
