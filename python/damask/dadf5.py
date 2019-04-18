# -*- coding: UTF-8 no BOM -*-
import h5py
import re
import numpy as np
import os

# ------------------------------------------------------------------
class DADF5():
  """Read and write to DADF5 files"""
  
# ------------------------------------------------------------------
  def __init__(self,
               filename,
               mode     = 'r',
              ):
    
    if mode not in ['a','r']:
      print('Invalid file access mode')
      with h5py.File(filename,mode):
        pass
      
    with h5py.File(filename,'r') as f:
      
      if f.attrs['DADF5-major'] != 0 or f.attrs['DADF5-minor'] != 1:
        raise TypeError('Unsupported DADF5 version {} '.format(f.attrs['DADF5-version']))
    
      self.structured = 'grid' in f['mapping'].attrs.keys()
    
      if self.structured:
        self.grid = f['mapping'].attrs['grid']
        self.size = f['mapping'].attrs['size']
        
      r=re.compile('inc[0-9]+')
      self.increments = [{'inc':  int(u[3:]),
                          'time': round(f[u].attrs['time/s'],12),
                          }   for u in f.keys() if r.match(u)]
      
      self.constituents    = np.unique(f['mapping/cellResults/constituent']['Name']).tolist()        # ToDo: I am not to happy with the name
      self.constituents    = [c.decode() for c in self.constituents]
      
      self.materialpoints  = np.unique(f['mapping/cellResults/materialpoint']['Name']).tolist()      # ToDo: I am not to happy with the name
      self.materialpoints  = [m.decode() for m in self.materialpoints]
      
      self.Nconstituents   = [i for i in range(np.shape(f['mapping/cellResults/constituent'])[1])]
      self.Nmaterialpoints = np.shape(f['mapping/cellResults/constituent'])[0]
      
      self.c_output_types  = []
      for c in self.constituents:
        for o in f['inc{:05}/constituent/{}'.format(self.increments[0]['inc'],c)].keys():
          self.c_output_types.append(o)
      self.c_output_types = list(set(self.c_output_types))                                          # make unique
      
    self.active= {'increments':     self.increments,
                  'constituents':   self.constituents,
                  'materialpoints': self.materialpoints,
                  'constituent':    self.Nconstituents,
                  'c_output_types': self.c_output_types}

    self.filename   = filename
    self.mode       = mode
    
  def list_data(self):
    """Shows information on all datasets in the file"""
    with h5py.File(self.filename,'r') as f:
      group_inc = 'inc{:05}'.format(self.active['increments'][0]['inc'])
      for c in self.active['constituents']:
        print('\n'+c)
        group_constituent = group_inc+'/constituent/'+c
        for t in self.active['c_output_types']:
          print('  {}'.format(t))
          group_output_types = group_constituent+'/'+t
          try:
            for x in f[group_output_types].keys():
              print('    {} ({})'.format(x,f[group_output_types+'/'+x].attrs['Description'].decode()))
          except:
            pass
    

  def get_dataset_location(self,label):
    """Returns the location of all active datasets with given label"""
    path = []
    with h5py.File(self.filename,'r') as f:
      for i in self.active['increments']:
        group_inc = 'inc{:05}'.format(i['inc'])
        for c in self.active['constituents']:
          group_constituent = group_inc+'/constituent/'+c
          for t in self.active['c_output_types']:
            try:
              f[group_constituent+'/'+t+'/'+label]
              path.append(group_constituent+'/'+t+'/'+label)
            except:
              pass
    return path
    
    
  def read_dataset(self,path,c):
    """
    Dataset for all points/cells
    
    
    If more than one path is given, the dataset is composed of the individual contributions
    """
    with h5py.File(self.filename,'r') as f:
      shape = (self.Nmaterialpoints,) + np.shape(f[path[0]])[1:]
      dataset = np.full(shape,np.nan)
      for pa in path:
        label   = pa.split('/')[2]
        p = np.where(f['mapping/cellResults/constituent'][:,c]['Name'] == str.encode(label))[0]
        u = (f['mapping/cellResults/constituent'][p,c]['Position'])
        dataset[p,:] = f[pa][u,:]                                                    
    
    return dataset
        
      
