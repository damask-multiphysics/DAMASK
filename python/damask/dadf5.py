# -*- coding: UTF-8 no BOM -*-
import h5py
import re

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
        print('Unsupported DADF5 version {} '.format(f.attrs['DADF5-version']))
    
      self.structured = 'grid' in f['mapping'].attrs.keys()
    
      if self.structured:
        self.grid = f['mapping'].attrs['grid']
        self.size = f['mapping'].attrs['size']
        
      r=re.compile('inc[0-9]+')
      self.increments = [{'group':  u,
                          'time':   f[u].attrs['time/s'],
                          'active': True
                          }   for u in f.keys() if r.match(u)]

    self.filename   = filename
    self.mode       = mode
