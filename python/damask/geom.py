import math
from io import StringIO

import numpy as np


class Geom():
  """Geometry definition for grid solvers"""

  def __init__(self,microstructure,size,homogenization=1,comments=[]):
    """New geometry definition from array of microstructures and size"""
    if len(microstructure.shape) != 3:
      raise ValueError('Invalid microstructure shape {}'.format(*microstructure.shape))
    elif microstructure.dtype not in [int,float]:
      raise TypeError('Invalid data type {} for microstructure'.format(microstructure.dtype))
    else:
      self.microstructure = microstructure
    
    if len(size) != 3 or any(np.array(size)<=0):
      raise ValueError('Invalid size {}'.format(*size))
    else:
      self.size = np.array(size)
      
    if not isinstance(homogenization,int) or homogenization < 1:
      raise TypeError('Invalid homogenization {}'.format(homogenization))
    else:
      self.homogenization = homogenization

    if not isinstance(comments,list):
      self.comments = [str(comments)]
    else:
      self.comments = [str(comment) for comment in comments]

  def __repr__(self):
    """Basic information on geometry definition"""
    return 'grid     a b c:      {}\n'.format(' x '.join(map(str,self.get_grid()))) + \
           'size     x y z:      {}\n'.format(' x '.join(map(str,self.size)))       + \
           'homogenization:      {}\n'.format(self.homogenization)                  + \
           '# microstructures:   {}\n'.format(len(np.unique(self.microstructure)))  + \
           'max microstructures: {}\n'.format(np.max(self.microstructure))


  def update(self,microstructure=None,size=None,rescale=False):
    """Updates microstructure and size"""
    grid_old    = self.get_grid()
    size_old    = self.get_size()
    unique_old  = len(np.unique(self.microstructure))
    max_old     = np.max(self.microstructure)
    
    if size is not None and rescale:
      raise ValueError('Either set size explicitly or rescale automatically')

    if microstructure is not None:
      if len(microstructure.shape) != 3:
        raise ValueError('Invalid microstructure shape {}'.format(*microstructure.shape))
      elif microstructure.dtype not in ['int','float']:
        raise TypeError('Invalid data type {} for microstructure'.format(microstructure.dtype))
      else:
        self.microstructure = microstructure
    
    if size is not None:
      if len(size) != 3 or any(np.array(size)<=0):
        raise ValueError('Invalid size {}'.format(*size))
      else:
        self.size = np.array(size)
    
    if rescale:
      self.size = self.size * self.get_grid()/grid_old
    
    message = ''
    if np.any(grid_old != self.get_grid()):
      message += 'grid     a b c:      {}\n'.format(' x '.join(map(str,self.get_grid())))
    if np.any(size_old != self.size):
      message += 'size     x y z:      {}\n'.format(' x '.join(map(str,self.size)))
    if unique_old != len(np.unique(self.microstructure)):
      message += '# microstructures:   {}\n'.format(len(np.unique(self.microstructure)))
    if max_old != np.max(self.microstructure):
      message += 'max microstructures: {}\n'.format(np.max(self.microstructure))
    
    if message != '': return message


  def add_comment(self,comment):
    if not isinstance(comment,list):
      self.comments = [str(comment)] + self.comments
    else:
      self.comments = [str(c) for c in comment] + self.comments

  def set_microstructure(self,microstructure):
    self.microstructure = np.copy(microstructure)

  def set_size(self,size):
    self.size = np.array(size)
    
    
  def get_microstructure(self):
    return np.copy(self.microstructure)

  def get_size(self):
    return np.copy(self.size)

  def get_grid(self):
    return np.array(self.microstructure.shape)

  def get_homogenization(self):
    return self.homogenization

  def get_comments(self):
    return self.comments[:]

  @classmethod
  def from_file(cls,fname):
    """Reads from *.geom file"""
    if isinstance(fname,str):
      f = open(fname)
      header_length,keyword = f.readline().split()
      if not keyword.startswith('head') or int(header_length) < 3:
        raise TypeError('Header length information missing or invalid')
      comments_old = [f.readline() for i in range(int(header_length))]
    else:
      fname.seek(0)
      header_length,keyword = fname.readline().split()
      if not keyword.startswith('head') or int(header_length) < 3:
        raise TypeError('Header length information missing or invalid')
      comments_old = [fname.readline() for i in range(int(header_length))]

    comments = [] 
    for i,line in enumerate(comments_old):
      if line.lower().strip().startswith('grid'):
        grid = np.array([int(line.split()[j]) for j in [2,4,6]])                                   # assume correct order (a,b,c)
      elif line.lower().strip().startswith('size'):
        size = np.array([float(line.split()[j]) for j in [2,4,6]])                                 # assume correct order (x,y,z)
      elif line.lower().strip().startswith('homogenization'):
        homogenization = int(line.split()[1])
      else:
        comments.append(line.rstrip().strip())

    if isinstance(fname,str):
      raw = f.readlines()
      f.close()
    else:
      raw = fname.readlines()

    microstructure = np.empty(grid.prod())                                                          # initialize as flat array
    i = 0
    for line in raw:
      items = line.split()
      if len(items) == 3:
        if   items[1].lower() == 'of':
          items = np.ones(int(items[0]))*float(items[2])
        elif items[1].lower() == 'to':
          items = np.linspace(int(items[0]),int(items[2]),
                              abs(int(items[2])-int(items[0]))+1,dtype=float)
        else:                          items = list(map(float,items))
      else:                            items = list(map(float,items))
      
      microstructure[i:i+len(items)] = items
      i += len(items)
    
    if i != grid.prod():
      raise TypeError('Invalid file: expected {} entries,found {}'.format(grid.prod(),i))
    
    microstructure = microstructure.reshape(grid,order='F')
    
    if '.' in raw[0]:                                                                               # contains float values
      pass
    else:                                                                                           # assume int values
      microstructure = microstructure.astype('int')
    
    return cls(microstructure.reshape(grid),size,homogenization,comments)

  def to_file(self,fname):
    """Saves to file"""
    grid = self.get_grid()
    header =  ['{} header'.format(len(self.comments)+3)]
    header += self.comments
    header.append('grid a {} b {} c {}'.format(*grid))
    header.append('size x {} y {} z {}'.format(*self.size))
    header.append('homogenization {}'.format(self.get_homogenization()))

    if self.microstructure.dtype == 'int':
      format_string='%{}i'.format(int(math.floor(math.log10(self.microstructure.max())+1)))
    else:
      format_string='%.18e'
    np.savetxt(fname, self.microstructure.reshape([grid[0],np.prod(grid[1:])],order='F').T,
               header='\n'.join(header), fmt=format_string, comments='')
               
  def show(self):
    """Show raw content (as in file)"""
    f=StringIO()
    self.to_file(f)
    f.seek(0)
    return ''.join(f.readlines())
