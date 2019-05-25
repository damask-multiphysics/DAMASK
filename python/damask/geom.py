import numpy as np
import math
from io import StringIO

class Geom():
  """Geometry definition for grid solvers"""

  def __init__(self,size,microstructure,homogenization=1,comments=[]):
    """New geometry definition from array of microstructures and size"""
    if len(size) != 3 or any(np.array(size)<=0):
      raise ValueError('invalid size')
    else:
      self.size = np.array(size)

    if len(microstructure.shape) != 3:
      raise ValueError('invalid microstructure')
    else:
      self.microstructure = microstructure

    if not isinstance(homogenization,int) or homogenization < 1:
      raise ValueError('invalid homogenization')
    else:
      self.homogenization = homogenization

    if not isinstance(comments,list):
      self.comments = [str(comments)]
    else:
      self.comments = [str(comment) for comment in comments]

  def __repr__(self):
    """Basic information on geometry definition"""
    return 'grid     a b c:      {}\n'.format(' x '.join(map(str,self.get_grid()))) + \
           'size     x y z:      {}\n'.format(' x '.join(map(str,self.get_size()))) + \
           'homogenization:      {}\n'.format(self.get_homogenization())            + \
           '# microstructures:   {}\n'.format(len(np.unique(self.microstructure)))  + \
           'max microstructures: {}\n'.format(np.max(self.microstructure))
    
  def add_comment(self,comment):
    if not isinstance(comment,list):
      self.comments += [str(comment)]
    else:
      self.comments += [str(c) for c in comment]

  def set_size(self,size):
    self.size = np.array(size)
    
  def get_size(self):
    return self.size

  def get_grid(self):
    return np.array(self.microstructure.shape)

  def get_homogenization(self):
    return self.homogenization

  @classmethod
  def from_file(cls,fname):
    if isinstance(fname,str):
      with open(fname) as f:
        header_length = int(f.readline().split()[0])
        comments_old = [f.readline() for i in range(header_length)]
    else:
      fname.seek(0)
      header_length = int(fname.readline().split()[0])
      comments_old = [fname.readline() for i in range(header_length)]

    comments = [] 
    for i,line in enumerate(comments_old):
      if line.lower().strip().startswith('grid'):
        grid = np.array([int(line.split()[j]) for j in [2,4,6]])         # assume correct order (a,b,c)
      elif line.lower().strip().startswith('size'):
        size = np.array([float(line.split()[j]) for j in [2,4,6]])       # assume correct order (x,y,z)
      elif line.lower().strip().startswith('homogenization'):
        homogenization = int(line.split()[1])
      else:
        comments.append(line.rstrip().strip())

    if isinstance(fname,str):
      with open(fname) as f:
        raw = f.readlines()[header_length+1:]
    else:
      raw = fname.readlines()

    microstructure = np.empty(grid.prod())                                                         # initialize as flat array
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
    
    microstructure = microstructure.reshape(grid,order='F')
    
    if np.any(np.mod(microstructure.flatten(),1)!=0.0):
      pass
    else:
      microstructure = microstructure.astype('int')
    
    return cls(size,microstructure.reshape(grid),homogenization,comments)

  def to_file(self,fname):
    grid = self.get_grid()
    header =  ['{} header'.format(len(self.comments)+3)]
    header += self.comments
    header.append('grid a {} b {} c {}'.format(*grid))
    header.append('size x {} y {} z {}'.format(*self.get_size()))
    header.append('homogenization {}'.format(self.get_homogenization()))

    if self.microstructure.dtype == 'int':
      format_string='%{}i'.format(int(math.floor(math.log10(self.microstructure.max())+1)))
    else:
      format_string='%.18e'
    np.savetxt(fname, self.microstructure.reshape([grid[0],np.prod(grid[1:])],order='F').T,
               header='\n'.join(header), fmt=format_string, comments='')
               
  def show(self):
    f=StringIO()
    self.to_file(f)
    f.seek(0)
    return ''.join(f.readlines())
