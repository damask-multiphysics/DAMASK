from queue import Queue
import re
import glob

import h5py
import numpy as np

from . import util
from . import version

# ------------------------------------------------------------------
class DADF5():
  """
  Read and write to DADF5 files.
  
  DADF5 files contain DAMASK results.
  """
  
# ------------------------------------------------------------------
  def __init__(self,filename):
    """
    Opens an existing DADF5 file.
    
    Parameters
    ----------
    filename : str
        name of the DADF5 file to be openend.
    
    """
    with h5py.File(filename,'r') as f:
      
      if f.attrs['DADF5-major'] != 0 or f.attrs['DADF5-minor'] != 2:
        raise TypeError('Unsupported DADF5 version {} '.format(f.attrs['DADF5-version']))
    
      self.structured = 'grid' in f['geometry'].attrs.keys()
    
      if self.structured:
        self.grid = f['geometry'].attrs['grid']
        self.size = f['geometry'].attrs['size']
        
      r=re.compile('inc[0-9]+')
      self.increments = [i for i in f.keys() if r.match(i)]
      self.times      = [round(f[i].attrs['time/s'],12) for i in self.increments]

      self.Nmaterialpoints, self.Nconstituents =   np.shape(f['mapping/cellResults/constituent'])
      self.materialpoints  = [m.decode() for m in np.unique(f['mapping/cellResults/materialpoint']['Name'])]
      self.constituents    = [c.decode() for c in np.unique(f['mapping/cellResults/constituent']  ['Name'])]

      self.con_physics  = []
      for c in self.constituents:
        self.con_physics += f['/'.join([self.increments[0],'constituent',c])].keys()
      self.con_physics = list(set(self.con_physics))                                          # make unique

      self.mat_physics = []
      for m in self.materialpoints:
        self.mat_physics += f['/'.join([self.increments[0],'materialpoint',m])].keys()
      self.mat_physics = list(set(self.mat_physics))                                          # make unique

    self.visible= {'increments':     self.increments,
                   'constituents':   self.constituents,
                   'materialpoints': self.materialpoints,
                   'constituent':    range(self.Nconstituents),                               # ToDo: stupid naming
                   'con_physics':    self.con_physics,
                   'mat_physics':    self.mat_physics}
                  
    self.filename   = filename


  def __manage_visible(self,datasets,what,action):
    """
    Manages the visibility of the groups.
    
    Parameters
    ----------
    datasets : list of str or Boolean
      name of datasets as list, supports ? and * wildcards.
      True is equivalent to [*], False is equivalent to []
    what : str
      attribute to change (must be in self.visible) 
    action : str
      select from 'set', 'add', and 'del' 

    """
    # allow True/False and string arguments
    if datasets is True:
      datasets = ['*']
    elif datasets is False:
      datasets = []
    choice = [datasets] if isinstance(datasets,str) else datasets
    
    valid = [e for e_ in [glob.fnmatch.filter(getattr(self,what),s) for s in choice] for e in e_]
    existing = set(self.visible[what])
    
    if   action == 'set':
      self.visible[what] = valid
    elif action == 'add':
      self.visible[what] = list(existing.union(valid))
    elif action == 'del':
      self.visible[what] = list(existing.difference_update(valid))
      
  
  def __time_to_inc(self,start,end):
    selected = []
    for i,time in enumerate(self.times):
      if start <= time < end:
        selected.append(self.increments[i])
    return selected


  def set_by_time(self,start,end):
    """
    Sets active time increments based on start and end time.

    Parameters
    ----------
    start : float
      start time (included)
    end : float
      end time (exclcuded)

    """
    self.__manage_visible(self.__time_to_inc(start,end),'increments','set')


  def add_by_time(self,start,end):
    """
    Adds to active time increments based on start and end time.

    Parameters
    ----------
    start : float
      start time (included)
    end : float
      end time (exclcuded)

    """
    self.__manage_visible(self.__time_to_inc(start,end),'increments','add')


  def del_by_time(self,start,end):
    """
    Delets from active time increments based on start and end time.

    Parameters
    ----------
    start : float
      start time (included)
    end : float
      end time (exclcuded)

    """
    self.__manage_visible(self.__time_to_inc(start,end),'increments','del')


  def iter_visible(self,what):
    """
    Iterates over visible items by setting each one visible.

    Parameters
    ----------
    what : str
      attribute to change (must be in self.visible) 

    """
    datasets = self.visible[what]
    last_datasets = datasets.copy()
    for dataset in datasets:
      if last_datasets != self.visible[what]:
        self.__manage_visible(datasets,what,'set')
        raise Exception
      self.__manage_visible(dataset,what,'set')
      last_datasets = self.visible[what]
      yield dataset
    self.__manage_visible(datasets,what,'set')
    
  
  def set_visible(self,what,datasets):
    """
    Sets active groups.
    
    Parameters
    ----------
    datasets : list of str or Boolean
      name of datasets as list, supports ? and * wildcards.
      True is equivalent to [*], False is equivalent to []
    what : str
      attribute to change (must be in self.visible) 

    """
    self.__manage_visible(datasets,what,'set')


  def add_visible(self,what,datasets):
    """
    Adds to active groups.
    
    Parameters
    ----------
    datasets : list of str or Boolean
      name of datasets as list, supports ? and * wildcards.
      True is equivalent to [*], False is equivalent to []
    what : str
      attribute to change (must be in self.visible) 

    """
    self.__manage_visible(datasets,what,'add')


  def del_visible(self,what,datasets):
    """
    Removes from active groupse.
    
    Parameters
    ----------
    datasets : list of str or Boolean
      name of datasets as list, supports ? and * wildcards.
      True is equivalent to [*], False is equivalent to []
    what : str
      attribute to change (must be in self.visible) 

    """
    self.__manage_visible(datasets,what,'del')


  def groups_with_datasets(self,datasets):
    """
    Get groups that contain all requested datasets.
    
    Only groups within inc?????/constituent/*_*/* inc?????/materialpoint/*_*/* 
    are considered as they contain the data.
    Single strings will be treated as list with one entry.
    
    Wild card matching is allowed, but the number of arguments need to fit.
    
    Parameters
    ----------
      datasets : iterable or str or boolean
    
    Examples
    --------
      datasets = False matches no group
      datasets = True matches all groups
      datasets = ['F','P'] matches a group with ['F','P','sigma']
      datasets = ['*','P'] matches a group with ['F','P']
      datasets = ['*'] does not match a group with ['F','P','sigma']
      datasets = ['*','*'] does not match a group with ['F','P','sigma']
      datasets = ['*','*','*'] matches a group with ['F','P','sigma']
      
    """
    if datasets is False: return []
    sets = [datasets] if isinstance(datasets,str) else datasets

    groups = []
    
    with h5py.File(self.filename,'r') as f:
      for i in self.iter_visible('increments'):
        for o,p in zip(['constituents','materialpoints'],['con_physics','mat_physics']):
          for oo in self.iter_visible(o):
            for pp in self.iter_visible(p):
              group = '/'.join([i,o[:-1],oo,pp])                              # o[:-1]: plural/singular issue
              if sets is True:
                groups.append(group)
              else:
                match = [e for e_ in [glob.fnmatch.filter(f[group].keys(),s) for s in sets] for e in e_]
                if len(set(match)) == len(sets) : groups.append(group)
    return groups


  def list_data(self):
    """Gives information on all active datasets in the file."""
    message = ''
    with h5py.File(self.filename,'r') as f:
      for i in self.iter_visible('increments'):
        message+='\n{}\n'.format(i)
        for o,p in zip(['constituents','materialpoints'],['con_physics','mat_physics']):
          for oo in self.iter_visible(o):
            message+='  {}\n'.format(oo)
            for pp in self.iter_visible(p):
              message+='    {}\n'.format(pp)
              group = '/'.join([i,o[:-1],oo,pp])                              # o[:-1]: plural/singular issue
              for d in f[group].keys():
                try:
                  message+='      {} ({})\n'.format(d,f['/'.join([group,d])].attrs['Description'].decode())
                except KeyError:
                  pass
    return message


  def get_dataset_location(self,label):
    """Returns the location of all active datasets with given label."""
    path = []
    with h5py.File(self.filename,'r') as f:
      for i in self.iter_visible('increments'):        
        for o,p in zip(['constituents','materialpoints'],['con_physics','mat_physics']):
          for oo in self.iter_visible(o):
            for pp in self.iter_visible(p):
              k = '/'.join([i,o[:-1],oo,pp,label])
              try:
                f[k]
                path.append(k)
              except KeyError as e:
                print('unable to locate constituents dataset: '+ str(e))
    return path
    
    
  def get_constituent_ID(self,c=0):
    """Pointwise constituent ID."""
    with h5py.File(self.filename,'r') as f:
      names = f['/mapping/cellResults/constituent']['Name'][:,c].astype('str')
    return np.array([int(n.split('_')[0]) for n in names.tolist()],dtype=np.int32)


  def get_crystal_structure(self):                                                                  # ToDo: extension to multi constituents/phase
    """Info about the crystal structure."""
    with h5py.File(self.filename,'r') as f:
      return f[self.get_dataset_location('orientation')[0]].attrs['Lattice'].astype('str')          # np.bytes_ to string


  def read_dataset(self,path,c):
    """
    Dataset for all points/cells.
    
    If more than one path is given, the dataset is composed of the individual contributions
    """
    with h5py.File(self.filename,'r') as f:
      shape = (self.Nmaterialpoints,) + np.shape(f[path[0]])[1:]
      if len(shape) == 1: shape = shape +(1,)
      dataset = np.full(shape,np.nan,dtype=np.dtype(f[path[0]]))
      for pa in path:
        label   = pa.split('/')[2]
        
        p = np.where(f['mapping/cellResults/constituent'][:,c]['Name'] == str.encode(label))[0]
        if len(p)>0:
          u = (f['mapping/cellResults/constituent'][p,c]['Position'])
          a = np.array(f[pa])
          if len(a.shape) == 1:
            a=a.reshape([a.shape[0],1])
          dataset[p,:] = a[u,:]

        p = np.where(f['mapping/cellResults/materialpoint']['Name'] == str.encode(label))[0]
        if len(p)>0:
          u = (f['mapping/cellResults/materialpoint'][p.tolist()]['Position'])
          a = np.array(f[pa])
          if len(a.shape) == 1:
            a=a.reshape([a.shape[0],1])
          dataset[p,:] = a[u,:]

    return dataset


  def cell_coordinates(self):
    """Initial coordinates of the cell centers."""
    if self.structured:
      delta = self.size/self.grid*0.5
      z, y, x = np.meshgrid(np.linspace(delta[2],self.size[2]-delta[2],self.grid[2]),
                            np.linspace(delta[1],self.size[1]-delta[1],self.grid[1]),
                            np.linspace(delta[0],self.size[0]-delta[0],self.grid[0]),
                           )
      return np.concatenate((x[:,:,:,None],y[:,:,:,None],y[:,:,:,None]),axis = 3).reshape([np.product(self.grid),3])
    else:
      with h5py.File(self.filename,'r') as f:
        return f['geometry/x_c'][()]


  def add_Cauchy(self,P='P',F='F'):
    """
    Adds Cauchy stress calculated from 1st Piola-Kirchhoff stress and deformation gradient.
    
    Resulting tensor is symmetrized as the Cauchy stress should be symmetric.
    """
    def Cauchy(F,P):
      sigma = np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F['data']),P['data'],F['data'])
      sigma = (sigma + np.transpose(sigma,(0,2,1)))*0.5                                             # enforce symmetry
      return  {
               'data'  : sigma,
               'label' : 'sigma',
               'meta' : {
                        'Unit' :        P['meta']['Unit'],
                        'Description' : 'Cauchy stress calculated from {} ({}) '.format(P['label'],P['meta']['Description'])+\
                                        'and deformation gradient {} ({})'.format(F['label'],F['meta']['Description']),
                        'Creator' :     'dadf5.py:add_Cauchy v{}'.format(version)
                        }
               }
      
    requested = [{'label':F,'arg':'F'},
                 {'label':P,'arg':'P'} ]
    
    self.__add_generic_pointwise(Cauchy,requested)


  def add_Mises(self,x):
    """Adds the equivalent Mises stress or strain of a tensor."""
    def Mises(x):
      
      if x['meta']['Unit'] == b'Pa': #ToDo: Should we use this? Then add_Cauchy  and add_strain_tensors also should perform sanity checks
        factor = 3.0/2.0
        t = 'stress'
      elif x['meta']['Unit'] == b'1':
        factor = 2.0/3.0
        t = 'strain'
      else:
        print(x['meta']['Unit'])
        raise ValueError
      
      d = x['data']
      dev = d - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[d.shape[0],3,3]),np.trace(d,axis1=1,axis2=2)/3.0)
      #dev_sym = (dev + np.einsum('ikj',dev))*0.5 # ToDo: this is not needed (only if the input is not symmetric, but then the whole concept breaks down)

      return  {
               'data' :  np.sqrt(np.einsum('ijk->i',dev**2)*factor), 
               'label' : '{}_vM'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Mises equivalent {} of {} ({})'.format(t,x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_Mises_stress v{}'.format(version)
                        }
               }
      
    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(Mises,requested)


  def add_norm(self,x,ord=None):
    """
    Adds norm of vector or tensor.
    
    See numpy.linalg.norm manual for details.
    """
    def norm(x,ord):

      o = ord
      if len(x['data'].shape) == 2:
        axis = 1
        t = 'vector'
        if o is None: o = 2
      elif len(x['data'].shape) == 3:
        axis = (1,2)
        t = 'tensor'
        if o is None: o = 'fro'
      else:
        raise ValueError

      return  {
               'data' : np.linalg.norm(x['data'],ord=o,axis=axis,keepdims=True),
               'label' : '|{}|_{}'.format(x['label'],o),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : '{}-Norm of {} {} ({})'.format(ord,t,x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_norm v{}'.format(version)
                        }
               }

    requested = [{'label':x,'arg':'x'}]

    self.__add_generic_pointwise(norm,requested,{'ord':ord})
    
    
  def add_absolute(self,x):
    """Adds absolute value."""
    def absolute(x):

      return  {
               'data' : np.abs(x['data']),
               'label' : '|{}|'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Absolute value of {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_abs v{}'.format(version)
                        }
               }

    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(absolute,requested)


  def add_determinant(self,x):
    """Adds the determinant component of a tensor."""
    def determinant(x):
      
      return  {
               'data' : np.linalg.det(x['data']),
               'label' : 'det({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Determinant of tensor {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_determinant v{}'.format(version)
                        }
               }
      
    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(determinant,requested)


  def add_spherical(self,x):
    """Adds the spherical component of a tensor."""
    def spherical(x):
      
      if not np.all(np.array(x['data'].shape[1:]) == np.array([3,3])):
        raise ValueError

      return  {
               'data' : np.trace(x['data'],axis1=1,axis2=2)/3.0,
               'label' : 'sph({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Spherical component of tensor {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_spherical v{}'.format(version)
                        }
               }

    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(spherical,requested)


  def add_deviator(self,x):
    """Adds the deviator of a tensor."""
    def deviator(x):
      d = x['data']
      
      if not np.all(np.array(d.shape[1:]) == np.array([3,3])):
        raise ValueError
      
      return  {
               'data' :  d - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[d.shape[0],3,3]),np.trace(d,axis1=1,axis2=2)/3.0), 
               'label' : 'dev({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Deviator of tensor {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_deviator v{}'.format(version)
                        }
               }
      
    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(deviator,requested)
    

  def add_calculation(self,formula,label,unit='n/a',description=None,vectorized=True):
    """
    General formula.
    
    Works currently only for vectorized expressions
    
    """
    if vectorized is not True:
      raise NotImplementedError
    
    def calculation(**kwargs):
      
      formula = kwargs['formula']
      for d in re.findall(r'#(.*?)#',formula):
        formula = re.sub('#{}#'.format(d),"kwargs['{}']['data']".format(d),formula)
        
      return  {
               'data' : eval(formula), 
               'label' : kwargs['label'],
               'meta' : {
                        'Unit' :        kwargs['unit'],
                        'Description' : '{}'.format(kwargs['description']),
                        'Creator' :     'dadf5.py:add_calculation v{}'.format(version)
                        }
               }
    
    requested    = [{'label':d,'arg':d} for d in re.findall(r'#(.*?)#',formula)]         # datasets used in the formula
    pass_through = {'formula':formula,'label':label,'unit':unit,'description':description} 
    
    self.__add_generic_pointwise(calculation,requested,pass_through)


  def add_strain_tensor(self,t,ord,defgrad='F'): #ToDo: Use t and ord
    """
    Adds the a strain tensor.
    
    Albrecht Bertram: Elasticity and Plasticity of Large Deformations An Introduction (3rd Edition, 2012), p. 102.
    """
    def strain_tensor(defgrad,t,ord):
      
      operator = { 
                   'V#ln':   lambda V: np.log(V),
                   'U#ln':   lambda U: np.log(U),
                   'V#Biot': lambda V: np.broadcast_to(np.ones(3),[V.shape[0],3]) - 1.0/V,
                   'U#Biot': lambda U: U - np.broadcast_to(np.ones(3),[U.shape[0],3]),
                   'V#Green':lambda V: np.broadcast_to(np.ones(3),[V.shape[0],3]) - 1.0/V**2,
                   'U#Green':lambda U: U**2 - np.broadcast_to(np.ones(3),[U.shape[0],3]), 
                 }

      (U,S,Vh) = np.linalg.svd(defgrad['data'])                                                             # singular value decomposition
      R_inv    = np.transpose(np.matmul(U,Vh),(0,2,1))                                                      # transposed rotation of polar decomposition
      U        = np.matmul(R_inv,defgrad['data'])                                                           # F = RU
      (D,V)    = np.linalg.eigh((U+np.transpose(U,(0,2,1)))*.5)                                             # eigen decomposition (of symmetric(ed) matrix)

      neg      = np.where(D < 0.0)                                                                          # find negative eigenvalues ...
      D[neg[0],neg[1]]   = D[neg[0],neg[1]]* -1                                                             # ... flip value ...
      V[neg[0],:,neg[1]] = V[neg[0],:,neg[1]]* -1                                                           # ... and vector
      
      d = operator['V#ln'](D)
      a = np.matmul(V,np.einsum('ij,ikj->ijk',d,V))
      
      return  {
               'data' : a, 
               'label' : 'ln(V)({})'.format(defgrad['label']),
               'meta' : {
                        'Unit' :        defgrad['meta']['Unit'],
                        'Description' : 'Strain tensor ln(V){} ({})'.format(defgrad['label'],defgrad['meta']['Description']),
                        'Creator' :     'dadf5.py:add_deviator v{}'.format(version)
                        }
               }

    requested = [{'label':defgrad,'arg':'defgrad'}]
    
    self.__add_generic_pointwise(strain_tensor,requested,{'t':t,'ord':ord})


  def __add_generic_pointwise(self,func,datasets_requested,extra_args={}):
    """
    General function to add pointwise data.
    
    Parameters
    ----------
      func : function
        Function that calculates a new dataset from one or more datasets per HDF5 group.
      datasets_requested : list of dictionaries
        Details of the datasets to be used: label (in HDF5 file) and arg (argument to which the data is parsed in func).
      extra_args : dictionary, optional
        Any extra arguments parsed to func.

    """
    def job(args):
      """Call function with input data + extra arguments, returns results + group."""
      args['results'].put({**args['func'](**args['in']),'group':args['group']})
    

    N_threads = 1 # ToDo: should be a parameter

    results = Queue(N_threads)
    pool    = util.ThreadPool(N_threads)
    N_added = N_threads + 1
    
    todo = []
    # ToDo: It would be more memory efficient to read only from file when required, i.e. do to it in pool.add_task
    for group in self.groups_with_datasets([d['label'] for d in datasets_requested]):
      with h5py.File(self.filename,'r') as f:
        datasets_in = {}
        for d in datasets_requested:
          loc  = f[group+'/'+d['label']]
          data = loc[()]
          meta = {k:loc.attrs[k] for k in loc.attrs.keys()}
          datasets_in[d['arg']] = {'data': data, 'meta' : meta, 'label' : d['label']}

      todo.append({'in':{**datasets_in,**extra_args},'func':func,'group':group,'results':results})
    
    pool.map(job, todo[:N_added])                                                                   # initialize

    N_not_calculated = len(todo)
    while N_not_calculated > 0:    
      result = results.get()
      with h5py.File(self.filename,'a') as f:                                                       # write to file
        dataset_out = f[result['group']].create_dataset(result['label'],data=result['data'])
        for k in result['meta'].keys():
          dataset_out.attrs[k] = result['meta'][k]
        N_not_calculated-=1
      
      if N_added < len(todo):                                                                       # add more jobs
        pool.add_task(job,todo[N_added])
        N_added +=1

    pool.wait_completion()
