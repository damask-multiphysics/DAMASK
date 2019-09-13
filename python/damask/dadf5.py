from queue import Queue
import re

import h5py
import numpy as np

from . import util

# ------------------------------------------------------------------
class DADF5():
  """
  Read and write to DADF5 files.
  
  DADF5 files contain DAMASK results.
  """
  
# ------------------------------------------------------------------
  def __init__(self,
               filename,
               mode     = 'a',
              ):
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
      self.increments = [{'inc':  int(u[3:]),
                          'time': round(f[u].attrs['time/s'],12),
                          }   for u in f.keys() if r.match(u)]
      
      self.Nmaterialpoints, self.Nconstituents =   np.shape(f['mapping/cellResults/constituent'])
      self.materialpoints  = [m.decode() for m in np.unique(f['mapping/cellResults/materialpoint']['Name'])]
      self.constituents    = [c.decode() for c in np.unique(f['mapping/cellResults/constituent']  ['Name'])]


      self.c_output_types  = []
      for c in self.constituents:
        for o in f['inc{:05}/constituent/{}'.format(self.increments[0]['inc'],c)].keys():
          self.c_output_types.append(o)
      self.c_output_types = list(set(self.c_output_types))                                          # make unique

      self.m_output_types = []
      for m in self.materialpoints:
        for o in f['inc{:05}/materialpoint/{}'.format(self.increments[0]['inc'],m)].keys():
          self.m_output_types.append(o)
      self.m_output_types = list(set(self.m_output_types))                                          # make unique

    #self.on_air
    self.active= {'increments':     self.increments,                                                # ToDo:simplify, activity only positions that translate into (no complex types)
                  'constituents':   self.constituents,
                  'materialpoints': self.materialpoints,
                  'constituent':    range(self.Nconstituents),                                      # ToDo: stupid naming
                  'c_output_types': self.c_output_types,
                  'm_output_types': self.m_output_types}

# ToDo: store increments, select icrements (trivial), position, and time

    self.filename   = filename
    self.mode       = mode


  def get_groups(self,l): #group_with_data(datasets)
    """
    Get groups that contain all requested datasets.
    
    Parameters
    ----------
    l : list of str
        Names of datasets that need to be located in the group.

    """    
    groups = []
    if type(l) is not list:
     raise TypeError('Candidates should be given as a list')
    with h5py.File(self.filename,'r') as f:
      for g in self.get_active_groups():
        if set(l).issubset(f[g].keys()): groups.append(g)
    return groups


  def get_active_groups(self): # rename: get_groups needed? merge with datasets and have [] and ['*']
    """
    Get groups that are currently considered for evaluation. 
    """
    groups = []
    for i,x in enumerate(self.active['increments']):
      group_inc = 'inc{:05}'.format(self.active['increments'][i]['inc']) #ToDo: Merge path only once at the end '/'.join(listE)
      for c in self.active['constituents']:
        group_constituent = group_inc+'/constituent/'+c
        for t in self.active['c_output_types']:
          group_output_types = group_constituent+'/'+t
          groups.append(group_output_types)
      for m in self.active['materialpoints']:
        group_materialpoint = group_inc+'/materialpoint/'+m
        for t in self.active['m_output_types']:
          group_output_types = group_materialpoint+'/'+t
          groups.append(group_output_types)
    return groups
    

  def list_data(self): # print_datasets and have [] and ['*'], loop over all increment, soll auf anderen basieren (get groups with sternchen)
    """Shows information on all active datasets in the file."""
    with h5py.File(self.filename,'r') as f:
      group_inc = 'inc{:05}'.format(self.active['increments'][0]['inc']) #ToDo: Merge path only once at the end '/'.join(listE)
      for c in self.active['constituents']:
        print('\n'+c)
        group_constituent = group_inc+'/constituent/'+c
        for t in self.active['c_output_types']:
          print('  {}'.format(t))
          group_output_types = group_constituent+'/'+t 
          try:
            for x in f[group_output_types].keys():
              print('    {} ({})'.format(x,f[group_output_types+'/'+x].attrs['Description'].decode()))
          except KeyError:
            pass
      for m in self.active['materialpoints']:
        group_materialpoint = group_inc+'/materialpoint/'+m
        for t in self.active['m_output_types']:
          print('  {}'.format(t))
          group_output_types = group_materialpoint+'/'+t
          try:
            for x in f[group_output_types].keys():
              print('    {} ({})'.format(x,f[group_output_types+'/'+x].attrs['Description'].decode()))
          except KeyError:
            pass
    

  def get_dataset_location(self,label): # names
    """Returns the location of all active datasets with given label.""" #ToDo: Merge path only once at the end '/'.join(listE)
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
            except KeyError as e:
              print('unable to locate constituents dataset: '+ str(e))
       
        for m in self.active['materialpoints']:
          group_materialpoint = group_inc+'/materialpoint/'+m
          for t in self.active['m_output_types']:
            try:
              f[group_materialpoint+'/'+t+'/'+label]
              path.append(group_materialpoint+'/'+t+'/'+label)
            except KeyError as e:
              print('unable to locate materialpoints dataset: '+ str(e))
              
    return path
    
    
  def read_dataset(self,path,c):
    """
    Dataset for all points/cells.
    
    If more than one path is given, the dataset is composed of the individual contributions
    """
    with h5py.File(self.filename,'r') as f:
      shape = (self.Nmaterialpoints,) + np.shape(f[path[0]])[1:]
      if len(shape) == 1: shape = shape +(1,)
      dataset = np.full(shape,np.nan)
      for pa in path:
        label   = pa.split('/')[2]
        try:
          p = np.where(f['mapping/cellResults/constituent'][:,c]['Name'] == str.encode(label))[0]
          u = (f['mapping/cellResults/constituent'][p,c]['Position'])
          a = np.array(f[pa])
          if len(a.shape) == 1:
            a=a.reshape([a.shape[0],1])
          dataset[p,:] = a[u,:]
        except KeyError as e:
          print('unable to read constituent: '+ str(e))
        try:
          p = np.where(f['mapping/cellResults/materialpoint']['Name'] == str.encode(label))[0]
          u = (f['mapping/cellResults/materialpoint'][p.tolist()]['Position'])
          a = np.array(f[pa])
          if len(a.shape) == 1:
            a=a.reshape([a.shape[0],1])
          dataset[p,:] = a[u,:]
        except KeyError as e:
          print('unable to read materialpoint: '+ str(e))

    return dataset


  def add_Cauchy(self,P='P',F='F'):
    """
    Adds Cauchy stress calculated from 1st Piola-Kirchhoff stress and deformation gradient.
    
    Todo
    ----
    The einsum formula is completely untested!
    
    """
    def Cauchy(F,P):
      return  {
               'data' : np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F['data']),F['data'],P['data']),
               'label' : 'sigma',
               'meta' : {
                        'Unit' :        P['meta']['Unit'],
                        'Description' : 'Cauchy stress calculated from {} ({}) '.format(P['label'],P['meta']['Description'])+\
                                        'and deformation gradient {} ({})'.format(F['label'],P['meta']['Description']),
                        'Creator' :     'dadf5.py:add_Cauchy vXXXXX'
                        }
               }
      
    requested = [{'label':F,'arg':'F'},
                 {'label':P,'arg':'P'} ]
    
    self.__add_generic_pointwise(Cauchy,requested)


  def add_Mises(self,x):
    """Adds the equivalent Mises stres of a tensor."""
    def deviator(x):
      
      if x['meta']['Unit'] == 'Pa':
        factor = 3.0/2.0
      elif x['meta']['Unit'] == '-':
        factor = 2.0/3.0
      else:
        ValueError
      
      d = x['data']
      dev = d - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[d.shape[0],3,3]),np.trace(d,axis1=1,axis2=2)/3.0)
      dev_sym = (dev + np.einsum('ikj',dev))*0.5

      return  {
               'data' :  np.sqrt(np.einsum('ijk->i',dev_sym**2)*factor), 
               'label' : 'dev({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Mises equivalent stress of {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_Mises_stress vXXXXX'
                        }
               }
      
    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(deviator,requested)


  def add_norm(self,x,ord=None):
    """
    Adds norm of vector or tensor or magnitude of a scalar.
    See numpy.linalg.norm manual for details.
    """
    def norm(x,ord):

      if   len(x['data'].shape) == 1:
        axis = 0
        t = 'scalar'
      elif len(x['data'].shape) == 2:
        axis = 1
        t = 'vector'
      elif len(x['data'].shape) == 3:
        axis = (1,2)
        t = 'tensor'
      else:
        raise ValueError

      return  {
               'data' : np.linalg.norm(x['data'],ord=ord,axis=axis,keepdims=True),
               'label' : 'norm({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Norm of {} {} ({})'.format(t,x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_norm vXXXXX'
                        }
               }

    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(norm,requested,{'ord':ord})


  def add_determinant(self,x):
    """Adds the determinant component of a tensor."""
    def determinant(x):
      
      return  {
               'data' : np.linalg.det(x['data']),
               'label' : 'det({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Determinant of tensor {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_determinant vXXXXX'
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
                        'Creator' :     'dadf5.py:add_spherical vXXXXX'
                        }
               }

    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(spherical,requested)


  def add_deviator(self,x):
    """Adds the deviator of a tensor."""
    def deviator(x):
      d = x['data']
      return  {
               'data' :  d - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[d.shape[0],3,3]),np.trace(d,axis1=1,axis2=2)/3.0), 
               'label' : 'dev({})'.format(x['label']),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : 'Deviator of tensor {} ({})'.format(x['label'],x['meta']['Description']),
                        'Creator' :     'dadf5.py:add_deviator vXXXXX'
                        }
               }
      
    requested = [{'label':x,'arg':'x'}]
    
    self.__add_generic_pointwise(deviator,requested)


  def add_strain_tensor(self,t,ord,defgrad='F'):
    """Adds the a strain tensor."""
    def strain_tensor(defgrad,t,ord):
      (U,S,Vh) = np.linalg.svd(defgrad['data'])                                                             # singular value decomposition
      R_inv    = np.einsum('ijk->ikj',np.matmul(U,Vh))                                                      # inverse rotation of polar decomposition
      U        = np.matmul(R_inv,defgrad['data'])                                                           # F = RU
      (D,V)    = np.linalg.eigh((U+np.einsum('ikj',U))*.5)                                                  # eigen decomposition (of symmetric(ed) matrix)

      neg      = np.where(D < 0.0)                                                                          # find negative eigenvalues ...
      D[neg[0],neg[1]]   = D[neg[0],neg[1]]* -1                                                             # ... flip value ...
      V[neg[0],:,neg[1]] = V[neg[0],:,neg[1]]* -1                                                           # ... and vector
      
      d = np.log(D)
      a = np.matmul(V,np.einsum('ij,ikj->ijk',d,V)) # this is wrong ... 
      for j in range(V.shape[0]):                                  # but this is slow ...
        a[j,:,:] = np.dot(V[j,:,:],np.dot(np.diag(d[j,:]),V[j,:,:].T))
      print(np.max(a))

      return  {
               'data' : a, 
               'label' : 'lnV({})'.format(defgrad['label']),
               'meta' : {
                        'Unit' :        defgrad['meta']['Unit'],
                        'Description' : 'Strain tensor {} ({})'.format(defgrad['label'],defgrad['meta']['Description']),
                        'Creator' :     'dadf5.py:add_deviator vXXXXX'
                        }
               }

    requested = [{'label':defgrad,'arg':'defgrad'}]
    
    self.__add_generic_pointwise(strain_tensor,requested,{'t':t,'ord':ord})


  def __add_generic_pointwise(self,func,datasets_requested,extra_args={}):
    """
    General function to add pointwise data.
    """

    def job(args):
      """
      Call function with input data + extra arguments, returns results + group.
      """
      args['results'].put({**args['func'](**args['in']),'group':args['group']})
    

    N_threads = 1 # ToDo: should be a parameter

    results = Queue(N_threads)
    pool    = util.ThreadPool(N_threads)
    N_added = N_threads + 1
    
    todo = []
    # ToDo: It would be more memory efficient to read only from file when required, i.e. do to it in pool.add_task
    for group in self.get_groups([d['label'] for d in datasets_requested]):
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
      with h5py.File(self.filename,self.mode) as f:                                                 # write to file
        dataset_out = f[result['group']].create_dataset(result['label'],data=result['data'])
        for k in result['meta'].keys():
          dataset_out.attrs[k] = result['meta'][k]
        N_not_calculated-=1
      
      if N_added < len(todo):                                                                       # add more jobs
        pool.add_task(job,todo[N_added])
        N_added +=1

    pool.wait_completion()
