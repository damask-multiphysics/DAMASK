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

    self.visible= {'increments':     self.increments, # ToDo:simplify, activity only positions that translate into (no complex types)
                  'constituents':   self.constituents,
                  'materialpoints': self.materialpoints,
                  'constituent':    range(self.Nconstituents),                                      # ToDo: stupid naming
                  'c_output_types': self.c_output_types,
                  'm_output_types': self.m_output_types}
                  
    self.filename   = filename

  def __visible_set(self,output,t,p):
    valid  = set(p)
    choice = [output] if isinstance(output,str) else output
    self.visible[t] = list(valid.intersection(choice))


  def __visible_add(self,output,t,p):
    choice   = [output] if isinstance(output,str) else output
    valid    = set(p).intersection(choice)
    existing = set(self.visible[t])
    self.visible[t] = list(existing.add(valid))


  def __visible_del(self,output,t):
    choice   = [output] if isinstance(output,str) else output
    existing = set(self.visible[t])
    self.visible[t] = list(existing.remove(choice))


  def __visible_iter(self,t):
    a = self.visible[t]
    last_a = a.copy()
    for i in a:
      if last_a != self.visible[t]:
        self.__visible_set(a,t,a)
        raise Exception
      self.__visible_set(i,t,a)
      last_a = self.visible[t]
      yield i
    self.__visible_set(a,t,a)
  
  
  def constituent_output_iter(self):
    return self.__visible_iter('c_output_types')
  
    
  def constituent_output_set(self,output):
    self.__visible_set(output,'c_output_types',self.c_output_types)


  def constituent_output_add(self,output):
    self.__visible_add(output,'c_output_types',self.c_output_types)


  def constituent_output_del(self,output):
    self.__visible_del(output,'c_output_types')
  
    
  def materialpoint_output_iter(self):
    return self.__visible_iter('m_output_types')
  
    
  def materialpoint_output_set(self,output):
    self.__visible_set(output,'m_output_types',self.m_output_types)


  def materialpoint_output_add(self,output):
    self.__visible_add(output,'m_output_types',self.m_output_types)


  def materialpoint_output_del(self,output):
    self.__visible_del(output,'m_output_types')
  
  
  def constituent_iter(self):
    return self.__visible_iter('constituents')
    
    
  def constituent_set(self,output):
    self.__visible_set(output,'constituents',self.constituents)


  def constituent_add(self,output):
    self.__visible_add(output,'constituents',self.constituents)


  def constituent_del(self,output):
    self.__visible_del(output,'constituents')
  
  
  def materialpoint_iter(self):
    return self.__visible_iter('materialpoints')
    

  def materialpoint_set(self,output):
    self.__visible_set(output,'materialpoints',self.materialpoints)


  def materialpoint_add(self,output):
    self.__visible_add(output,'materialpoints',self.materialpoints)


  def materialpoint_del(self,output):
    self.__visible_del(output,'materialpoints')
    


# ToDo: store increments, select icrements (trivial), position, and time



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
    """Get groups that are currently considered for evaluation."""
    groups = []
    for i in self.visible['increments']:
      group_inc = 'inc{:05}'.format(i['inc'])                               #ToDo: Merge path only once at the end '/'.join(listE)
      for c in self.visible['constituents']:
        for t in self.visible['c_output_types']:
          groups.append('/'.join([group_inc,'constituent',c,t]))
      for m in self.visible['materialpoints']:
        for t in self.visible['m_output_types']:
          groups.append('/'.join([group_inc,'materialpoint',m,t]))
    return groups
    

  def list_data(self): # print_datasets and have [] and ['*'], loop over all increment, soll auf anderen basieren (get groups with sternchen)
    """Shows information on all active datasets in the file."""
    with h5py.File(self.filename,'r') as f:
      group_inc = 'inc{:05}'.format(self.visible['increments'][0]['inc']) #ToDo: Merge path only once at the end '/'.join(listE)
      for c in self.visible['constituents']:
        print('\n'+c)
        group_constituent = group_inc+'/constituent/'+c
        for t in self.visible['c_output_types']:
          print('  {}'.format(t))
          group_output_types = group_constituent+'/'+t 
          try:
            for x in f[group_output_types].keys():
              print('    {} ({})'.format(x,f[group_output_types+'/'+x].attrs['Description'].decode()))
          except KeyError:
            pass
      for m in self.visible['materialpoints']:
        group_materialpoint = group_inc+'/materialpoint/'+m
        for t in self.visible['m_output_types']:
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
      for i in self.visible['increments']:
        group_inc = 'inc{:05}'.format(i['inc'])
        
        for c in self.visible['constituents']:
          for t in self.visible['c_output_types']:
            try:
              p = '/'.join([group_inc,'constituent',c,t,label])
              f[p]
              path.append(p)
            except KeyError as e:
              print('unable to locate constituents dataset: '+ str(e))
       
        for m in self.visible['materialpoints']:
          for t in self.visible['m_output_types']:
            try:
              p = '/'.join([group_inc,'materialpoint',m,t,label])
              f[p]
              path.append(p)
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


  def add_Cauchy(self,P='P',F='F'):
    """
    Adds Cauchy stress calculated from 1st Piola-Kirchhoff stress and deformation gradient.
    
    Resulting tensor is symmetrized as the Cauchy stress should be symmetric.
    """
    def Cauchy(F,P):
      sigma = np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F['data']),P['data'],F['data'])
      sigma = (sigma + np.einsum('ikj',sigma))*0.5                                     # enforce symmetry
      return  {
               'data'  : sigma,
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
    """Adds the equivalent Mises stress or strain of a tensor."""
    def deviator(x):
      
      if x['meta']['Unit'] == 'Pa': #ToDo: Should we use this? Then add_Cauchy  and add_strain_tensors also should perform sanity checks
        factor = 3.0/2.0
      elif x['meta']['Unit'] == '-':
        factor = 2.0/3.0
      else:
        ValueError
      
      d = x['data']
      dev = d - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[d.shape[0],3,3]),np.trace(d,axis1=1,axis2=2)/3.0)
      #dev_sym = (dev + np.einsum('ikj',dev))*0.5 # ToDo: this is not needed (only if the input is not symmetric, but then the whole concept breaks down)

      return  {
               'data' :  np.sqrt(np.einsum('ijk->i',dev**2)*factor), 
               'label' : 'Mises({})'.format(x['label']),
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

      o = ord
      if   len(x['data'].shape) == 1:
        axis = 0
        t = 'scalar'
        if o is None: o = 2
      elif len(x['data'].shape) == 2:
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
               'label' : '|{}|_{}'.format(x['label'],ord),
               'meta' : {
                        'Unit' :        x['meta']['Unit'],
                        'Description' : '{}-Norm of {} {} ({})'.format(ord,t,x['label'],x['meta']['Description']),
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
      
      if not np.all(np.array(d.shape[1:]) == np.array([3,3])):
        raise ValueError
      
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


  def add_strain_tensor(self,t,ord,defgrad='F'): #ToDo: Use t and ord
    """Adds the a strain tensor."""
    def strain_tensor(defgrad,t,ord):
  #    def operator(stretch,strain,eigenvalues):
  #"""Albrecht Bertram: Elasticity and Plasticity of Large Deformations An Introduction (3rd Edition, 2012), p. 102"""
  #return {
  #  'V#ln':    np.log(eigenvalues)                                 ,
  #  'U#ln':    np.log(eigenvalues)                                 ,
  #  'V#Biot':  ( np.ones(3,'d') - 1.0/eigenvalues )                ,
  #  'U#Biot':  ( eigenvalues - np.ones(3,'d') )                    ,
  #  'V#Green': ( np.ones(3,'d') - 1.0/eigenvalues/eigenvalues) *0.5,
  #  'U#Green': ( eigenvalues*eigenvalues - np.ones(3,'d'))     *0.5,
  #      }[stretch+'#'+strain]
      (U,S,Vh) = np.linalg.svd(defgrad['data'])                                                             # singular value decomposition
      R_inv    = np.einsum('ikj',np.matmul(U,Vh))                                                           # inverse rotation of polar decomposition
      U        = np.matmul(R_inv,defgrad['data'])                                                           # F = RU
      (D,V)    = np.linalg.eigh((U+np.einsum('ikj',U))*.5)                                                  # eigen decomposition (of symmetric(ed) matrix)

      neg      = np.where(D < 0.0)                                                                          # find negative eigenvalues ...
      D[neg[0],neg[1]]   = D[neg[0],neg[1]]* -1                                                             # ... flip value ...
      V[neg[0],:,neg[1]] = V[neg[0],:,neg[1]]* -1                                                           # ... and vector
      
      d = np.log(D)
      a = np.matmul(V,np.einsum('ij,ikj->ijk',d,V))
      
      return  {
               'data' : a, 
               'label' : 'ln(V)({})'.format(defgrad['label']),
               'meta' : {
                        'Unit' :        defgrad['meta']['Unit'],
                        'Description' : 'Strain tensor ln(V){} ({})'.format(defgrad['label'],defgrad['meta']['Description']),
                        'Creator' :     'dadf5.py:add_deviator vXXXXX'
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
      with h5py.File(self.filename,'a') as f:                                                       # write to file
        dataset_out = f[result['group']].create_dataset(result['label'],data=result['data'])
        for k in result['meta'].keys():
          dataset_out.attrs[k] = result['meta'][k]
        N_not_calculated-=1
      
      if N_added < len(todo):                                                                       # add more jobs
        pool.add_task(job,todo[N_added])
        N_added +=1

    pool.wait_completion()
