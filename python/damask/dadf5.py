from queue import Queue
import re
import glob

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
      self.time_information = [{'inc':  int(u[3:]),
                          'time': round(f[u].attrs['time/s'],12),
                          }   for u in f.keys() if r.match(u)]
                          
      self.increments = self.time_information.copy() # unify later
      
      self.Nmaterialpoints, self.Nconstituents =   np.shape(f['mapping/cellResults/constituent'])
      self.materialpoints  = [m.decode() for m in np.unique(f['mapping/cellResults/materialpoint']['Name'])]
      self.constituents    = [c.decode() for c in np.unique(f['mapping/cellResults/constituent']  ['Name'])]


      self.con_physics  = []
      for c in self.constituents:
        self.con_physics += f['inc{:05}/constituent/{}'.format(self.increments[0]['inc'],c)].keys()
      self.con_physics = list(set(self.con_physics))                                          # make unique

      self.mat_physics = []
      for m in self.materialpoints:
        self.mat_physics += f['inc{:05}/materialpoint/{}'.format(self.increments[0]['inc'],m)].keys()
      self.mat_physics = list(set(self.mat_physics))                                          # make unique

    self.visible= {'increments':    self.increments, # ToDo:simplify, activity only positions that translate into (no complex types)
                  'constituents':   self.constituents,
                  'materialpoints': self.materialpoints,
                  'constituent':    range(self.Nconstituents),                                      # ToDo: stupid naming
                  'con_physics': self.con_physics,
                  'mat_physics': self.mat_physics}
                  
    self.filename   = filename

  def __visible_set(self,output,t):
    """Sets visible."""
    # allow True/False and string arguments
    if output is True:
      output = ['*']
    elif output is False:
      output = []
    choice = [output] if isinstance(output,str) else output
    
    valid = [e for e_ in [glob.fnmatch.filter(getattr(self,t),s) for s in choice] for e in e_]
    
    self.visible[t] = valid


  def __visible_add(self,output,t):
    """Adds to visible."""
    # allow True/False and string arguments
    if output is True:
      output = ['*']
    elif output is False:
      output = []
    choice = [output] if isinstance(output,str) else output
    
    existing = set(self.visible[t])
    valid    = [e for e_ in [glob.fnmatch.filter(getattr(self,t),s) for s in choice] for e in e_]

    self.visible[t] = list(existing.union(valid))


  def __visible_del(self,output,t):
    """Deletes from visible."""
    # allow True/False and string arguments
    if output is True:
      output = ['*']
    elif output is False:
      output = []
    choice = [output] if isinstance(output,str) else output
    
    existing = set(self.visible[t])
    valid    = [e for e_ in [glob.fnmatch.filter(existing,s) for s in choice] for e in e_]

    self.visible[t] = list(existing.difference_update(valid))


  def __visible_iter(self,t):
    a = self.visible[t]
    last_a = a.copy()
    for i in a:
      if last_a != self.visible[t]:
        self.__visible_set(a,t)
        raise Exception
      self.__visible_set(i,t)
      last_a = self.visible[t]
      yield i
    self.__visible_set(a,t)


# ToDo: store increments, select icrements (trivial), position, and time
  def increment_set_by_time(self,start,end):
    for t in self.time_information:
      if start<= t['time']< end:
        print(t)


  def increment_set_by_position(self,start,end):
    for t in self.time_information[start:end]:
      print(t)


  def increment_set(self,start,end):
    for t in self.time_information:
      if start<= t['inc']< end:
        print(t)
      
      
  def constituent_output_iter(self):
    return self.__visible_iter('con_physics')
  
    
  def constituent_output_set(self,output):
    self.__visible_set(output,'con_physics')


  def constituent_output_add(self,output):
    self.__visible_add(output,'con_physics')


  def constituent_output_del(self,output):
    self.__visible_del(output,'con_physics')
  
    
  def materialpoint_output_iter(self):
    return self.__visible_iter('mat_physics')
  
    
  def materialpoint_output_set(self,output):
    self.__visible_set(output,'mat_physics')


  def materialpoint_output_add(self,output):
    self.__visible_add(output,'mat_physics')


  def materialpoint_output_del(self,output):
    self.__visible_del(output,'mat_physics')
  
  
  def constituent_iter(self):
    return self.__visible_iter('constituents')
    
    
  def constituent_set(self,output):
    self.__visible_set(output,'constituents')


  def constituent_add(self,output):
    self.__visible_add(output,'constituents')


  def constituent_del(self,output):
    self.__visible_del(output,'constituents')
  
  
  def materialpoint_iter(self):
    return self.__visible_iter('materialpoints')
    

  def materialpoint_set(self,output):
    self.__visible_set(output,'materialpoints')


  def materialpoint_add(self,output):
    self.__visible_add(output,'materialpoints')


  def materialpoint_del(self,output):
    self.__visible_del(output,'materialpoints')


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
      for i in self.visible['increments']:
        group_inc = 'inc{:05}'.format(i['inc'])                               #ToDo: Merge path only once at the end '/'.join(listE)
        for c in self.constituent_iter():
          for t in self.constituent_output_iter():
            group = '/'.join([group_inc,'constituent',c,t])
            if sets is True:
              groups.append(group)
            else:
              match = [e for e_ in [glob.fnmatch.filter(f[group].keys(),s) for s in sets] for e in e_]
              if len(set(match)) == len(sets) : groups.append(group)
        for m in self.materialpoint_iter():
          for t in self.materialpoint_output_iter():
            group = '/'.join([group_inc,'materialpoint',m,t])
            if sets is True:
              groups.append(group)
            else:
              match = [e for e_ in [glob.fnmatch.filter(f[group].keys(),s) for s in sets] for e in e_]
              if len(set(match)) == len(sets) : groups.append(group)
    return groups


  def list_data(self): # print_datasets and have [] and ['*'], loop over all increment, soll auf anderen basieren (get groups with sternchen)
    """Shows information on all active datasets in the file."""
    with h5py.File(self.filename,'r') as f:
      group_inc = 'inc{:05}'.format(self.visible['increments'][0]['inc']) #ToDo: Merge path only once at the end '/'.join(listE)
      for c in self.visible['constituents']:
        print('\n'+c)
        group_constituent = group_inc+'/constituent/'+c
        for t in self.visible['con_physics']:
          print('  {}'.format(t))
          group_output_types = group_constituent+'/'+t 
          try:
            for x in f[group_output_types].keys():
              print('    {} ({})'.format(x,f[group_output_types+'/'+x].attrs['Description'].decode()))
          except KeyError:
            pass
      for m in self.visible['materialpoints']:
        group_materialpoint = group_inc+'/materialpoint/'+m
        for t in self.visible['mat_physics']:
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
          for t in self.visible['con_physics']:
            try:
              p = '/'.join([group_inc,'constituent',c,t,label])
              f[p]
              path.append(p)
            except KeyError as e:
              print('unable to locate constituents dataset: '+ str(e))
       
        for m in self.visible['materialpoints']:
          for t in self.visible['mat_physics']:
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
                                        'and deformation gradient {} ({})'.format(F['label'],F['meta']['Description']),
                        'Creator' :     'dadf5.py:add_Cauchy vXXXXX'
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
                        'Creator' :     'dadf5.py:add_Mises_stress vXXXXX'
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
                        'Creator' :     'dadf5.py:add_norm vXXXXX'
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
                        'Creator' :     'dadf5.py:add_abs vXXXXX'
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
                        'Creator' :     'dadf5.py:add_calculation vXXXXX'
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
                   'U#Biot': lambda U: U**2 - np.broadcast_to(np.ones(3),[U.shape[0],3]), 
                 }

      (U,S,Vh) = np.linalg.svd(defgrad['data'])                                                             # singular value decomposition
      R_inv    = np.einsum('ikj',np.matmul(U,Vh))                                                           # inverse rotation of polar decomposition
      U        = np.matmul(R_inv,defgrad['data'])                                                           # F = RU
      (D,V)    = np.linalg.eigh((U+np.einsum('ikj',U))*.5)                                                  # eigen decomposition (of symmetric(ed) matrix)

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
