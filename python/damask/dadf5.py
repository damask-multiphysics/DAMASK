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
               mode     = 'r',
              ):
    """
    Opens an existing DADF5 file.
    
    Parameters
    ----------
    filename : str
        name of the DADF5 file to be openend.
    mode : str, optional
        filemode for opening, either 'r' or 'a'.
    
    """
    if mode not in ['a','r']:
      print('Invalid file access mode')
    else:
      with h5py.File(filename,mode):
        pass
      
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

      self.m_output_types = []
      for m in self.materialpoints:
        for o in f['inc{:05}/materialpoint/{}'.format(self.increments[0]['inc'],m)].keys():
          self.m_output_types.append(o)
      self.m_output_types = list(set(self.m_output_types))                                          # make unique
      
    self.active= {'increments':     self.increments,
                  'constituents':   self.constituents,
                  'materialpoints': self.materialpoints,
                  'constituent':    self.Nconstituents,
                  'c_output_types': self.c_output_types,
                  'm_output_types': self.m_output_types}

    self.filename   = filename
    self.mode       = mode


  def get_candidates(self,l):
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


  def get_active_groups(self):
    """
    Get groups that are currently considered for evaluation. 
    """
    groups = []
    for i,x in enumerate(self.active['increments']):
      group_inc = 'inc{:05}'.format(self.active['increments'][i]['inc'])
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
    

  def list_data(self):
    """Shows information on all active datasets in the file."""
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
    

  def get_dataset_location(self,label):
    """Returns the location of all active datasets with given label."""
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
      return np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F),F,P)
      
    args   = [{'label':F,'shape':[3,3],'unit':'-'},
              {'label':P,'shape':[3,3],'unit':'Pa'} ]
    result = {'label':'sigma',
              'unit':'Pa',
              'Description': 'Cauchy stress calculated from 1st Piola-Kirchhoff stress and deformation gradient'}
    
    self.add_generic_pointwise_vectorized(Cauchy,args,None,result)

    
  def add_Mises_stress(self,stress='sigma'):
    """Adds equivalent von Mises stress."""
    def Mises_stress(stress):
      dev = stress - np.trace(stress)/3.0*np.eye(3)
      symdev = 0.5*(dev+dev.T)
      return np.sqrt(np.sum(symdev*symdev.T)*3.0/2.0)
      
    args   = [{'label':stress,'shape':[3,3],'unit':'Pa'}]
    result = {'label':'Mises({})'.format(stress),
              'unit':'Pa',
              'Description': 'Equivalent Mises stress'}
    
    self.add_generic_pointwise(Mises_stress,args,result)
    
    
  def add_norm(self,x,ord=None):
    """
    Adds norm of vector or tensor or magnitude of a scalar.

    Todo
    ----
      The output unit should be the input unit.
      The ord parameter should be taken into account.
      The whole thing should be vectorized. This requires to parse optional arguments to func.

    """
    args   = [{'label':x,'shape':None,'unit':None}]
    result = {'label':'norm_{}({})'.format(str(ord),x),
              'unit':'n/a',
              'Description': 'Norm of vector or tensor or magnitude of a scalar. See numpy.linalg.norm manual for details'}
    
    self.add_generic_pointwise(np.linalg.norm,args,result)
  
  
  def add_determinant(self,a):
    """Adds the determinant of a tensor."""
    # ToDo: The output unit should be the input unit
    args   = [{'label':a,'shape':[3,3],'unit':None}]
    result = {'label':'det({})'.format(a),
              'unit':'n/a',
              'Description': 'Determinant of a tensor'}
    
    self.add_generic_pointwise_vectorized(np.linalg.det,args,result)


  def add_spherical(self,a):
    """Adds the spherical component of a tensor."""
    def spherical(m):
      return (m[0,0]+m[1,1]+m[2,2])/3.0

    # ToDo: The output unit should be the input unit  
    args   = [{'label':a,'shape':[3,3],'unit':None}]
    result = {'label':'sph({})'.format(a),
              'unit':'n/a',
              'Description': 'Spherical component of a tensor'}

    self.add_generic_pointwise(spherical,args,result)
  
  
  def add_deviator(self,a):
    """Adds the deviator of a tensor."""
    def deviator(m):
      return m - np.eye(3)*(m[0,0]+m[1,1]+m[2,2])/3.0

    # ToDo: The output unit should be the input unit  
    args   = [{'label':a,'shape':[3,3],'unit':'Pa'}]
    result = {'label':'dev({})'.format(a),
              'unit':'n/a',
              'Description': 'Deviatoric component of a tensor'}

    self.add_generic_pointwise(deviator,args,result)

    
    
  def add_strain_tensors(self,defgrad='F'):
    """Adds a strain definition."""
    def strain(defgrad):
      (U,S,Vh) = np.linalg.svd(defgrad)                                                             # singular value decomposition
      R_inv    = np.dot(U,Vh).T                                                                     # inverse rotation of polar decomposition
      U        = np.dot(R_inv,defgrad)                                                              # F = RU
      U        = np.where(abs(U) < 1e-12, 0, U)                                                     # kill nasty noisy data
      (D,V) = np.linalg.eig(U)                                                                      # eigen decomposition (of symmetric matrix)
      neg   = np.where(D < 0.0)                                                                     # find negative eigenvalues ...
      D[neg]   *= -1.                                                                               # ... flip value ...
      V[:,neg] *= -1.                                                                               # ... and vector
      for i,eigval in enumerate(D):
        if np.dot(V[:,i],V[:,(i+1)%3]) != 0.0:                                                      # check each vector for orthogonality
          V[:,(i+1)%3] = np.cross(V[:,(i+2)%3],V[:,i])                                              # correct next vector
          V[:,(i+1)%3] /= np.sqrt(np.dot(V[:,(i+1)%3],V[:,(i+1)%3].conj()))                         # and renormalize (hyperphobic?)
      d = np.log(D)                                                                                 # operate on eigenvalues of U o r V
      return np.dot(V,np.dot(np.diag(d),V.T)).real                                                  # build tensor back from eigenvalue/vector basis

    # ToDo: The output unit should be the input unit
    args   = [{'label':defgrad,'shape':[3,3],'unit':None}]
    result = {'label':'strain({})'.format(defgrad),
              'unit':'-',
              'Description': 'strain (ln(V)) of a deformation gradient'}
    
    self.add_generic_pointwise(strain,args,result)
    

  def get_fitting(self,data):
    groups = []
    if type(data) is not list:
      print('mist')
    with h5py.File(self.filename,'r') as f:
      for g in self.get_candidates([l['label'] for l in data]):
        print(g)
        fits = True
        for d in data:    # ToDo: check for unit
          if d['shape'] is not None:
            fits = fits and np.all(np.array(f[g+'/'+d['label']].shape[1:]) == np.array(d['shape']))
        if fits: groups.append(g)
    return groups

    
  def add_generic_pointwise(self,func,args,result):
    """
    General function to add pointwise data.
    
    function 'func' first needs to have data arguments before other arguments
    Works for functions that are pointwise defined.
    """
    groups = self.get_fitting(args)
    
    def job(args):
      out         = args['out']
      datasets_in = args['dat']
      func        = args['fun']
      for i in range(out.shape[0]):
        arg = tuple([d[i,] for d in datasets_in]) 
        out[i,] = func(*arg)
      args['results'].put({'out':out,'group':args['group']})
    
    Nthreads = 4 # ToDo: should be a parameter
    results  = Queue(Nthreads+1)
    
    todo = []
    
    for g in groups:
      with h5py.File(self.filename,'r') as f:
        datasets_in = [f[g+'/'+u['label']][()] for u in args]
        
      # figure out dimension of results
      testArg = tuple([d[0,] for d in datasets_in])                                                 # to call function with first point
      out = np.empty([datasets_in[0].shape[0]] + list(func(*testArg).shape))                        # shape is Npoints x shape of the results for one point
      todo.append({'dat':datasets_in,'fun':func,'out':out,'group':g,'results':results})
    
    # Instantiate a thread pool with worker threads
    pool = util.ThreadPool(Nthreads)
    missingResults = len(todo)


    # Add the jobs in bulk to the thread pool. Alternatively you could use
    # `pool.add_task` to add single jobs. The code will block here, which
    # makes it possible to cancel the thread pool with an exception when
    # the currently running batch of workers is finished

    pool.map(job, todo[:Nthreads+1])
    i = 0
    while missingResults > 0:    
        r=results.get() # noqa
        print(r['group'])
        with h5py.File(self.filename,'r+') as f:
           dataset_out = f[r['group']].create_dataset(result['label'],data=r['out'])
           dataset_out.attrs['Unit'] = result['unit']
           dataset_out.attrs['Description'] = result['Description']
           dataset_out.attrs['Creator'] = 'dadf5.py v{}'.format('n/a')
        missingResults-=1
        try:
          pool.add_task(job,todo[Nthreads+1+i])
        except IndexError:
          pass
        i+=1   

    pool.wait_completion()
    
    
  def add_generic_pointwise_vectorized(self,func,args,args2=None,result=None):
    """
    General function to add pointwise data.
    
    function 'func' first needs to have data arguments before other arguments
    Works for vectorized functions.
    """
    groups = self.get_fitting(args)
    
    def job(args):
      """
      A job. It has different args!
      """
      print('args for job',args)
      out         = args['out']
      datasets_in = args['dat']
      func        = args['fun']
#      try:
 #       out         = func(*datasets_in,*args['fun_args'])
  #    except:
      out         = func(*datasets_in)
      args['results'].put({'out':out,'group':args['group']})
    
    Nthreads = 4 # ToDo: should be a parameter
    results  = Queue(Nthreads+1)
    
    todo = []
    
    for g in groups:
      with h5py.File(self.filename,'r') as f:
        datasets_in = [f[g+'/'+u['label']][()] for u in args]

      if args2 is not None:
        todo.append({'dat':datasets_in,'fun':func,'group':g,'results':results,'func_args':args,'out':None})
      else:  
        todo.append({'dat':datasets_in,'fun':func,'group':g,'results':results,'out':None})
    
    # Instantiate a thread pool with worker threads
    pool = util.ThreadPool(Nthreads)
    missingResults = len(todo)


    pool.map(job, todo[:Nthreads+1])
    i = 0
    while missingResults > 0:    
        r=results.get() # noqa
        print(r['group'])
        with h5py.File(self.filename,'r+') as f:
           dataset_out = f[r['group']].create_dataset(result['label'],data=r['out'])
           dataset_out.attrs['Unit'] = result['unit']
           dataset_out.attrs['Description'] = result['Description']
           dataset_out.attrs['Creator'] = 'dadf5.py v{}'.format('n/a')
        missingResults-=1
        try:
          pool.add_task(job,todo[Nthreads+1+i])
        except IndexError:
          pass
        i+=1   

    pool.wait_completion()
