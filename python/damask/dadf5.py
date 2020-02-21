import multiprocessing
import re
import glob
import os
from functools import partial

import vtk
from vtk.util import numpy_support
import h5py
import numpy as np

from . import util
from . import version
from . import mechanics
from . import Rotation
from . import Orientation
from . import Environment
from . import grid_filters

class DADF5():
    """
    Read and write to DADF5 files.

    DADF5 files contain DAMASK results.
    """

    def __init__(self,fname):
        """
        Opens an existing DADF5 file.

        Parameters
        ----------
        fname : str
            name of the DADF5 file to be openend.

        """
        with h5py.File(fname,'r') as f:

          try:
              self.version_major = f.attrs['DADF5_version_major']
              self.version_minor = f.attrs['DADF5_version_minor']
          except KeyError:
              self.version_major = f.attrs['DADF5-major']
              self.version_minor = f.attrs['DADF5-minor']

          if self.version_major != 0 or not 2 <= self.version_minor <= 6:
              raise TypeError('Unsupported DADF5 version {}.{} '.format(self.version_major,
                                                                        self.version_minor))

          self.structured = 'grid' in f['geometry'].attrs.keys()

          if self.structured:
              self.grid   = f['geometry'].attrs['grid']
              self.size   = f['geometry'].attrs['size']
              self.origin = f['geometry'].attrs['origin'] if self.version_major == 0 and self.version_minor >= 5 else \
                            np.zeros(3)

          r=re.compile('inc[0-9]+')
          increments_unsorted = {int(i[3:]):i for i in f.keys() if r.match(i)}
          self.increments     = [increments_unsorted[i] for i in sorted(increments_unsorted)]
          self.times          = [round(f[i].attrs['time/s'],12) for i in self.increments]

          self.Nmaterialpoints, self.Nconstituents =   np.shape(f['mapping/cellResults/constituent'])
          self.materialpoints  = [m.decode() for m in np.unique(f['mapping/cellResults/materialpoint']['Name'])]
          self.constituents    = [c.decode() for c in np.unique(f['mapping/cellResults/constituent']  ['Name'])]

          self.con_physics  = []
          for c in self.constituents:
              self.con_physics += f['/'.join([self.increments[0],'constituent',c])].keys()
          self.con_physics = list(set(self.con_physics))                                                # make unique

          self.mat_physics = []
          for m in self.materialpoints:
              self.mat_physics += f['/'.join([self.increments[0],'materialpoint',m])].keys()
          self.mat_physics = list(set(self.mat_physics))                                                # make unique

        self.selection= {'increments':     self.increments,
                         'constituents':   self.constituents,
                         'materialpoints': self.materialpoints,
                         'con_physics':    self.con_physics,
                         'mat_physics':    self.mat_physics}

        self.fname = fname


    def __manage_visible(self,datasets,what,action):
        """
        Manages the visibility of the groups.

        Parameters
        ----------
        datasets : list of str or Boolean
          name of datasets as list, supports ? and * wildcards.
          True is equivalent to [*], False is equivalent to []
        what : str
          attribute to change (must be in self.selection)
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
        existing = set(self.selection[what])

        if   action == 'set':
          self.selection[what] = valid
        elif action == 'add':
          self.selection[what] = list(existing.union(valid))
        elif action == 'del':
          self.selection[what] = list(existing.difference_update(valid))


    def __time_to_inc(self,start,end):
        selected = []
        for i,time in enumerate(self.times):
          if start <= time <= end:
            selected.append(self.increments[i])
        return selected


    def set_by_time(self,start,end):
        """
        Set active increments based on start and end time.

        Parameters
        ----------
        start : float
          start time (included)
        end : float
          end time (included)

        """
        self.__manage_visible(self.__time_to_inc(start,end),'increments','set')


    def add_by_time(self,start,end):
        """
        Add to active increments based on start and end time.

        Parameters
        ----------
        start : float
          start time (included)
        end : float
          end time (included)

        """
        self.__manage_visible(self.__time_to_inc(start,end),'increments','add')


    def del_by_time(self,start,end):
        """
        Delete from active increments based on start and end time.

        Parameters
        ----------
        start : float
          start time (included)
        end : float
          end time (included)

        """
        self.__manage_visible(self.__time_to_inc(start,end),'increments','del')


    def set_by_increment(self,start,end):
        """
        Set active time increments based on start and end increment.

        Parameters
        ----------
        start : int
          start increment (included)
        end : int
          end increment (included)

        """
        if self.version_minor >= 4:
          self.__manage_visible([    'inc{}'.format(i) for i in range(start,end+1)],'increments','set')
        else:
          self.__manage_visible(['inc{:05d}'.format(i) for i in range(start,end+1)],'increments','set')


    def add_by_increment(self,start,end):
        """
        Add to active time increments based on start and end increment.

        Parameters
        ----------
        start : int
          start increment (included)
        end : int
          end increment (included)

        """
        if self.version_minor >= 4:
          self.__manage_visible([    'inc{}'.format(i) for i in range(start,end+1)],'increments','add')
        else:
          self.__manage_visible(['inc{:05d}'.format(i) for i in range(start,end+1)],'increments','add')


    def del_by_increment(self,start,end):
        """
        Delet from active time increments based on start and end increment.

        Parameters
        ----------
        start : int
          start increment (included)
        end : int
          end increment (included)

        """
        if self.version_minor >= 4:
          self.__manage_visible([    'inc{}'.format(i) for i in range(start,end+1)],'increments','del')
        else:
          self.__manage_visible(['inc{:05d}'.format(i) for i in range(start,end+1)],'increments','del')


    def iter_visible(self,what):
        """
        Iterate over visible items by setting each one visible.

        Parameters
        ----------
        what : str
          attribute to change (must be in self.selection)

        """
        datasets = self.selection[what]
        last_datasets = datasets.copy()
        for dataset in datasets:
          if last_datasets != self.selection[what]:
            self.__manage_visible(datasets,what,'set')
            raise Exception
          self.__manage_visible(dataset,what,'set')
          last_datasets = self.selection[what]
          yield dataset
        self.__manage_visible(datasets,what,'set')


    def set_visible(self,what,datasets):
        """
        Set active groups.

        Parameters
        ----------
        datasets : list of str or Boolean
          name of datasets as list, supports ? and * wildcards.
          True is equivalent to [*], False is equivalent to []
        what : str
          attribute to change (must be in self.selection)

        """
        self.__manage_visible(datasets,what,'set')


    def add_visible(self,what,datasets):
        """
        Add to active groups.

        Parameters
        ----------
        datasets : list of str or Boolean
          name of datasets as list, supports ? and * wildcards.
          True is equivalent to [*], False is equivalent to []
        what : str
          attribute to change (must be in self.selection)

        """
        self.__manage_visible(datasets,what,'add')


    def del_visible(self,what,datasets):
        """
        Delete from active groupe.

        Parameters
        ----------
        datasets : list of str or Boolean
          name of datasets as list, supports ? and * wildcards.
          True is equivalent to [*], False is equivalent to []
        what : str
          attribute to change (must be in self.selection)

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

        with h5py.File(self.fname,'r') as f:
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
        """Return information on all active datasets in the file."""
        message = ''
        with h5py.File(self.fname,'r') as f:
          for i in self.iter_visible('increments'):
            message+='\n{} ({}s)\n'.format(i,self.times[self.increments.index(i)])
            for o,p in zip(['constituents','materialpoints'],['con_physics','mat_physics']):
              for oo in self.iter_visible(o):
                message+='  {}\n'.format(oo)
                for pp in self.iter_visible(p):
                  message+='    {}\n'.format(pp)
                  group = '/'.join([i,o[:-1],oo,pp])                              # o[:-1]: plural/singular issue
                  for d in f[group].keys():
                    try:
                      dataset = f['/'.join([group,d])]
                      message+='      {} / ({}): {}\n'.\
                                format(d,dataset.attrs['Unit'].decode(),dataset.attrs['Description'].decode())
                    except KeyError:
                      pass
        return message


    def get_dataset_location(self,label):
        """Return the location of all active datasets with given label."""
        path = []
        with h5py.File(self.fname,'r') as f:
          for i in self.iter_visible('increments'):
            k = '/'.join([i,'geometry',label])
            try:
              f[k]
              path.append(k)
            except KeyError as e:
              pass
            for o,p in zip(['constituents','materialpoints'],['con_physics','mat_physics']):
              for oo in self.iter_visible(o):
                for pp in self.iter_visible(p):
                  k = '/'.join([i,o[:-1],oo,pp,label])
                  try:
                    f[k]
                    path.append(k)
                  except KeyError as e:
                    pass
        return path


    def get_constituent_ID(self,c=0):
        """Pointwise constituent ID."""
        with h5py.File(self.fname,'r') as f:
          names = f['/mapping/cellResults/constituent']['Name'][:,c].astype('str')
        return np.array([int(n.split('_')[0]) for n in names.tolist()],dtype=np.int32)


    def get_crystal_structure(self):                                                                # ToDo: extension to multi constituents/phase
        """Info about the crystal structure."""
        with h5py.File(self.fname,'r') as f:
          return f[self.get_dataset_location('orientation')[0]].attrs['Lattice'].astype('str')      # np.bytes_ to string


    def read_dataset(self,path,c=0,plain=False):
        """
        Dataset for all points/cells.

        If more than one path is given, the dataset is composed of the individual contributions.
        """
        with h5py.File(self.fname,'r') as f:
          shape = (self.Nmaterialpoints,) + np.shape(f[path[0]])[1:]
          if len(shape) == 1: shape = shape +(1,)
          dataset = np.full(shape,np.nan,dtype=np.dtype(f[path[0]]))
          for pa in path:
            label = pa.split('/')[2]

            if (pa.split('/')[1] == 'geometry'):
              dataset = np.array(f[pa])
              continue

            p = np.where(f['mapping/cellResults/constituent'][:,c]['Name'] == str.encode(label))[0]
            if len(p)>0:
              u = (f['mapping/cellResults/constituent']['Position'][p,c])
              a = np.array(f[pa])
              if len(a.shape) == 1:
                a=a.reshape([a.shape[0],1])
              dataset[p,:] = a[u,:]

            p = np.where(f['mapping/cellResults/materialpoint']['Name'] == str.encode(label))[0]
            if len(p)>0:
              u = (f['mapping/cellResults/materialpoint']['Position'][p.tolist()])
              a = np.array(f[pa])
              if len(a.shape) == 1:
                a=a.reshape([a.shape[0],1])
              dataset[p,:] = a[u,:]

        if plain and dataset.dtype.names is not None:
          return dataset.view(('float64',len(dataset.dtype.names)))
        else:
          return dataset


    def cell_coordinates(self):
        """Return initial coordinates of the cell centers."""
        if self.structured:
            return grid_filters.cell_coord0(self.grid,self.size,self.origin)
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_c'][()]


    @staticmethod
    def _add_absolute(x):
        return {
                'data':  np.abs(x['data']),
                'label': '|{}|'.format(x['label']),
                'meta':  {
                          'Unit':        x['meta']['Unit'],
                          'Description': 'Absolute value of {} ({})'.format(x['label'],x['meta']['Description']),
                          'Creator':     'dadf5.py:add_abs v{}'.format(version)
                          }
                 }
    def add_absolute(self,x):
        """
        Add absolute value.

        Parameters
        ----------
        x : str
          Label of scalar, vector, or tensor dataset to take absolute value of.

        """
        self.__add_generic_pointwise(self._add_absolute,{'x':x})


    @staticmethod
    def _add_calculation(**kwargs):
        formula = kwargs['formula']
        for d in re.findall(r'#(.*?)#',formula):
            formula = formula.replace('#{}#'.format(d),"kwargs['{}']['data']".format(d))

        return {
                'data':  eval(formula),
                'label': kwargs['label'],
                'meta':  {
                          'Unit':        kwargs['unit'],
                          'Description': '{} (formula: {})'.format(kwargs['description'],kwargs['formula']),
                          'Creator':     'dadf5.py:add_calculation v{}'.format(version)
                          }
                 }
    def add_calculation(self,label,formula,unit='n/a',description=None,vectorized=True):
        """
        Add result of a general formula.

        Parameters
        ----------
        label : str
          Label of resulting dataset.
        formula : str
          Formula to calculate resulting dataset. Existing datasets are referenced by ‘#TheirLabel#‘.
        unit : str, optional
          Physical unit of the result.
        description : str, optional
          Human-readable description of the result.
        vectorized : bool, optional
          Indicate whether the formula can be used in vectorized form. Defaults to ‘True’.

        """
        if not vectorized:
          raise NotImplementedError

        dataset_mapping  = {d:d for d in set(re.findall(r'#(.*?)#',formula))}                       # datasets used in the formula
        args             = {'formula':formula,'label':label,'unit':unit,'description':description}
        self.__add_generic_pointwise(self._add_calculation,dataset_mapping,args)


    @staticmethod
    def _add_Cauchy(P,F):
        return {
                'data':  mechanics.Cauchy(P['data'],F['data']),
                'label': 'sigma',
                'meta':  {
                          'Unit':        P['meta']['Unit'],
                          'Description': 'Cauchy stress calculated from {} ({}) '.format(P['label'],
                                                                                         P['meta']['Description'])+\
                                         'and {} ({})'.format(F['label'],F['meta']['Description']),
                          'Creator':     'dadf5.py:add_Cauchy v{}'.format(version)
                          }
                }
    def add_Cauchy(self,P='P',F='F'):
        """
        Add Cauchy stress calculated from first Piola-Kirchhoff stress and deformation gradient.

        Parameters
        ----------
        P : str, optional
          Label of the dataset containing the first Piola-Kirchhoff stress. Defaults to ‘P’.
        F : str, optional
          Label of the dataset containing the deformation gradient. Defaults to ‘F’.

        """
        self.__add_generic_pointwise(self._add_Cauchy,{'P':P,'F':F})


    @staticmethod
    def _add_determinant(T):
        return {
                'data':  np.linalg.det(T['data']),
                'label': 'det({})'.format(T['label']),
                'meta':  {
                          'Unit':        T['meta']['Unit'],
                          'Description': 'Determinant of tensor {} ({})'.format(T['label'],T['meta']['Description']),
                          'Creator':     'dadf5.py:add_determinant v{}'.format(version)
                          }
                }
    def add_determinant(self,T):
        """
        Add the determinant of a tensor.

        Parameters
        ----------
        T : str
          Label of tensor dataset.

        """
        self.__add_generic_pointwise(self._add_determinant,{'T':T})


    @staticmethod
    def _add_deviator(T):
        if not T['data'].shape[1:] == (3,3):
            raise ValueError

        return {
                'data':  mechanics.deviatoric_part(T['data']),
                'label': 's_{}'.format(T['label']),
                'meta':  {
                          'Unit':        T['meta']['Unit'],
                          'Description': 'Deviator of tensor {} ({})'.format(T['label'],T['meta']['Description']),
                          'Creator':     'dadf5.py:add_deviator v{}'.format(version)
                          }
                 }
    def add_deviator(self,T):
        """
        Add the deviatoric part of a tensor.

        Parameters
        ----------
        T : str
          Label of tensor dataset.

        """
        self.__add_generic_pointwise(self._add_deviator,{'T':T})


    @staticmethod
    def _add_eigenvalue(S):
        return {
                'data': mechanics.eigenvalues(S['data']),
                'label': 'lambda({})'.format(S['label']),
                'meta' : {
                          'Unit':         S['meta']['Unit'],
                          'Description': 'Eigenvalues of {} ({})'.format(S['label'],S['meta']['Description']),
                          'Creator':     'dadf5.py:add_eigenvalues v{}'.format(version)
                         }
                }
    def add_eigenvalues(self,S):
        """
        Add eigenvalues of symmetric tensor.

        Parameters
        ----------
        S : str
          Label of symmetric tensor dataset.

        """
        self.__add_generic_pointwise(self._add_eigenvalue,{'S':S})


    @staticmethod
    def _add_eigenvector(S):
        return {
                'data': mechanics.eigenvectors(S['data']),
                'label': 'v({})'.format(S['label']),
                'meta' : {
                          'Unit':        '1',
                          'Description': 'Eigenvectors of {} ({})'.format(S['label'],S['meta']['Description']),
                          'Creator':     'dadf5.py:add_eigenvectors v{}'.format(version)
                         }
                }
    def add_eigenvectors(self,S):
        """
        Add eigenvectors of symmetric tensor.

        Parameters
        ----------
        S : str
          Label of symmetric tensor dataset.

        """
        self.__add_generic_pointwise(self._add_eigenvector,{'S':S})


    @staticmethod
    def _add_IPFcolor(q,l):
        d      = np.array(l)
        d_unit = d/np.linalg.norm(d)
        m      = util.scale_to_coprime(d)
        colors = np.empty((len(q['data']),3),np.uint8)

        lattice   = q['meta']['Lattice']

        for i,q in enumerate(q['data']):
             o = Orientation(np.array([q['w'],q['x'],q['y'],q['z']]),lattice).reduced()
             colors[i] = np.uint8(o.IPFcolor(d_unit)*255)

        return {
                'data': colors,
                'label': 'IPFcolor_[{} {} {}]'.format(*m),
                'meta' : {
                          'Unit':        'RGB (8bit)',
                          'Lattice':     lattice,
                          'Description': 'Inverse Pole Figure (IPF) colors for direction/plane [{} {} {})'.format(*m),
                          'Creator':     'dadf5.py:add_IPFcolor v{}'.format(version)
                         }
               }
    def add_IPFcolor(self,q,l):
        """
        Add RGB color tuple of inverse pole figure (IPF) color.

        Parameters
        ----------
        q : str
          Label of the dataset containing the crystallographic orientation as quaternions.
        l : numpy.array of shape (3)
          Lab frame direction for inverse pole figure.

        """
        self.__add_generic_pointwise(self._add_IPFcolor,{'q':q},{'l':l})


    @staticmethod
    def _add_maximum_shear(S):
        return {
                'data':  mechanics.maximum_shear(S['data']),
                'label': 'max_shear({})'.format(S['label']),
                'meta':  {
                          'Unit':        S['meta']['Unit'],
                          'Description': 'Maximum shear component of {} ({})'.format(S['label'],S['meta']['Description']),
                          'Creator':     'dadf5.py:add_maximum_shear v{}'.format(version)
                          }
                 }
    def add_maximum_shear(self,S):
        """
        Add maximum shear components of symmetric tensor.

        Parameters
        ----------
        S : str
          Label of symmetric tensor dataset.

        """
        self.__add_generic_pointwise(self._add_maximum_shear,{'S':S})


    @staticmethod
    def _add_Mises(S):
        t = 'strain' if S['meta']['Unit'] == '1' else \
            'stress'

        return {
                'data':  mechanics.Mises_strain(S['data']) if t=='strain' else mechanics.Mises_stress(S['data']),
                'label': '{}_vM'.format(S['label']),
                'meta':  {
                          'Unit':        S['meta']['Unit'],
                          'Description': 'Mises equivalent {} of {} ({})'.format(t,S['label'],S['meta']['Description']),
                          'Creator':     'dadf5.py:add_Mises v{}'.format(version)
                          }
                }
    def add_Mises(self,S):
        """
        Add the equivalent Mises stress or strain of a symmetric tensor.

        Parameters
        ----------
        S : str
          Label of symmetric tensorial stress or strain dataset.

        """
        self.__add_generic_pointwise(self._add_Mises,{'S':S})


    @staticmethod
    def _add_norm(x,ord):
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

        return {
                'data':  np.linalg.norm(x['data'],ord=o,axis=axis,keepdims=True),
                'label': '|{}|_{}'.format(x['label'],o),
                'meta':  {
                          'Unit':        x['meta']['Unit'],
                          'Description': '{}-norm of {} {} ({})'.format(o,t,x['label'],x['meta']['Description']),
                          'Creator':     'dadf5.py:add_norm v{}'.format(version)
                          }
                 }
    def add_norm(self,x,ord=None):
        """
        Add the norm of vector or tensor.

        Parameters
        ----------
        x : str
          Label of vector or tensor dataset.
        ord : {non-zero int, inf, -inf, ‘fro’, ‘nuc’}, optional
          Order of the norm. inf means NumPy’s inf object. For details refer to numpy.linalg.norm.

        """
        self.__add_generic_pointwise(self._add_norm,{'x':x},{'ord':ord})


    @staticmethod
    def _add_PK2(P,F):
        return {
                'data':  mechanics.PK2(P['data'],F['data']),
                'label': 'S',
                'meta':  {
                          'Unit':        P['meta']['Unit'],
                          'Description': '2. Kirchhoff stress calculated from {} ({}) '.format(P['label'],
                                                                                               P['meta']['Description'])+\
                                         'and {} ({})'.format(F['label'],F['meta']['Description']),
                          'Creator':     'dadf5.py:add_PK2 v{}'.format(version)
                          }
                }
    def add_PK2(self,P='P',F='F'):
        """
        Add second Piola-Kirchhoff calculated from first Piola-Kirchhoff stress and deformation gradient.

        Parameters
        ----------
        P : str, optional
          Label first  Piola-Kirchhoff stress dataset. Defaults to ‘P’.
        F : str, optional
          Label of deformation gradient dataset. Defaults to ‘F’.

        """
        self.__add_generic_pointwise(self._add_PK2,{'P':P,'F':F})


    @staticmethod
    def _add_pole(q,p,polar):
        pole      = np.array(p)
        unit_pole = pole/np.linalg.norm(pole)
        m         = util.scale_to_coprime(pole)
        coords    = np.empty((len(q['data']),2))

        for i,q in enumerate(q['data']):
            o = Rotation(np.array([q['w'],q['x'],q['y'],q['z']]))
            rotatedPole = o*unit_pole                                                               # rotate pole according to crystal orientation
            (x,y) = rotatedPole[0:2]/(1.+abs(unit_pole[2]))                                         # stereographic projection
            coords[i] = [np.sqrt(x*x+y*y),np.arctan2(y,x)] if polar else [x,y]

        return {
                'data': coords,
                'label': 'p^{}_[{} {} {})'.format(u'rφ' if polar else 'xy',*m),
                'meta' : {
                          'Unit': '1',
                          'Description': '{} coordinates of stereographic projection of pole (direction/plane) in crystal frame'\
                                         .format('Polar' if polar else 'Cartesian'),
                          'Creator' : 'dadf5.py:add_pole v{}'.format(version)
                         }
               }
    def add_pole(self,q,p,polar=False):
        """
        Add coordinates of stereographic projection of given pole in crystal frame.

        Parameters
        ----------
        q : str
          Label of the dataset containing the crystallographic orientation as quaternions.
        p : numpy.array of shape (3)
          Crystallographic direction or plane.
        polar : bool, optional
          Give pole in polar coordinates. Defaults to False.

        """
        self.__add_generic_pointwise(self._add_pole,{'q':q},{'p':p,'polar':polar})


    @staticmethod
    def _add_rotational_part(F):
        if not F['data'].shape[1:] == (3,3):
            raise ValueError
        return {
                'data':  mechanics.rotational_part(F['data']),
                'label': 'R({})'.format(F['label']),
                'meta':  {
                          'Unit':        F['meta']['Unit'],
                          'Description': 'Rotational part of {} ({})'.format(F['label'],F['meta']['Description']),
                          'Creator':     'dadf5.py:add_rotational_part v{}'.format(version)
                          }
                 }
    def add_rotational_part(self,F):
        """
        Add rotational part of a deformation gradient.

        Parameters
        ----------
        F : str, optional
          Label of deformation gradient dataset.

        """
        self.__add_generic_pointwise(self._add_rotational_part,{'F':F})


    @staticmethod
    def _add_spherical(T):
        if not T['data'].shape[1:] == (3,3):
            raise ValueError

        return {
                'data':  mechanics.spherical_part(T['data']),
                'label': 'p_{}'.format(T['label']),
                'meta':  {
                          'Unit':        T['meta']['Unit'],
                          'Description': 'Spherical component of tensor {} ({})'.format(T['label'],T['meta']['Description']),
                          'Creator':     'dadf5.py:add_spherical v{}'.format(version)
                          }
                 }
    def add_spherical(self,T):
        """
        Add the spherical (hydrostatic) part of a tensor.

        Parameters
        ----------
        T : str
          Label of tensor dataset.

        """
        self.__add_generic_pointwise(self._add_spherical,{'T':T})


    @staticmethod
    def _add_strain_tensor(F,t,m):
        if not F['data'].shape[1:] == (3,3):
            raise ValueError

        return {
                'data':  mechanics.strain_tensor(F['data'],t,m),
                'label': 'epsilon_{}^{}({})'.format(t,m,F['label']),
                'meta':  {
                          'Unit':        F['meta']['Unit'],
                          'Description': 'Strain tensor of {} ({})'.format(F['label'],F['meta']['Description']),
                          'Creator':     'dadf5.py:add_strain_tensor v{}'.format(version)
                          }
                 }
    def add_strain_tensor(self,F='F',t='V',m=0.0):
        """
        Add strain tensor of a deformation gradient.

        For details refer to damask.mechanics.strain_tensor

        Parameters
        ----------
        F : str, optional
          Label of deformation gradient dataset. Defaults to ‘F’.
        t : {‘V’, ‘U’}, optional
          Type of the polar decomposition, ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.
          Defaults to ‘V’.
        m : float, optional
          Order of the strain calculation. Defaults to ‘0.0’.

        """
        self.__add_generic_pointwise(self._add_strain_tensor,{'F':F},{'t':t,'m':m})


    @staticmethod
    def _add_stretch_tensor(F,t):
        if not F['data'].shape[1:] == (3,3):
            raise ValueError

        return {
                'data':  mechanics.left_stretch(F['data']) if t == 'V' else mechanics.right_stretch(F['data']),
                'label': '{}({})'.format(t,F['label']),
                'meta':  {
                          'Unit':        F['meta']['Unit'],
                          'Description': '{} stretch tensor of {} ({})'.format('Left' if t == 'V' else 'Right',
                                                                               F['label'],F['meta']['Description']),
                          'Creator':     'dadf5.py:add_stretch_tensor v{}'.format(version)
                          }
                 }
    def add_stretch_tensor(self,F='F',t='V'):
        """
        Add stretch tensor of a deformation gradient.

        Parameters
        ----------
        F : str, optional
          Label of deformation gradient dataset. Defaults to ‘F’.
        t : {‘V’, ‘U’}, optional
          Type of the polar decomposition, ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.
          Defaults to ‘V’.

        """
        self.__add_generic_pointwise(self._add_stretch_tensor,{'F':F},{'t':t})


    def job(self,group,func,datasets,args,lock):
        try:
            d = self._read(group,datasets,lock)
            r = func(**d,**args)
            return [group,r]
        except Exception as err:
            print('Error during calculation: {}.'.format(err))
            return None


    def _read(self,group,datasets,lock):
      datasets_in = {}
      lock.acquire()
      with h5py.File(self.fname,'r') as f:
        for k,v in datasets.items():
          loc  = f[group+'/'+v]
          datasets_in[k]={'data':loc[()],
                          'label':v,
                          'meta':{k2:v2.decode() for k2,v2 in loc.attrs.items()}}
      lock.release()
      return datasets_in

    def __add_generic_pointwise(self,func,datasets,args={}):

      env = Environment()
      N_threads = int(env.options['DAMASK_NUM_THREADS'])
      pool = multiprocessing.Pool(N_threads)
      m = multiprocessing.Manager()
      lock = m.Lock()

      groups = self.groups_with_datasets(datasets.values())
      default_arg = partial(self.job,func=func,datasets=datasets,args=args,lock=lock)
      util.progressBar(iteration=0,total=len(groups)-1)
      for i,result in enumerate(pool.imap_unordered(default_arg,groups)):
          util.progressBar(iteration=i,total=len(groups)-1)
          if not result: continue
          lock.acquire()
          with h5py.File(self.fname, 'a') as f:
              try:
                  dataset = f[result[0]].create_dataset(result[1]['label'],data=result[1]['data'])
                  for l,v in result[1]['meta'].items():
                      dataset.attrs[l]=v.encode()
              except OSError as err:
                  print('Could not add dataset: {}.'.format(err))
          lock.release()
      pool.close()
      pool.join()


    def to_vtk(self,labels,mode='cell'):
        """
        Export to vtk cell/point data.

        Parameters
        ----------
        labels : str or list of
          Labels of the datasets to be exported.
        mode : str, either 'cell' or 'point'
          Export in cell format or point format.
          Defaults to 'cell'.

        """
        if mode.lower()=='cell':

          if self.structured:

            coordArray = [vtk.vtkDoubleArray(),vtk.vtkDoubleArray(),vtk.vtkDoubleArray()]
            for dim in [0,1,2]:
              for c in np.linspace(0,self.size[dim],1+self.grid[dim]):
                coordArray[dim].InsertNextValue(c)

            vtk_geom = vtk.vtkRectilinearGrid()
            vtk_geom.SetDimensions(*(self.grid+1))
            vtk_geom.SetXCoordinates(coordArray[0])
            vtk_geom.SetYCoordinates(coordArray[1])
            vtk_geom.SetZCoordinates(coordArray[2])

          else:

            nodes = vtk.vtkPoints()
            with h5py.File(self.fname,'r') as f:
              nodes.SetData(numpy_support.numpy_to_vtk(f['/geometry/x_n'][()],deep=True))

              vtk_geom = vtk.vtkUnstructuredGrid()
              vtk_geom.SetPoints(nodes)
              vtk_geom.Allocate(f['/geometry/T_c'].shape[0])

              if self.version_major == 0 and self.version_minor <= 5:
                vtk_type = vtk.VTK_HEXAHEDRON
                n_nodes = 8
              else:
                if   f['/geometry/T_c'].attrs['VTK_TYPE'] == b'TRIANGLE':
                  vtk_type = vtk.VTK_TRIANGLE
                  n_nodes = 3
                elif f['/geometry/T_c'].attrs['VTK_TYPE'] == b'QUAD':
                  vtk_type = vtk.VTK_QUAD
                  n_nodes = 4
                elif f['/geometry/T_c'].attrs['VTK_TYPE'] == b'TETRA':                                  # not tested
                  vtk_type = vtk.VTK_TETRA
                  n_nodes = 4
                elif f['/geometry/T_c'].attrs['VTK_TYPE'] == b'HEXAHEDRON':
                  vtk_type = vtk.VTK_HEXAHEDRON
                  n_nodes = 8

              for i in f['/geometry/T_c']:
                vtk_geom.InsertNextCell(vtk_type,n_nodes,i-1)

        elif mode.lower()=='point':
          Points   = vtk.vtkPoints()
          Vertices = vtk.vtkCellArray()
          for c in self.cell_coordinates():
            pointID = Points.InsertNextPoint(c)
            Vertices.InsertNextCell(1)
            Vertices.InsertCellPoint(pointID)

          vtk_geom = vtk.vtkPolyData()
          vtk_geom.SetPoints(Points)
          vtk_geom.SetVerts(Vertices)
          vtk_geom.Modified()

        N_digits = int(np.floor(np.log10(int(self.increments[-1][3:]))))+1

        for i,inc in enumerate(self.iter_visible('increments')):
          vtk_data = []

          materialpoints_backup = self.selection['materialpoints'].copy()
          self.set_visible('materialpoints',False)
          for label in (labels if isinstance(labels,list) else [labels]):
            for p in self.iter_visible('con_physics'):
              if p != 'generic':
                for c in self.iter_visible('constituents'):
                  x = self.get_dataset_location(label)
                  if len(x) == 0:
                    continue
                  array = self.read_dataset(x,0)
                  shape = [array.shape[0],np.product(array.shape[1:])]
                  vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),
                                  deep=True,array_type= vtk.VTK_DOUBLE))
                  vtk_data[-1].SetName('1_'+x[0].split('/',1)[1]) #ToDo: hard coded 1!
                  vtk_geom.GetCellData().AddArray(vtk_data[-1])

              else:
                x = self.get_dataset_location(label)
                if len(x) == 0:
                  continue
                array = self.read_dataset(x,0)
                shape = [array.shape[0],np.product(array.shape[1:])]
                vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),
                                deep=True,array_type= vtk.VTK_DOUBLE))
                ph_name = re.compile(r'(?<=(constituent\/))(.*?)(?=(generic))')                         # identify  phase name
                dset_name = '1_' + re.sub(ph_name,r'',x[0].split('/',1)[1])                             # removing phase name
                vtk_data[-1].SetName(dset_name)
                vtk_geom.GetCellData().AddArray(vtk_data[-1])

          self.set_visible('materialpoints',materialpoints_backup)

          constituents_backup = self.selection['constituents'].copy()
          self.set_visible('constituents',False)
          for label in (labels if isinstance(labels,list) else [labels]):
            for p in self.iter_visible('mat_physics'):
              if p != 'generic':
                for m in self.iter_visible('materialpoints'):
                  x = self.get_dataset_location(label)
                  if len(x) == 0:
                    continue
                  array = self.read_dataset(x,0)
                  shape = [array.shape[0],np.product(array.shape[1:])]
                  vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),
                                  deep=True,array_type= vtk.VTK_DOUBLE))
                  vtk_data[-1].SetName('1_'+x[0].split('/',1)[1]) #ToDo: why 1_?
                  vtk_geom.GetCellData().AddArray(vtk_data[-1])
              else:
                x = self.get_dataset_location(label)
                if len(x) == 0:
                  continue
                array = self.read_dataset(x,0)
                shape = [array.shape[0],np.product(array.shape[1:])]
                vtk_data.append(numpy_support.numpy_to_vtk(num_array=array.reshape(shape),
                                deep=True,array_type= vtk.VTK_DOUBLE))
                vtk_data[-1].SetName('1_'+x[0].split('/',1)[1])
                vtk_geom.GetCellData().AddArray(vtk_data[-1])
          self.set_visible('constituents',constituents_backup)

          if mode.lower()=='cell':
            writer = vtk.vtkXMLRectilinearGridWriter() if self.structured else \
                     vtk.vtkXMLUnstructuredGridWriter()
            x = self.get_dataset_location('u_n')
            vtk_data.append(numpy_support.numpy_to_vtk(num_array=self.read_dataset(x,0),
                            deep=True,array_type=vtk.VTK_DOUBLE))
            vtk_data[-1].SetName('u')
            vtk_geom.GetPointData().AddArray(vtk_data[-1])
          elif mode.lower()=='point':
            writer = vtk.vtkXMLPolyDataWriter()


          file_out = '{}_inc{}.{}'.format(os.path.splitext(os.path.basename(self.fname))[0],
                                          inc[3:].zfill(N_digits),
                                          writer.GetDefaultFileExtension())

          writer.SetCompressorTypeToZLib()
          writer.SetDataModeToBinary()
          writer.SetFileName(file_out)
          writer.SetInputData(vtk_geom)

          writer.Write()
