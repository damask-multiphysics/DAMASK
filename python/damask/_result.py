import multiprocessing as mp
import re
import glob
import os
import datetime
import xml.etree.ElementTree as ET
import xml.dom.minidom
from pathlib import Path
from functools import partial
from collections import defaultdict
from collections.abc import Iterable

import h5py
import numpy as np
import numpy.ma as ma
from numpy.lib import recfunctions as rfn

import damask
from . import VTK
from . import Orientation
from . import grid_filters
from . import mechanics
from . import tensor
from . import util

h5py3 = h5py.__version__[0] == '3'

def _read(handle,path):
    metadata = {k:(v if h5py3 else v.decode()) for k,v in handle[path].attrs.items()}
    dtype = np.dtype(handle[path].dtype,metadata=metadata)
    return np.array(handle[path],dtype=dtype)


class Result:
    """
    Manipulate and read DADF5 files.

    DADF5 (DAMASK HDF5) files contain DAMASK results.
    The group/folder structure reflects the input data
    in material.yaml.
    """

    def __init__(self,fname):
        """
        New result view bound to a HDF5 file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Name of the DADF5 file to be opened.

        """
        with h5py.File(fname,'r') as f:

            self.version_major = f.attrs['DADF5_version_major']
            self.version_minor = f.attrs['DADF5_version_minor']

            if self.version_major != 0 or not 7 <= self.version_minor <= 12:
                raise TypeError(f'Unsupported DADF5 version {self.version_major}.{self.version_minor}')

            self.structured = 'grid' in f['geometry'].attrs.keys() or \
                              'cells' in f['geometry'].attrs.keys()

            if self.structured:
                try:
                    self.cells  = f['geometry'].attrs['cells']
                except KeyError:
                    self.cells  = f['geometry'].attrs['grid']
                self.size   = f['geometry'].attrs['size']
                self.origin = f['geometry'].attrs['origin']

            r=re.compile('inc[0-9]+' if self.version_minor < 12 else 'increment_[0-9]+')
            increments_unsorted = {int(i[10:]):i for i in f.keys() if r.match(i)}
            self.increments     = [increments_unsorted[i] for i in sorted(increments_unsorted)]
            self.times          = [round(f[i].attrs['time/s'],12) for i in self.increments] if self.version_minor < 12 else \
                                  [round(f[i].attrs['t/s'],12) for i in self.increments]

            grp = 'mapping' if self.version_minor < 12 else 'cell_to'

            self.N_materialpoints, self.N_constituents = np.shape(f[f'{grp}/phase'])

            self.homogenizations  = [m.decode() for m in np.unique(f[f'{grp}/homogenization']
                                                                    ['Name' if self.version_minor < 12 else 'label'])]
            self.phases           = [c.decode() for c in np.unique(f[f'{grp}/phase']
                                                                    ['Name' if self.version_minor < 12 else 'label'])]

            self.out_type_ph = []
            for c in self.phases:
                self.out_type_ph += f['/'.join([self.increments[0],'phase',c])].keys()
            self.out_type_ph = list(set(self.out_type_ph))                                          # make unique

            self.out_type_ho = []
            for m in self.homogenizations:
                self.out_type_ho += f['/'.join([self.increments[0],'homogenization',m])].keys()
            self.out_type_ho = list(set(self.out_type_ho))                                          # make unique

        self.visible = {'increments':      self.increments,
                        'phases':          self.phases,
                        'homogenizations': self.homogenizations,
                        'out_type_ph':     self.out_type_ph,
                        'out_type_ho':     self.out_type_ho
                       }

        self.fname = Path(fname).absolute()

        self._allow_modification = False


    def __repr__(self):
        """Show summary of file content."""
        visible_increments = self.visible['increments']

        self.view('increments',visible_increments[0:1])
        first = self.list_data()

        self.view('increments',visible_increments[-1:])
        last  = '' if len(visible_increments) < 2 else self.list_data()

        self.view('increments',visible_increments)

        in_between = '' if len(visible_increments) < 3 else \
                     ''.join([f'\n{inc}\n  ...\n' for inc in visible_increments[1:-1]])

        return util.srepr(first + in_between + last)


    def _manage_view(self,action,what,datasets):
        """
        Manages the visibility of the groups.

        Parameters
        ----------
        action : str
            Select from 'set', 'add', and 'del'.
        what : str
            Attribute to change (must be from self.visible).
        datasets : list of str or bool
            Name of datasets as list; supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to [].

        """
        def natural_sort(key):
            convert = lambda text: int(text) if text.isdigit() else text
            return [ convert(c) for c in re.split('([0-9]+)', key) ]

        # allow True/False and string arguments
        if  datasets is True:
            datasets = ['*']
        elif datasets is False:
            datasets = []
        choice = datasets if hasattr(datasets,'__iter__') and not isinstance(datasets,str) else \
                [datasets]

        inc = 'inc' if self.version_minor < 12 else 'increment_' # compatibility hack
        if   what == 'increments':
            choice = [c if isinstance(c,str) and c.startswith(inc) else
                      f'{inc}{c}' for c in choice]
        elif what == 'times':
            what = 'increments'
            if choice == ['*']:
                choice = self.increments
            else:
                iterator = map(float,choice)
                choice = []
                for c in iterator:
                    idx = np.searchsorted(self.times,c)
                    if idx >= len(self.times): continue
                    if   np.isclose(c,self.times[idx]):
                        choice.append(self.increments[idx])
                    elif np.isclose(c,self.times[idx+1]):
                        choice.append(self.increments[idx+1])

        valid = [e for e_ in [glob.fnmatch.filter(getattr(self,what),s) for s in choice] for e in e_]
        existing = set(self.visible[what])

        if   action == 'set':
            self.visible[what] = valid
        elif action == 'add':
            add = existing.union(valid)
            add_sorted = sorted(add, key=natural_sort)
            self.visible[what] = add_sorted
        elif action == 'del':
            diff = existing.difference(valid)
            diff_sorted = sorted(diff, key=natural_sort)
            self.visible[what] = diff_sorted


    def _get_attribute(self,path,attr):
        """
        Get the attribute of a dataset.

        Parameters
        ----------
        Path : str
            Path to the dataset.
        attr : str
            Name of the attribute to get.

        Returns
        -------
        attr at path, str or None.
            The requested attribute, None if not found.

        """
        with h5py.File(self.fname,'r') as f:
            try:
                return f[path].attrs[attr] if h5py3 else f[path].attrs[attr].decode()
            except KeyError:
                return None


    def allow_modification(self):
        """Allow to overwrite existing data."""
        print(util.warn('Warning: Modification of existing datasets allowed!'))
        self._allow_modification = True

    def disallow_modification(self):
        """Disallow to overwrite existing data (default case)."""
        self._allow_modification = False


    def increments_in_range(self,start,end):
        """
        Select all increments within a given range.

        Parameters
        ----------
        start : int or str
            Start increment.
        end : int or str
            End increment.

        """
        # compatibility hack
        ln = 3 if self.version_minor < 12 else 10
        selected = []
        for i,inc in enumerate([int(i[ln:]) for i in self.increments]):
            s,e = map(lambda x: int(x[ln:] if isinstance(x,str) and x.startswith('inc') else x), (start,end))
            if s <= inc <= e:
                selected.append(self.increments[i])
        return selected


    def times_in_range(self,start,end):
        """
        Select all increments within a given time range.

        Parameters
        ----------
        start : float
            Time of start increment.
        end : float
            Time of end increment.

        """
        selected = []
        for i,time in enumerate(self.times):
            if start <= time <= end:
                selected.append(self.times[i])
        return selected


    def iterate(self,what):
        """
        Iterate over visible items and view them independently.

        Parameters
        ----------
        what : str
            Attribute to change (must be from self.visible).

        """
        datasets = self.visible[what]
        last_view = datasets.copy()
        for dataset in datasets:
            if last_view != self.visible[what]:
                self._manage_view('set',what,datasets)
                raise Exception
            self._manage_view('set',what,dataset)
            last_view = self.visible[what]
            yield dataset
        self._manage_view('set',what,datasets)


    def view(self,what,datasets):
        """
        Set view.

        Parameters
        ----------
        what : str
            Attribute to change (must be from self.visible).
        datasets : list of str or bool
            Name of datasets as list; supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to [].

        """
        self._manage_view('set',what,datasets)


    def view_more(self,what,datasets):
        """
        Add to view.

        Parameters
        ----------
        what : str
            Attribute to change (must be from self.visible).
        datasets : list of str or bool
            Name of datasets as list; supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to [].

        """
        self._manage_view('add',what,datasets)


    def view_less(self,what,datasets):
        """
        Delete from view.

        Parameters
        ----------
        what : str
            Attribute to change (must be from self.visible).
        datasets : list of str or bool
            Name of datasets as list; supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to [].

        """
        self._manage_view('del',what,datasets)


    def rename(self,name_old,name_new):
        """
        Rename dataset.

        Parameters
        ----------
        name_old : str
            Name of the dataset to be renamed.
        name_new : str
            New name of the dataset.

        """
        if self._allow_modification:
            with h5py.File(self.fname,'a') as f:
                for path_old in self.get_dataset_location(name_old):
                    path_new = os.path.join(os.path.dirname(path_old),name_new)
                    f[path_new] = f[path_old]
                    f[path_new].attrs['Renamed'] = f'Original name: {name_old}' if h5py3 else \
                                                   f'Original name: {name_old}'.encode()
                    del f[path_old]
        else:
            raise PermissionError('Rename operation not permitted')


    def groups_with_datasets(self,datasets):
        """
        Return groups that contain all requested datasets.

        Only groups within
          - inc*/phase/*/*
          - inc*/homogenization/*/*
          - inc*/geometry/*

        are considered as they contain user-relevant data.
        Single strings will be treated as list with one entry.

        Wild card matching is allowed, but the number of arguments needs to fit.

        Parameters
        ----------
            datasets : iterable or str or bool

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

        sets = datasets if isinstance(datasets,bool) or (hasattr(datasets,'__iter__') and not isinstance(datasets,str)) else \
              [datasets]

        groups = []

        with h5py.File(self.fname,'r') as f:
            for i in self.iterate('increments'):
                for o,p in zip(['phases','homogenizations'],['out_type_ph','out_type_ho']):
                    for oo in self.iterate(o):
                        for pp in self.iterate(p):
                            group = '/'.join([i,o[:-1],oo,pp])                                      # o[:-1]: plural/singular issue
                            if sets is True:
                                groups.append(group)
                            else:
                                if group in f.keys():
                                    match = [e for e_ in [glob.fnmatch.filter(f[group].keys(),s) for s in sets] for e in e_]
                                    if len(set(match)) == len(sets): groups.append(group)
        return groups


    def list_data(self):
        """Return information on all active datasets in the file."""
        # compatibility hack
        de = 'Description' if self.version_minor < 12 else 'description'
        un = 'Unit'        if self.version_minor < 12 else 'unit'
        message = ''
        with h5py.File(self.fname,'r') as f:
            for i in self.iterate('increments'):
                message += f'\n{i} ({self.times[self.increments.index(i)]}s)\n'
                for o,p in zip(['phases','homogenizations'],['out_type_ph','out_type_ho']):
                    message += f'  {o[:-1]}\n'
                    for oo in self.iterate(o):
                        message += f'    {oo}\n'
                        for pp in self.iterate(p):
                            message += f'      {pp}\n'
                            group = '/'.join([i,o[:-1],oo,pp])                                      # o[:-1]: plural/singular issue
                            for d in f[group].keys():
                                try:
                                    dataset = f['/'.join([group,d])]
                                    if un in dataset.attrs:
                                        unit = f" / {dataset.attrs[un]}" if h5py3 else \
                                               f" / {dataset.attrs[un].decode()}"
                                    else:
                                        unit = ''
                                    description = dataset.attrs[de] if h5py3 else \
                                                  dataset.attrs[de].decode()
                                    message += f'        {d}{unit}: {description}\n'
                                except KeyError:
                                    pass
        return message


    def get_dataset_location(self,label):
        """Return the location of all active datasets with given label."""
        path = []
        with h5py.File(self.fname,'r') as f:
            for i in self.iterate('increments'):
                k = '/'.join([i,'geometry',label])
                try:
                    f[k]
                    path.append(k)
                except KeyError:
                    pass
                for o,p in zip(['phases','homogenizations'],['out_type_ph','out_type_ho']):
                    for oo in self.iterate(o):
                        for pp in self.iterate(p):
                            k = '/'.join([i,o[:-1],oo,pp,label])
                            try:
                                f[k]
                                path.append(k)
                            except KeyError:
                                pass
        return path


    def enable_user_function(self,func):
        globals()[func.__name__]=func
        print(f'Function {func.__name__} enabled in add_calculation.')


    def read_dataset(self,path,c=0,plain=False):
        """
        Dataset for all points/cells.

        If more than one path is given, the dataset is composed of the individual contributions.

        Parameters
        ----------
        path : list of strings
            The name of the datasets to consider.
        c : int, optional
            The constituent to consider. Defaults to 0.
        plain: boolean, optional
            Convert into plain numpy datatype.
            Only relevant for compound datatype, e.g. the orientation.
            Defaults to False.

        """
        # compatibility hack
        name   = 'Name' if self.version_minor < 12 else 'label'
        member = 'Position' if self.version_minor < 12 else 'entry'
        grp    = 'mapping' if self.version_minor < 12 else 'cell_to'
        with h5py.File(self.fname,'r') as f:
            shape = (self.N_materialpoints,) + np.shape(f[path[0]])[1:]
            if len(shape) == 1: shape = shape +(1,)
            dataset = np.full(shape,np.nan,dtype=np.dtype(f[path[0]]))
            for pa in path:
                label = pa.split('/')[2]

                if pa.split('/')[1] == 'geometry':
                    dataset = np.array(f[pa])
                    continue

                p = np.where(f[f'{grp}/phase'][:,c][name] == str.encode(label))[0]
                if len(p)>0:
                    u = (f[f'{grp}/phase'][member][p,c])
                    a = np.array(f[pa])
                    if len(a.shape) == 1:
                        a=a.reshape([a.shape[0],1])
                    dataset[p,:] = a[u,:]

                p = np.where(f[f'{grp}/homogenization'][name] == str.encode(label))[0]
                if len(p)>0:
                    u = (f[f'{grp}/homogenization'][member][p.tolist()])
                    a = np.array(f[pa])
                    if len(a.shape) == 1:
                        a=a.reshape([a.shape[0],1])
                    dataset[p,:] = a[u,:]

        if plain and dataset.dtype.names is not None:
            return dataset.view(('float64',len(dataset.dtype.names)))
        else:
            return dataset

    @property
    def coordinates0_point(self):
        """Return initial coordinates of the cell centers."""
        if self.structured:
            return grid_filters.coordinates0_point(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_c'][()]

    @property
    def coordinates0_node(self):
        """Return initial coordinates of the cell centers."""
        if self.structured:
            return grid_filters.coordinates0_node(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_n'][()]


    @staticmethod
    def _add_absolute(x):
        return {
                'data':  np.abs(x['data']),
                'label': f'|{x["label"]}|',
                'meta':  {
                          'unit':        x['meta']['unit'],
                          'description': f"absolute value of {x['label']} ({x['meta']['description']})",
                          'creator':     'add_absolute'
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
        self._add_generic_pointwise(self._add_absolute,{'x':x})


    @staticmethod
    def _add_calculation(**kwargs):
        formula = kwargs['formula']
        for d in re.findall(r'#(.*?)#',formula):
            formula = formula.replace(f'#{d}#',f"kwargs['{d}']['data']")

        return {
                'data':  eval(formula),
                'label': kwargs['label'],
                'meta':  {
                          'unit':        kwargs['unit'],
                          'description': f"{kwargs['description']} (formula: {kwargs['formula']})",
                          'creator':     'add_calculation'
                          }
                 }
    def add_calculation(self,label,formula,unit='n/a',description=None):
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

        """
        dataset_mapping  = {d:d for d in set(re.findall(r'#(.*?)#',formula))}                       # datasets used in the formula
        args             = {'formula':formula,'label':label,'unit':unit,'description':description}
        self._add_generic_pointwise(self._add_calculation,dataset_mapping,args)


    @staticmethod
    def _add_stress_Cauchy(P,F):
        return {
                'data':  mechanics.stress_Cauchy(P['data'],F['data']),
                'label': 'sigma',
                'meta':  {
                          'unit':        P['meta']['unit'],
                          'description': "Cauchy stress calculated "
                                         f"from {P['label']} ({P['meta']['description']})"
                                         f" and {F['label']} ({F['meta']['description']})",
                          'creator':     'add_stress_Cauchy'
                          }
                }
    def add_stress_Cauchy(self,P='P',F='F'):
        """
        Add Cauchy stress calculated from first Piola-Kirchhoff stress and deformation gradient.

        Parameters
        ----------
        P : str, optional
            Label of the dataset containing the first Piola-Kirchhoff stress. Defaults to ‘P’.
        F : str, optional
            Label of the dataset containing the deformation gradient. Defaults to ‘F’.

        """
        self._add_generic_pointwise(self._add_stress_Cauchy,{'P':P,'F':F})


    @staticmethod
    def _add_determinant(T):
        return {
                'data':  np.linalg.det(T['data']),
                'label': f"det({T['label']})",
                'meta':  {
                          'unit':        T['meta']['unit'],
                          'description': f"determinant of tensor {T['label']} ({T['meta']['description']})",
                          'creator':     'add_determinant'
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
        self._add_generic_pointwise(self._add_determinant,{'T':T})


    @staticmethod
    def _add_deviator(T):
        return {
                'data':  tensor.deviatoric(T['data']),
                'label': f"s_{T['label']}",
                'meta':  {
                          'unit':        T['meta']['unit'],
                          'description': f"deviator of tensor {T['label']} ({T['meta']['description']})",
                          'creator':     'add_deviator'
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
        self._add_generic_pointwise(self._add_deviator,{'T':T})


    @staticmethod
    def _add_eigenvalue(T_sym,eigenvalue):
        if   eigenvalue == 'max':
            label,p = 'maximum',2
        elif eigenvalue == 'mid':
            label,p = 'intermediate',1
        elif eigenvalue == 'min':
            label,p = 'minimum',0

        return {
                'data': tensor.eigenvalues(T_sym['data'])[:,p],
                'label': f"lambda_{eigenvalue}({T_sym['label']})",
                'meta' : {
                          'unit':        T_sym['meta']['unit'],
                          'description': f"{label} eigenvalue of {T_sym['label']} ({T_sym['meta']['description']})",
                          'creator':     'add_eigenvalue'
                         }
                }
    def add_eigenvalue(self,T_sym,eigenvalue='max'):
        """
        Add eigenvalues of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Label of symmetric tensor dataset.
        eigenvalue : str, optional
            Eigenvalue. Select from ‘max’, ‘mid’, ‘min’. Defaults to ‘max’.

        """
        self._add_generic_pointwise(self._add_eigenvalue,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


    @staticmethod
    def _add_eigenvector(T_sym,eigenvalue):
        if   eigenvalue == 'max':
            label,p = 'maximum',2
        elif eigenvalue == 'mid':
            label,p = 'intermediate',1
        elif eigenvalue == 'min':
            label,p = 'minimum',0
        return {
                'data': tensor.eigenvectors(T_sym['data'])[:,p],
                'label': f"v_{eigenvalue}({T_sym['label']})",
                'meta' : {
                          'unit':        '1',
                          'description': f"eigenvector corresponding to {label} eigenvalue"
                                         f" of {T_sym['label']} ({T_sym['meta']['description']})",
                          'creator':     'add_eigenvector'
                         }
               }
    def add_eigenvector(self,T_sym,eigenvalue='max'):
        """
        Add eigenvector of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Label of symmetric tensor dataset.
        eigenvalue : str, optional
            Eigenvalue to which the eigenvector corresponds. Select from
            ‘max’, ‘mid’, ‘min’. Defaults to ‘max’.

        """
        self._add_generic_pointwise(self._add_eigenvector,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


    @staticmethod
    def _add_IPF_color(l,q):
        m = util.scale_to_coprime(np.array(l))
        try:
            lattice = {'fcc':'cF','bcc':'cI','hex':'hP'}[q['meta']['lattice']]
        except KeyError:
            lattice =  q['meta']['lattice']
        try:
            o = Orientation(rotation = (rfn.structured_to_unstructured(q['data'])),lattice=lattice)
        except ValueError:
            o = Orientation(rotation = q['data'],lattice=lattice)

        return {
                'data': np.uint8(o.IPF_color(l)*255),
                'label': 'IPFcolor_[{} {} {}]'.format(*m),
                'meta' : {
                          'unit':        '8-bit RGB',
                          'lattice':     q['meta']['lattice'],
                          'description': 'Inverse Pole Figure (IPF) colors along sample direction [{} {} {}]'.format(*m),
                          'creator':     'add_IPF_color'
                         }
               }
    def add_IPF_color(self,l,q='O'):
        """
        Add RGB color tuple of inverse pole figure (IPF) color.

        Parameters
        ----------
        l : numpy.array of shape (3)
            Lab frame direction for inverse pole figure.
        q : str
            Label of the dataset containing the crystallographic orientation as quaternions.
            Defaults to 'O'.

        """
        self._add_generic_pointwise(self._add_IPF_color,{'q':q},{'l':l})


    @staticmethod
    def _add_maximum_shear(T_sym):
        return {
                'data':  mechanics.maximum_shear(T_sym['data']),
                'label': f"max_shear({T_sym['label']})",
                'meta':  {
                          'unit':        T_sym['meta']['unit'],
                          'description': f"maximum shear component of {T_sym['label']} ({T_sym['meta']['description']})",
                          'creator':     'add_maximum_shear'
                          }
                 }
    def add_maximum_shear(self,T_sym):
        """
        Add maximum shear components of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Label of symmetric tensor dataset.

        """
        self._add_generic_pointwise(self._add_maximum_shear,{'T_sym':T_sym})


    @staticmethod
    def _add_equivalent_Mises(T_sym,kind):
        k = kind
        if k is None:
            if T_sym['meta']['unit'] == '1':
                k = 'strain'
            elif T_sym['meta']['unit'] == 'Pa':
                k = 'stress'
        if k not in ['stress', 'strain']:
            raise ValueError('invalid von Mises kind {kind}')

        return {
                'data':  (mechanics.equivalent_strain_Mises if k=='strain' else \
                          mechanics.equivalent_stress_Mises)(T_sym['data']),
                'label': f"{T_sym['label']}_vM",
                'meta':  {
                          'unit':        T_sym['meta']['unit'],
                          'description': f"Mises equivalent {k} of {T_sym['label']} ({T_sym['meta']['description']})",
                          'creator':     'add_Mises'
                          }
                }
    def add_equivalent_Mises(self,T_sym,kind=None):
        """
        Add the equivalent Mises stress or strain of a symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Label of symmetric tensorial stress or strain dataset.
        kind : {'stress', 'strain', None}, optional
            Kind of the von Mises equivalent. Defaults to None, in which case
            it is selected based on the unit of the dataset ('1' -> strain, 'Pa' -> stress').

        """
        self._add_generic_pointwise(self._add_equivalent_Mises,{'T_sym':T_sym},{'kind':kind})


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
                'label': f"|{x['label']}|_{o}",
                'meta':  {
                          'unit':        x['meta']['unit'],
                          'description': f"{o}-norm of {t} {x['label']} ({x['meta']['description']})",
                          'creator':     'add_norm'
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
        self._add_generic_pointwise(self._add_norm,{'x':x},{'ord':ord})


    @staticmethod
    def _add_stress_second_Piola_Kirchhoff(P,F):
        return {
                'data':  mechanics.stress_second_Piola_Kirchhoff(P['data'],F['data']),
                'label': 'S',
                'meta':  {
                          'unit':        P['meta']['unit'],
                          'description': "second Piola-Kirchhoff stress calculated "
                                         f"from {P['label']} ({P['meta']['description']})"
                                         f" and {F['label']} ({F['meta']['description']})",
                          'creator':     'add_stress_second_Piola_Kirchhoff'
                          }
                }
    def add_stress_second_Piola_Kirchhoff(self,P='P',F='F'):
        """
        Add second Piola-Kirchhoff stress calculated from first Piola-Kirchhoff stress and deformation gradient.

        Parameters
        ----------
        P : str, optional
            Label of first Piola-Kirchhoff stress dataset. Defaults to ‘P’.
        F : str, optional
            Label of deformation gradient dataset. Defaults to ‘F’.

        """
        self._add_generic_pointwise(self._add_stress_second_Piola_Kirchhoff,{'P':P,'F':F})


# The add_pole functionality needs discussion.
# The new Crystal object can perform such a calculation but the outcome depends on the lattice parameters
# as well as on whether a direction or plane is concerned (see the DAMASK_examples/pole_figure notebook).
# Below code appears to be too simplistic.

    # @staticmethod
    # def _add_pole(q,p,polar):
    #     pole      = np.array(p)
    #     unit_pole = pole/np.linalg.norm(pole)
    #     m         = util.scale_to_coprime(pole)
    #     rot       = Rotation(q['data'].view(np.double).reshape(-1,4))
    #
    #     rotatedPole = rot @ np.broadcast_to(unit_pole,rot.shape+(3,))                               # rotate pole according to crystal orientation
    #     xy = rotatedPole[:,0:2]/(1.+abs(unit_pole[2]))                                              # stereographic projection
    #     coords = xy if not polar else \
    #              np.block([np.sqrt(xy[:,0:1]*xy[:,0:1]+xy[:,1:2]*xy[:,1:2]),np.arctan2(xy[:,1:2],xy[:,0:1])])
    #     return {
    #             'data': coords,
    #             'label': 'p^{}_[{} {} {})'.format(u'rφ' if polar else 'xy',*m),
    #             'meta' : {
    #                       'unit':        '1',
    #                       'description': '{} coordinates of stereographic projection of pole (direction/plane) in crystal frame'\
    #                                      .format('Polar' if polar else 'Cartesian'),
    #                       'creator':     'add_pole'
    #                      }
    #            }
    # def add_pole(self,q,p,polar=False):
    #     """
    #     Add coordinates of stereographic projection of given pole in crystal frame.
    #
    #     Parameters
    #     ----------
    #     q : str
    #         Label of the dataset containing the crystallographic orientation as quaternions.
    #     p : numpy.array of shape (3)
    #         Crystallographic direction or plane.
    #     polar : bool, optional
    #         Give pole in polar coordinates. Defaults to False.
    #
    #     """
    #     self._add_generic_pointwise(self._add_pole,{'q':q},{'p':p,'polar':polar})


    @staticmethod
    def _add_rotation(F):
        return {
                'data':  mechanics.rotation(F['data']).as_matrix(),
                'label': f"R({F['label']})",
                'meta':  {
                          'unit':        F['meta']['unit'],
                          'description': f"rotational part of {F['label']} ({F['meta']['description']})",
                          'creator':     'add_rotation'
                          }
                 }
    def add_rotation(self,F):
        """
        Add rotational part of a deformation gradient.

        Parameters
        ----------
        F : str, optional
            Label of deformation gradient dataset.

        """
        self._add_generic_pointwise(self._add_rotation,{'F':F})


    @staticmethod
    def _add_spherical(T):
        return {
                'data':  tensor.spherical(T['data'],False),
                'label': f"p_{T['label']}",
                'meta':  {
                          'unit':        T['meta']['unit'],
                          'description': f"spherical component of tensor {T['label']} ({T['meta']['description']})",
                          'creator':     'add_spherical'
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
        self._add_generic_pointwise(self._add_spherical,{'T':T})


    @staticmethod
    def _add_strain(F,t,m):
        return {
                'data':  mechanics.strain(F['data'],t,m),
                'label': f"epsilon_{t}^{m}({F['label']})",
                'meta':  {
                          'unit':        F['meta']['unit'],
                          'description': f"strain tensor of {F['label']} ({F['meta']['description']})",
                          'creator':     'add_strain'
                          }
                 }
    def add_strain(self,F='F',t='V',m=0.0):
        """
        Add strain tensor of a deformation gradient.

        For details refer to damask.mechanics.strain

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
        self._add_generic_pointwise(self._add_strain,{'F':F},{'t':t,'m':m})


    @staticmethod
    def _add_stretch_tensor(F,t):
        return {
                'data':  (mechanics.stretch_left if t.upper() == 'V' else mechanics.stretch_right)(F['data']),
                'label': f"{t}({F['label']})",
                'meta':  {
                          'unit':        F['meta']['unit'],
                          'description': '{} stretch tensor of {} ({})'.format('left' if t.upper() == 'V' else 'right',
                                                                               F['label'],F['meta']['description']),
                          'creator':     'add_stretch_tensor'
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
        self._add_generic_pointwise(self._add_stretch_tensor,{'F':F},{'t':t})


    def _job(self,group,func,datasets,args,lock):
        """Execute job for _add_generic_pointwise."""
        try:
            datasets_in = {}
            lock.acquire()
            with h5py.File(self.fname,'r') as f:
                for arg,label in datasets.items():
                    loc  = f[group+'/'+label]
                    datasets_in[arg]={'data' :loc[()],
                                      'label':label,
                                      'meta': {k:(v if h5py3 else v.decode()) for k,v in loc.attrs.items()}}
            lock.release()
            r = func(**datasets_in,**args)
            return [group,r]
        except Exception as err:
            print(f'Error during calculation: {err}.')
            return None


    def _add_generic_pointwise(self,func,datasets,args={}):
        """
        General function to add pointwise data.

        Parameters
        ----------
        func : function
            Callback function that calculates a new dataset from one or
            more datasets per HDF5 group.
        datasets : dictionary
            Details of the datasets to be used: label (in HDF5 file) and
            arg (argument to which the data is parsed in func).
        args : dictionary, optional
            Arguments parsed to func.

        """
        chunk_size = 1024**2//8
        pool = mp.Pool(int(os.environ.get('OMP_NUM_THREADS',1)))
        lock = mp.Manager().Lock()

        groups = self.groups_with_datasets(datasets.values())
        if len(groups) == 0:
            print('No matching dataset found, no data was added.')
            return

        default_arg = partial(self._job,func=func,datasets=datasets,args=args,lock=lock)

        for result in util.show_progress(pool.imap_unordered(default_arg,groups),len(groups)):
            if not result:
                continue
            lock.acquire()
            with h5py.File(self.fname, 'a') as f:
                try:
                    if self._allow_modification and result[0]+'/'+result[1]['label'] in f:
                        dataset = f[result[0]+'/'+result[1]['label']]
                        dataset[...] = result[1]['data']
                        dataset.attrs['overwritten'] = True
                    else:
                        if result[1]['data'].size >= chunk_size*2:
                            shape  = result[1]['data'].shape
                            chunks = (chunk_size//np.prod(shape[1:]),)+shape[1:]
                            dataset = f[result[0]].create_dataset(result[1]['label'],data=result[1]['data'],
                                                                  maxshape=shape, chunks=chunks,
                                                                  compression='gzip', compression_opts=6,
                                                                  shuffle=True,fletcher32=True)
                        else:
                            dataset = f[result[0]].create_dataset(result[1]['label'],data=result[1]['data'])

                    now = datetime.datetime.now().astimezone()
                    dataset.attrs['created'] = now.strftime('%Y-%m-%d %H:%M:%S%z') if h5py3 else \
                                               now.strftime('%Y-%m-%d %H:%M:%S%z').encode()

                    for l,v in result[1]['meta'].items():
                        dataset.attrs[l.lower()]=v if h5py3 else v.encode()
                    creator = dataset.attrs['creator'] if h5py3 else \
                              dataset.attrs['creator'].decode()
                    dataset.attrs['creator'] = f"damask.Result.{creator} v{damask.version}" if h5py3 else \
                                               f"damask.Result.{creator} v{damask.version}".encode()

                except (OSError,RuntimeError) as err:
                    print(f'Could not add dataset: {err}.')
            lock.release()

        pool.close()
        pool.join()


    def save_XDMF(self):
        """
        Write XDMF file to directly visualize data in DADF5 file.

        The view is not taken into account, i.e. the content of the
        whole file will be included.
        """
        # compatibility hack
        u = 'Unit' if self.version_minor < 12 else 'unit'
        if self.N_constituents != 1 or len(self.phases) != 1 or not self.structured:
            raise TypeError('XDMF output requires homogeneous grid')


        attribute_type_map = defaultdict(lambda:'Matrix', ( ((),'Scalar'), ((3,),'Vector'), ((3,3),'Tensor')) )

        def number_type_map(dtype):
            if dtype in np.sctypes['int']:   return 'Int'
            if dtype in np.sctypes['uint']:  return 'UInt'
            if dtype in np.sctypes['float']: return 'Float'


        xdmf=ET.Element('Xdmf')
        xdmf.attrib={'Version':  '2.0',
                     'xmlns:xi': 'http://www.w3.org/2001/XInclude'}

        domain=ET.SubElement(xdmf, 'Domain')

        collection = ET.SubElement(domain, 'Grid')
        collection.attrib={'GridType':       'Collection',
                           'CollectionType': 'Temporal'}

        time = ET.SubElement(collection, 'Time')
        time.attrib={'TimeType': 'List'}

        time_data = ET.SubElement(time, 'DataItem')
        time_data.attrib={'Format':     'XML',
                          'NumberType': 'Float',
                          'Dimensions': f'{len(self.times)}'}
        time_data.text = ' '.join(map(str,self.times))

        attributes = []
        data_items = []

        for inc in self.increments:

            grid=ET.SubElement(collection,'Grid')
            grid.attrib = {'GridType': 'Uniform',
                           'Name':      inc}

            topology=ET.SubElement(grid, 'Topology')
            topology.attrib={'TopologyType': '3DCoRectMesh',
                             'Dimensions':   '{} {} {}'.format(*self.cells+1)}

            geometry=ET.SubElement(grid, 'Geometry')
            geometry.attrib={'GeometryType':'Origin_DxDyDz'}

            origin=ET.SubElement(geometry, 'DataItem')
            origin.attrib={'Format':     'XML',
                           'NumberType': 'Float',
                           'Dimensions': '3'}
            origin.text="{} {} {}".format(*self.origin)

            delta=ET.SubElement(geometry, 'DataItem')
            delta.attrib={'Format':     'XML',
                          'NumberType': 'Float',
                          'Dimensions': '3'}
            delta.text="{} {} {}".format(*(self.size/self.cells))


            with h5py.File(self.fname,'r') as f:
                attributes.append(ET.SubElement(grid, 'Attribute'))
                attributes[-1].attrib={'Name':          'u / m',
                                       'Center':        'Node',
                                       'AttributeType': 'Vector'}
                data_items.append(ET.SubElement(attributes[-1], 'DataItem'))
                data_items[-1].attrib={'Format':     'HDF',
                                       'Precision':  '8',
                                       'Dimensions': '{} {} {} 3'.format(*(self.cells+1))}
                data_items[-1].text=f'{os.path.split(self.fname)[1]}:/{inc}/geometry/u_n'

                for o,p in zip(['phases','homogenizations'],['out_type_ph','out_type_ho']):
                    for oo in getattr(self,o):
                        for pp in getattr(self,p):
                            g = '/'.join([inc,o[:-1],oo,pp])
                            for l in f[g]:
                                name = '/'.join([g,l])
                                shape = f[name].shape[1:]
                                dtype = f[name].dtype

                                if dtype not in np.sctypes['int']+np.sctypes['uint']+np.sctypes['float']: continue
                                unit = f[name].attrs[u] if h5py3 else f[name].attrs[u].decode()

                                attributes.append(ET.SubElement(grid, 'Attribute'))
                                attributes[-1].attrib={'Name':          name.split('/',2)[2]+f' / {unit}',
                                                       'Center':       'Cell',
                                                       'AttributeType': attribute_type_map[shape]}
                                data_items.append(ET.SubElement(attributes[-1], 'DataItem'))
                                data_items[-1].attrib={'Format':     'HDF',
                                                       'NumberType': number_type_map(dtype),
                                                       'Precision':  f'{dtype.itemsize}',
                                                       'Dimensions': '{} {} {} {}'.format(*self.cells,1 if shape == () else
                                                                                                      np.prod(shape))}
                                data_items[-1].text=f'{os.path.split(self.fname)[1]}:{name}'

        with open(self.fname.with_suffix('.xdmf').name,'w',newline='\n') as f:
            f.write(xml.dom.minidom.parseString(ET.tostring(xdmf).decode()).toprettyxml())


    def save_VTK(self,labels=[],mode='cell'):
        """
        Export to vtk cell/point data.

        Parameters
        ----------
        labels : str or list of, optional
            Labels of the datasets to be exported.
        mode : str, either 'cell' or 'point'
            Export in cell format or point format.
            Defaults to 'cell'.

        """
        if mode.lower()=='cell':

            if self.structured:
                v = VTK.from_rectilinear_grid(self.cells,self.size,self.origin)
            else:
                with h5py.File(self.fname,'r') as f:
                    v = VTK.from_unstructured_grid(f['/geometry/x_n'][()],
                                                   f['/geometry/T_c'][()]-1,
                                                   f['/geometry/T_c'].attrs['VTK_TYPE'] if h5py3 else \
                                                   f['/geometry/T_c'].attrs['VTK_TYPE'].decode())
        elif mode.lower()=='point':
            v = VTK.from_poly_data(self.coordinates0_point)

        # compatibility hack
        ln = 3 if self.version_minor < 12 else 10

        N_digits = int(np.floor(np.log10(max(1,int(self.increments[-1][ln:])))))+1

        for inc in util.show_progress(self.iterate('increments'),len(self.visible['increments'])):

            viewed_backup_ho = self.visible['homogenizations'].copy()
            self.view('homogenizations',False)
            for label in (labels if isinstance(labels,list) else [labels]):
                for o in self.iterate('out_type_ph'):
                    for c in range(self.N_constituents):
                        prefix = '' if self.N_constituents == 1 else f'constituent{c}/'
                        if o not in ['mechanics', 'mechanical']:                                    # compatibility hack
                            for _ in self.iterate('phases'):
                                path = self.get_dataset_location(label)
                                if len(path) == 0:
                                    continue
                                array = self.read_dataset(path,c)
                                v.add(array,prefix+path[0].split('/',1)[1]+f' / {self._get_attribute(path[0],"unit")}')
                        else:
                            paths = self.get_dataset_location(label)
                            if len(paths) == 0:
                                continue
                            array = self.read_dataset(paths,c)
                            if self.version_minor < 12:
                                ph_name = re.compile(r'(?<=(phase\/))(.*?)(?=(mechanics))')         # identify  phase name
                            else:
                                ph_name = re.compile(r'(?<=(phase\/))(.*?)(?=(mechanical))')        # identify  phase name
                            dset_name = prefix+re.sub(ph_name,r'',paths[0].split('/',1)[1])         # remove phase name
                            v.add(array,dset_name+f' / {self._get_attribute(paths[0],"unit")}')
            self.view('homogenizations',viewed_backup_ho)

            viewed_backup_ph = self.visible['phases'].copy()
            self.view('phases',False)
            for label in (labels if isinstance(labels,list) else [labels]):
                for _ in self.iterate('out_type_ho'):
                    paths = self.get_dataset_location(label)
                    if len(paths) == 0:
                        continue
                    array = self.read_dataset(paths)
                    v.add(array,paths[0].split('/',1)[1]+f' / {self._get_attribute(paths[0],"unit")}')
            self.view('phases',viewed_backup_ph)

            u = self.read_dataset(self.get_dataset_location('u_n' if mode.lower() == 'cell' else 'u_p'))
            v.add(u,'u')

            v.save(f'{self.fname.stem}_inc{inc[ln:].zfill(N_digits)}')


    def read(self,labels,compress=True,strip=True):
        """
        Export data from file per phase/homogenization.

        The returned data structure reflects the group/folder structure
        in the DADF5 file.

        Parameters
        ----------
        labels : str or list of, optional
            Labels of the datasets to be read.
        compress : bool
            Squeeze out dictionaries that are not needed for a unique
            structure. This might be beneficial in the case of single
            constituent or single phase simulations or if only one
            time increment is considered. Defaults to 'True'.
        strip : bool
            Remove branches that contain no dataset. Defaults to
            'True'.

        """
        r = {}
        labels_ = set([labels] if isinstance(labels,str) else labels)

        with h5py.File(self.fname,'r') as f:
            for inc in util.show_progress(self.visible['increments']):
                r[inc] = {'phase':{},'homogenization':{},'geometry':{}}

                for la in labels_.intersection(f['/'.join((inc,'geometry'))].keys()):
                    r[inc]['geometry'][la] = _read(f,'/'.join((inc,'geometry',la)))

                for ty in ['phase','homogenization']:
                    for na in self.visible[ty+'s']:
                        r[inc][ty][na] = {}
                        for field in f['/'.join((inc,ty,na))].keys():
                            r[inc][ty][na][field] = {}
                            for la in labels_.intersection(f['/'.join((inc,ty,na,field))].keys()):
                                r[inc][ty][na][field][la] = _read(f,'/'.join((inc,ty,na,field,la)))

        if strip:    r = util.dict_strip(r)
        if compress: r = util.dict_compress(r)

        return r


    def place(self,labels,compress=True,strip=True,constituents=None,fill_float=0.0,fill_int=0):
        """
        Export data from file suitable sorted for spatial operations.

        The returned data structure reflects the group/folder structure
        in the DADF5 file. In the case of multi phase simulations, the
        data is merged from the individual phases/homogenizations.

        In the cases of a single constituent and single phase simulation
        this function is equivalent to `read`.

        Parameters
        ----------
        labels : str or list of, optional
            Labels of the datasets to be read.
        compress : bool
            Squeeze out dictionaries that are not needed for a unique
            structure. This might be beneficial in the case of single
            phase simulations or if only one time increment is
            considered. Defaults to 'True'.
        strip : bool
            Remove branches that contain no dataset. Defaults to
            'True'.
        constituents : int or list of, optional
            Constituents to consider. Defaults to 'None', in which case
            all constituents are considered.
        fill_float : float
            Fill value for non existent entries of floating point type.
            Defaults to 0.0.
        fill_int : int
            Fill value for non existent entries of integer type.
            Defaults to 0.
        """
        r = {}

        labels_ = set([labels] if isinstance(labels,str) else labels)
        if constituents is None:
            constituents_ = range(self.N_constituents)
        else:
            constituents_ = constituents if isinstance(constituents,Iterable) else [constituents]

        suffixes = [''] if self.N_constituents == 1 or len(constituents_) == 1 else \
                   [f'#{c}' for c in constituents_]

        grp  = 'mapping' if self.version_minor < 12 else 'cell_to'
        name = 'Name' if self.version_minor < 12 else 'label'
        member = 'member' if self.version_minor < 12 else 'entry'

        with h5py.File(self.fname,'r') as f:

            at_cell_ph = []
            in_data_ph = []
            for c in constituents_:
                at_cell_ph.append({label: np.where(f[os.path.join(grp,'phase')][:,c][name] == label.encode())[0] \
                                          for label in self.visible['phases']})
                in_data_ph.append({label: f[os.path.join(grp,'phase')][member][at_cell_ph[c][label]][...,0] \
                                          for label in self.visible['phases']})

            at_cell_ho = {label: np.where(f[os.path.join(grp,'homogenization')][:][name] == label.encode())[0] \
                                 for label in self.visible['homogenizations']}
            in_data_ho = {label: f[os.path.join(grp,'homogenization')][member][at_cell_ho[label]] \
                                 for label in self.visible['homogenizations']}

            for inc in util.show_progress(self.visible['increments']):
                r[inc] = {'phase':{},'homogenization':{},'geometry':{}}

                for la in labels_.intersection(f[os.path.join(inc,'geometry')].keys()):
                    r[inc]['geometry'][la] = _read(f,os.path.join(inc,'geometry',la))

                for ph in self.visible['phases']:
                    for field in f[os.path.join(inc,'phase',ph)].keys():
                        if field not in r[inc]['phase'].keys():
                            r[inc]['phase'][field] = {}

                        for la in labels_.intersection(f[os.path.join(inc,'phase',ph,field)].keys()):
                            data = ma.array(_read(f,os.path.join(inc,'phase',ph,field,la)))

                            if la+suffixes[0] not in r[inc]['phase'][field].keys():
                                container = np.empty((self.N_materialpoints,)+data.shape[1:],dtype=data.dtype)
                                fill_value = fill_float if data.dtype in np.sctypes['float'] else \
                                             fill_int
                                for c,suffix in zip(constituents_, suffixes):
                                    r[inc]['phase'][field][la+suffix] = \
                                        ma.array(container,fill_value=fill_value,mask=True)

                            for c,suffix in zip(constituents_, suffixes):
                                r[inc]['phase'][field][la+suffix][at_cell_ph[c][ph]] = data[in_data_ph[c][ph]]

                for ho in self.visible['homogenizations']:
                    for field in f[os.path.join(inc,'homogenization',ho)].keys():
                        if field not in r[inc]['homogenization'].keys():
                            r[inc]['homogenization'][field] = {}

                        for la in labels_.intersection(f[os.path.join(inc,'homogenization',ho,field)].keys()):
                            data = ma.array(_read(f,os.path.join(inc,'homogenization',ho,field,la)))

                            if la not in r[inc]['homogenization'][field].keys():
                                container = np.empty((self.N_materialpoints,)+data.shape[1:],dtype=data.dtype)
                                fill_value = fill_float if data.dtype in np.sctypes['float'] else \
                                             fill_int
                                r[inc]['homogenization'][field][la] = \
                                    ma.array(container,fill_value=fill_value,mask=True)

                            r[inc]['homogenization'][field][la][at_cell_ho[ho]] = data[in_data_ho[ho]]

        if strip:    r = util.dict_strip(r)
        if compress: r = util.dict_compress(r)
        return r
