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

import h5py
import numpy as np
from numpy.lib import recfunctions as rfn

import damask
from . import VTK
from . import Table
from . import Orientation
from . import grid_filters
from . import mechanics
from . import tensor
from . import util

h5py3 = h5py.__version__[0] == '3'

class Result:
    """
    Read and write to DADF5 files.

    DADF5 (DAMASK HDF5) files contain DAMASK results.
    """

    def __init__(self,fname):
        """
        Open an existing DADF5 file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Name of the DADF5 file to be opened.

        """
        with h5py.File(fname,'r') as f:

            self.version_major = f.attrs['DADF5_version_major']
            self.version_minor = f.attrs['DADF5_version_minor']

            if self.version_major != 0 or not 7 <= self.version_minor <= 10:
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

            r=re.compile('inc[0-9]+')
            increments_unsorted = {int(i[3:]):i for i in f.keys() if r.match(i)}
            self.increments     = [increments_unsorted[i] for i in sorted(increments_unsorted)]
            self.times          = [round(f[i].attrs['time/s'],12) for i in self.increments]

            self.N_materialpoints, self.N_constituents = np.shape(f['mapping/phase'])

            self.homogenizations  = [m.decode() for m in np.unique(f['mapping/homogenization']['Name'])]
            self.phases           = [c.decode() for c in np.unique(f['mapping/phase']['Name'])]

            self.out_type_ph = []
            for c in self.phases:
                self.out_type_ph += f['/'.join([self.increments[0],'phase',c])].keys()
            self.out_type_ph = list(set(self.out_type_ph))                                          # make unique

            self.out_type_ho = []
            for m in self.homogenizations:
                self.out_type_ho += f['/'.join([self.increments[0],'homogenization',m])].keys()
            self.out_type_ho = list(set(self.out_type_ho))                                          # make unique

        self.selection = {'increments':      self.increments,
                          'phases':          self.phases,
                          'homogenizations': self.homogenizations,
                          'out_type_ph':     self.out_type_ph,
                          'out_type_ho':     self.out_type_ho
                         }

        self.fname = Path(fname).absolute()

        self._allow_modification = False


    def __repr__(self):
        """Show summary of file content."""
        all_selected_increments = self.selection['increments']

        self.pick('increments',all_selected_increments[0:1])
        first = self.list_data()

        self.pick('increments',all_selected_increments[-1:])
        last  = '' if len(all_selected_increments) < 2 else self.list_data()

        self.pick('increments',all_selected_increments)

        in_between = '' if len(all_selected_increments) < 3 else \
                     ''.join([f'\n{inc}\n  ...\n' for inc in all_selected_increments[1:-2]])

        return util.srepr(first + in_between + last)


    def _manage_selection(self,action,what,datasets):
        """
        Manages the visibility of the groups.

        Parameters
        ----------
        action : str
            Select from 'set', 'add', and 'del'.
        what : str
            Attribute to change (must be from self.selection).
        datasets : list of str or bool
            Name of datasets as list, supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to []

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

        if   what == 'increments':
            choice = [c if isinstance(c,str) and c.startswith('inc') else
                      f'inc{c}' for c in choice]
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
        existing = set(self.selection[what])

        if   action == 'set':
            self.selection[what] = valid
        elif action == 'add':
            add = existing.union(valid)
            add_sorted = sorted(add, key=natural_sort)
            self.selection[what] = add_sorted
        elif action == 'del':
            diff = existing.difference(valid)
            diff_sorted = sorted(diff, key=natural_sort)
            self.selection[what] = diff_sorted


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
        """Disllow to overwrite existing data (default case)."""
        self._allow_modification = False


    def incs_in_range(self,start,end):
        """
        Select all increments within a given range.

        Parameters
        ----------
        start : int or str
            Start increment.
        end : int or str
            End increment.

        """
        selected = []
        for i,inc in enumerate([int(i[3:]) for i in self.increments]):
            s,e = map(lambda x: int(x[3:] if isinstance(x,str) and x.startswith('inc') else x), (start,end))
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
        Iterate over selection items by setting each one selected.

        Parameters
        ----------
        what : str
            Attribute to change (must be from self.selection).

        """
        datasets = self.selection[what]
        last_selection = datasets.copy()
        for dataset in datasets:
            if last_selection != self.selection[what]:
                self._manage_selection('set',what,datasets)
                raise Exception
            self._manage_selection('set',what,dataset)
            last_selection = self.selection[what]
            yield dataset
        self._manage_selection('set',what,datasets)


    def pick(self,what,datasets):
        """
        Set selection.

        Parameters
        ----------
        what : str
            attribute to change (must be from self.selection)
        datasets : list of str or bool
            name of datasets as list, supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to []

        """
        self._manage_selection('set',what,datasets)


    def pick_more(self,what,datasets):
        """
        Add to selection.

        Parameters
        ----------
        what : str
            attribute to change (must be from self.selection)
        datasets : list of str or bool
            name of datasets as list, supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to []

        """
        self._manage_selection('add',what,datasets)


    def pick_less(self,what,datasets):
        """
        Delete from selection.

        Parameters
        ----------
        what : str
            attribute to change (must be from self.selection)
        datasets : list of str or bool
            name of datasets as list, supports ? and * wildcards.
            True is equivalent to [*], False is equivalent to []

        """
        self._manage_selection('del',what,datasets)


    def rename(self,name_old,name_new):
        """
        Rename datasets.

        Parameters
        ----------
        name_old : str
            name of the datasets to be renamed
        name_new : str
            new name of the datasets

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


    def place(self,datasets,constituent=0,tagged=False,split=True):
        """
        Distribute datasets onto geometry and return Table or (split) dictionary of Tables.

        Must not mix nodal end cell data.

        Only data within
        - inc*/phase/*/*
        - inc*/homogenization/*/*
        - inc*/geometry/*
        are considered.

        Parameters
        ----------
          datasets : iterable or str
          constituent : int
              Constituent to consider for phase data
          tagged : bool
              tag Table.column name with '#constituent'
              defaults to False
          split : bool
              split Table by increment and return dictionary of Tables
              defaults to True

        """
        sets = datasets if hasattr(datasets,'__iter__') and not isinstance(datasets,str) else \
              [datasets]
        tag = f'#{constituent}' if tagged else ''
        tbl = {} if split else None
        inGeom = {}
        inData = {}
        with h5py.File(self.fname,'r') as f:
            for dataset in sets:
                for group in self.groups_with_datasets(dataset):
                    path = os.path.join(group,dataset)
                    inc,prop,name,cat,item = (path.split('/') + ['']*5)[:5]
                    key = '/'.join([prop,name+tag])
                    if key not in inGeom:
                        if prop == 'geometry':
                            inGeom[key] = inData[key] = np.arange(self.N_materialpoints)
                        elif prop == 'phase':
                            inGeom[key] = np.where(f['mapping/phase'][:,constituent]['Name'] == str.encode(name))[0]
                            inData[key] =          f['mapping/phase'][inGeom[key],constituent]['Position']
                        elif prop == 'homogenization':
                            inGeom[key] = np.where(f['mapping/homogenization']['Name'] == str.encode(name))[0]
                            inData[key] =          f['mapping/homogenization'][inGeom[key].tolist()]['Position']
                    shape = np.shape(f[path])
                    data = np.full((self.N_materialpoints,) + (shape[1:] if len(shape)>1 else (1,)),
                                   np.nan,
                                   dtype=np.dtype(f[path]))
                    data[inGeom[key]] = (f[path] if len(shape)>1 else np.expand_dims(f[path],1))[inData[key]]
                    path = (os.path.join(*([prop,name]+([cat] if cat else [])+([item] if item else []))) if split else path)+tag
                    if split:
                        try:
                            tbl[inc].add(path,data)
                        except KeyError:
                            tbl[inc] = Table(data.reshape(self.N_materialpoints,-1),{path:data.shape[1:]})
                    else:
                        try:
                            tbl.add(path,data)
                        except AttributeError:
                            tbl = Table(data.reshape(self.N_materialpoints,-1),{path:data.shape[1:]})

        return tbl


    def groups_with_datasets(self,datasets):
        """
        Return groups that contain all requested datasets.

        Only groups within
          - inc*/phase/*/*
          - inc*/homogenization/*/*
          - inc*/geometry/*

        are considered as they contain user-relevant data.
        Single strings will be treated as list with one entry.

        Wild card matching is allowed, but the number of arguments need to fit.

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
                                    if 'Unit' in dataset.attrs:
                                        unit = f" / {dataset.attrs['Unit']}" if h5py3 else \
                                               f" / {dataset.attrs['Unit'].decode()}"
                                    else:
                                        unit = ''
                                    description = dataset.attrs['Description'] if h5py3 else \
                                                  dataset.attrs['Description'].decode()
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
        with h5py.File(self.fname,'r') as f:
            shape = (self.N_materialpoints,) + np.shape(f[path[0]])[1:]
            if len(shape) == 1: shape = shape +(1,)
            dataset = np.full(shape,np.nan,dtype=np.dtype(f[path[0]]))
            for pa in path:
                label = pa.split('/')[2]

                if pa.split('/')[1] == 'geometry':
                    dataset = np.array(f[pa])
                    continue

                p = np.where(f['mapping/phase'][:,c]['Name'] == str.encode(label))[0]
                if len(p)>0:
                    u = (f['mapping/phase']['Position'][p,c])
                    a = np.array(f[pa])
                    if len(a.shape) == 1:
                        a=a.reshape([a.shape[0],1])
                    dataset[p,:] = a[u,:]

                p = np.where(f['mapping/homogenization']['Name'] == str.encode(label))[0]
                if len(p)>0:
                    u = (f['mapping/homogenization']['Position'][p.tolist()])
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
                          'Unit':        x['meta']['Unit'],
                          'Description': f"Absolute value of {x['label']} ({x['meta']['Description']})",
                          'Creator':     'add_absolute'
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
                          'Unit':        kwargs['unit'],
                          'Description': f"{kwargs['description']} (formula: {kwargs['formula']})",
                          'Creator':     'add_calculation'
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
                          'Unit':        P['meta']['Unit'],
                          'Description': "Cauchy stress calculated "
                                         f"from {P['label']} ({P['meta']['Description']})"
                                         f" and {F['label']} ({F['meta']['Description']})",
                          'Creator':     'add_stress_Cauchy'
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
                          'Unit':        T['meta']['Unit'],
                          'Description': f"Determinant of tensor {T['label']} ({T['meta']['Description']})",
                          'Creator':     'add_determinant'
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
                          'Unit':        T['meta']['Unit'],
                          'Description': f"Deviator of tensor {T['label']} ({T['meta']['Description']})",
                          'Creator':     'add_deviator'
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
            label,p = 'Maximum',2
        elif eigenvalue == 'mid':
            label,p = 'Intermediate',1
        elif eigenvalue == 'min':
            label,p = 'Minimum',0

        return {
                'data': tensor.eigenvalues(T_sym['data'])[:,p],
                'label': f"lambda_{eigenvalue}({T_sym['label']})",
                'meta' : {
                          'Unit':         T_sym['meta']['Unit'],
                          'Description': f"{label} eigenvalue of {T_sym['label']} ({T_sym['meta']['Description']})",
                          'Creator':     'add_eigenvalue'
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
                          'Unit':        '1',
                          'Description': f"Eigenvector corresponding to {label} eigenvalue"
                                         f" of {T_sym['label']} ({T_sym['meta']['Description']})",
                          'Creator':     'add_eigenvector'
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
    def _add_IPF_color(q,l):
        m = util.scale_to_coprime(np.array(l))
        try:
            lattice = {'fcc':'cF','bcc':'cI','hex':'hP'}[q['meta']['Lattice']]
        except KeyError:
            lattice =  q['meta']['Lattice']
        o = Orientation(rotation = (rfn.structured_to_unstructured(q['data'])),lattice=lattice)

        return {
                'data': np.uint8(o.IPF_color(l)*255),
                'label': 'IPFcolor_[{} {} {}]'.format(*m),
                'meta' : {
                          'Unit':        '8-bit RGB',
                          'Lattice':     q['meta']['Lattice'],
                          'Description': 'Inverse Pole Figure (IPF) colors along sample direction [{} {} {}]'.format(*m),
                          'Creator':     'add_IPF_color'
                         }
               }
    def add_IPF_color(self,q,l):
        """
        Add RGB color tuple of inverse pole figure (IPF) color.

        Parameters
        ----------
        q : str
            Label of the dataset containing the crystallographic orientation as quaternions.
        l : numpy.array of shape (3)
            Lab frame direction for inverse pole figure.

        """
        self._add_generic_pointwise(self._add_IPF_color,{'q':q},{'l':l})


    @staticmethod
    def _add_maximum_shear(T_sym):
        return {
                'data':  mechanics.maximum_shear(T_sym['data']),
                'label': f"max_shear({T_sym['label']})",
                'meta':  {
                          'Unit':        T_sym['meta']['Unit'],
                          'Description': f"Maximum shear component of {T_sym['label']} ({T_sym['meta']['Description']})",
                          'Creator':     'add_maximum_shear'
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
            if T_sym['meta']['Unit'] == '1':
                k = 'strain'
            elif T_sym['meta']['Unit'] == 'Pa':
                k = 'stress'
        if k not in ['stress', 'strain']:
            raise ValueError('invalid von Mises kind {kind}')

        return {
                'data':  (mechanics.equivalent_strain_Mises if k=='strain' else \
                          mechanics.equivalent_stress_Mises)(T_sym['data']),
                'label': f"{T_sym['label']}_vM",
                'meta':  {
                          'Unit':        T_sym['meta']['Unit'],
                          'Description': f"Mises equivalent {k} of {T_sym['label']} ({T_sym['meta']['Description']})",
                          'Creator':     'add_Mises'
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
                          'Unit':        x['meta']['Unit'],
                          'Description': f"{o}-norm of {t} {x['label']} ({x['meta']['Description']})",
                          'Creator':     'add_norm'
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
                          'Unit':        P['meta']['Unit'],
                          'Description': "2. Piola-Kirchhoff stress calculated "
                                         f"from {P['label']} ({P['meta']['Description']})"
                                         f" and {F['label']} ({F['meta']['Description']})",
                          'Creator':     'add_stress_second_Piola_Kirchhoff'
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
    #                       'Unit':        '1',
    #                       'Description': '{} coordinates of stereographic projection of pole (direction/plane) in crystal frame'\
    #                                      .format('Polar' if polar else 'Cartesian'),
    #                       'Creator':     'add_pole'
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
                          'Unit':        F['meta']['Unit'],
                          'Description': f"Rotational part of {F['label']} ({F['meta']['Description']})",
                          'Creator':     'add_rotation'
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
                          'Unit':        T['meta']['Unit'],
                          'Description': f"Spherical component of tensor {T['label']} ({T['meta']['Description']})",
                          'Creator':     'add_spherical'
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
                          'Unit':        F['meta']['Unit'],
                          'Description': f"Strain tensor of {F['label']} ({F['meta']['Description']})",
                          'Creator':     'add_strain'
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
                          'Unit':        F['meta']['Unit'],
                          'Description': '{} stretch tensor of {} ({})'.format('Left' if t.upper() == 'V' else 'Right',
                                                                               F['label'],F['meta']['Description']),
                          'Creator':     'add_stretch_tensor'
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
        num_threads = damask.environment.options['DAMASK_NUM_THREADS']
        pool = mp.Pool(int(num_threads) if num_threads is not None else None)
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
                        dataset.attrs['Overwritten'] = 'Yes' if h5py3 else \
                                                       'Yes'.encode()
                    else:
                        dataset = f[result[0]].create_dataset(result[1]['label'],data=result[1]['data'])

                    now = datetime.datetime.now().astimezone()
                    dataset.attrs['Created'] = now.strftime('%Y-%m-%d %H:%M:%S%z') if h5py3 else \
                                               now.strftime('%Y-%m-%d %H:%M:%S%z').encode()

                    for l,v in result[1]['meta'].items():
                        dataset.attrs[l]=v if h5py3 else v.encode()
                    creator = dataset.attrs['Creator'] if h5py3 else \
                              dataset.attrs['Creator'].decode()
                    dataset.attrs['Creator'] = f"damask.Result.{creator} v{damask.version}" if h5py3 else \
                                               f"damask.Result.{creator} v{damask.version}".encode()

                except (OSError,RuntimeError) as err:
                    print(f'Could not add dataset: {err}.')
            lock.release()

        pool.close()
        pool.join()


    def save_XDMF(self):
        """
        Write XDMF file to directly visualize data in DADF5 file.

        This works only for scalar, 3-vector and 3x3-tensor data.
        Selection is not taken into account.
        """
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
                                unit = f[name].attrs['Unit'] if h5py3 else f[name].attrs['Unit'].decode()

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

        with open(self.fname.with_suffix('.xdmf').name,'w') as f:
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

        N_digits = int(np.floor(np.log10(max(1,int(self.increments[-1][3:])))))+1

        for inc in util.show_progress(self.iterate('increments'),len(self.selection['increments'])):

            picked_backup_ho = self.selection['homogenizations'].copy()
            self.pick('homogenizations',False)
            for label in (labels if isinstance(labels,list) else [labels]):
                for o in self.iterate('out_type_ph'):
                    for c in range(self.N_constituents):
                        prefix = '' if self.N_constituents == 1 else f'constituent{c}/'
                        if o != 'mechanics':
                            for _ in self.iterate('phases'):
                                path = self.get_dataset_location(label)
                                if len(path) == 0:
                                    continue
                                array = self.read_dataset(path,c)
                                v.add(array,prefix+path[0].split('/',1)[1]+f' / {self._get_attribute(path[0],"Unit")}')
                        else:
                            paths = self.get_dataset_location(label)
                            if len(paths) == 0:
                                continue
                            array = self.read_dataset(paths,c)
                            ph_name = re.compile(r'(?<=(phase\/))(.*?)(?=(mechanics))')              # identify  phase name
                            dset_name = prefix+re.sub(ph_name,r'',paths[0].split('/',1)[1])          # remove phase name
                            v.add(array,dset_name+f' / {self._get_attribute(paths[0],"Unit")}')
            self.pick('homogenizations',picked_backup_ho)

            picked_backup_ph = self.selection['phases'].copy()
            self.pick('phases',False)
            for label in (labels if isinstance(labels,list) else [labels]):
                for _ in self.iterate('out_type_ho'):
                    paths = self.get_dataset_location(label)
                    if len(paths) == 0:
                        continue
                    array = self.read_dataset(paths)
                    v.add(array,paths[0].split('/',1)[1]+f' / {self._get_attribute(paths[0],"Unit")}')
            self.pick('phases',picked_backup_ph)

            u = self.read_dataset(self.get_dataset_location('u_n' if mode.lower() == 'cell' else 'u_p'))
            v.add(u,'u')

            v.save(f'{self.fname.stem}_inc{inc[3:].zfill(N_digits)}')
