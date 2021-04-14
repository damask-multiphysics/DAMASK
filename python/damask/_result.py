import multiprocessing as mp
import re
import fnmatch
import os
import copy
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


def _read(dataset):
    """Read a dataset and its metadata into a numpy.ndarray."""
    metadata = {k:(v.decode() if not h5py3 and type(v) is bytes else v) for k,v in dataset.attrs.items()}
    dtype = np.dtype(dataset.dtype,metadata=metadata)
    return np.array(dataset,dtype=dtype)

def _match(requested,existing):
    """Find matches among two sets of labels."""
    def flatten_list(list_of_lists):
        return [e for e_ in list_of_lists for e in e_]

    if requested is True:
        requested = '*'
    elif requested is False or requested is None:
        requested = []

    requested_ = requested if hasattr(requested,'__iter__') and not isinstance(requested,str) else \
                [requested]

    return sorted(set(flatten_list([fnmatch.filter(existing,r) for r in requested_])),
                  key=util.natural_sort)

def _empty_like(dataset,N_materialpoints,fill_float,fill_int):
    """Create empty numpy.ma.MaskedArray."""
    return ma.array(np.empty((N_materialpoints,)+dataset.shape[1:],dataset.dtype),
                    fill_value = fill_float if dataset.dtype in np.sctypes['float'] else fill_int,
                    mask = True)

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
            self.increments = sorted([i for i in f.keys() if r.match(i)],key=util.natural_sort)
            self.times      = [round(f[i].attrs['time/s' if self.version_minor < 12 else
                                                't/s'],12) for i in self.increments]

            grp = 'mapping' if self.version_minor < 12 else 'cell_to'

            self.N_materialpoints, self.N_constituents = np.shape(f[f'{grp}/phase'])

            self.homogenizations  = [m.decode() for m in np.unique(f[f'{grp}/homogenization']
                                                                    ['Name' if self.version_minor < 12 else 'label'])]
            self.homogenizations  = sorted(self.homogenizations,key=util.natural_sort)
            self.phases           = [c.decode() for c in np.unique(f[f'{grp}/phase']
                                                                    ['Name' if self.version_minor < 12 else 'label'])]
            self.phases           = sorted(self.phases,key=util.natural_sort)

            self.fields = []
            for c in self.phases:
                self.fields += f['/'.join([self.increments[0],'phase',c])].keys()
            for m in self.homogenizations:
                self.fields += f['/'.join([self.increments[0],'homogenization',m])].keys()
            self.fields = sorted(set(self.fields),key=util.natural_sort)                            # make unique

        self.visible = {'increments':      self.increments,
                        'phases':          self.phases,
                        'homogenizations': self.homogenizations,
                        'fields':          self.fields,
                       }

        self.fname = Path(fname).absolute()

        self._allow_modification = False


    def __copy__(self):
        """Create deep copy."""
        return copy.deepcopy(self)

    copy = __copy__


    def __repr__(self):
        """Show summary of file content."""
        visible_increments = self.visible['increments']

        first = self.view('increments',visible_increments[0:1]).list_data()

        last  = '' if len(visible_increments) < 2 else \
                self.view('increments',visible_increments[-1:]).list_data()

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
        datasets : (list of) int (for increments), (list of) float (for times), (list of) str, or bool
            Name of datasets; supports '?' and '*' wildcards.
            True is equivalent to '*', False is equivalent to [].

        """
        # allow True/False and string arguments
        if  datasets is True:
            datasets = '*'
        elif datasets is False or datasets is None:
            datasets = []
        choice = list(datasets).copy() if hasattr(datasets,'__iter__') and not isinstance(datasets,str) else \
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

        valid = _match(choice,getattr(self,what))
        existing = set(self.visible[what])

        dup = self.copy()
        if   action == 'set':
            dup.visible[what] = sorted(set(valid), key=util.natural_sort)
        elif action == 'add':
            add = existing.union(valid)
            dup.visible[what] = sorted(add, key=util.natural_sort)
        elif action == 'del':
            diff = existing.difference(valid)
            dup.visible[what] = sorted(diff, key=util.natural_sort)

        return dup


    def allow_modification(self):
        """Allow to overwrite existing data."""
        print(util.warn('Warning: Modification of existing datasets allowed!'))
        dup = self.copy()
        dup._allow_modification = True
        return dup

    def disallow_modification(self):
        """Disallow to overwrite existing data (default case)."""
        dup = self.copy()
        dup._allow_modification = False
        return dup


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


    def view(self,what,datasets):
        """
        Set view.

        Parameters
        ----------
        what : {'increments', 'times', 'phases', 'homogenizations', 'fields'}
            Attribute to change.
        datasets : (list of) int (for increments), (list of) float (for times), (list of) str, or bool
            Name of datasets; supports '?' and '*' wildcards.
            True is equivalent to '*', False is equivalent to [].

        """
        return self._manage_view('set',what,datasets)


    def view_more(self,what,datasets):
        """
        Add to view.

        Parameters
        ----------
        what : {'increments', 'times', 'phases', 'homogenizations', 'fields'}
            Attribute to change.
        datasets : (list of) int (for increments), (list of) float (for times), (list of) str, or bool
            Name of datasets; supports '?' and '*' wildcards.
            True is equivalent to '*', False is equivalent to [].

        """
        return self._manage_view('add',what,datasets)


    def view_less(self,what,datasets):
        """
        Delete from view.

        Parameters
        ----------
        what : {'increments', 'times', 'phases', 'homogenizations', 'fields'}
            Attribute to change.
        datasets : (list of) int (for increments), (list of) float (for times), (list of) str, or bool
            Name of datasets; supports '?' and '*' wildcards.
            True is equivalent to '*', False is equivalent to [].

        """
        return self._manage_view('del',what,datasets)


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
        if not self._allow_modification:
            raise PermissionError('Rename operation not permitted')

        with h5py.File(self.fname,'a') as f:
            for inc in self.visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            path_old = '/'.join([inc,ty,label,field,name_old])
                            path_new = '/'.join([inc,ty,label,field,name_new])
                            if path_old in f.keys():
                                f[path_new] = f[path_old]
                                f[path_new].attrs['renamed'] = f'original name: {name_old}' if h5py3 else \
                                                               f'original name: {name_old}'.encode()
                                del f[path_old]


    def list_data(self):
        """Return information on all active datasets in the file."""
        # compatibility hack
        de = 'Description' if self.version_minor < 12 else 'description'
        un = 'Unit'        if self.version_minor < 12 else 'unit'
        msg = ''
        with h5py.File(self.fname,'r') as f:
            for inc in self.visible['increments']:
                msg = ''.join([msg,f'\n{inc} ({self.times[self.increments.index(inc)]}s)\n'])
                for ty in ['phase','homogenization']:
                    msg = '  '.join([msg,f'{ty}\n'])
                    for label in self.visible[ty+'s']:
                        msg = '    '.join([msg,f'{label}\n'])
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            msg = '      '.join([msg,f'{field}\n'])
                            for d in f['/'.join([inc,ty,label,field])].keys():
                                dataset = f['/'.join([inc,ty,label,field,d])]
                                unit = f' / {dataset.attrs[un]}' if h5py3 else \
                                       f' / {dataset.attrs[un].decode()}'
                                description = dataset.attrs[de] if h5py3 else \
                                              dataset.attrs[de].decode()
                                msg = '        '.join([msg,f'{d}{unit}: {description}\n'])

        return msg


    def enable_user_function(self,func):
        globals()[func.__name__]=func
        print(f'Function {func.__name__} enabled in add_calculation.')


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

    @property
    def geometry0(self):
        if self.structured:
            return VTK.from_rectilinear_grid(self.cells,self.size,self.origin)
        else:
            with h5py.File(self.fname,'r') as f:
                return VTK.from_unstructured_grid(f['/geometry/x_n'][()],
                                                  f['/geometry/T_c'][()]-1,
                                                  f['/geometry/T_c'].attrs['VTK_TYPE'] if h5py3 else \
                                                  f['/geometry/T_c'].attrs['VTK_TYPE'].decode())


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
            Formula to calculate resulting dataset. Existing datasets are referenced by '#TheirLabel#'.
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
            Label of the dataset containing the first Piola-Kirchhoff stress. Defaults to 'P'.
        F : str, optional
            Label of the dataset containing the deformation gradient. Defaults to 'F'.

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
            Eigenvalue. Select from 'max', 'mid', 'min'. Defaults to 'max'.

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
            Eigenvalue to which the eigenvector corresponds.
            Select from 'max', 'mid', 'min'. Defaults to 'max'.

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
                'label': 'IPFcolor_({} {} {})'.format(*m),
                'meta' : {
                          'unit':        '8-bit RGB',
                          'lattice':     q['meta']['lattice'],
                          'description': 'Inverse Pole Figure (IPF) colors along sample direction ({} {} {})'.format(*m),
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
            raise ValueError(f'invalid von Mises kind {kind}')

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
            it is selected based on the unit of the dataset ('1' -> strain, 'Pa' -> stress).

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
        ord : {non-zero int, inf, -inf, 'fro', 'nuc'}, optional
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
            Label of first Piola-Kirchhoff stress dataset. Defaults to 'P'.
        F : str, optional
            Label of deformation gradient dataset. Defaults to 'F'.

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

        For details, see damask.mechanics.strain.

        Parameters
        ----------
        F : str, optional
            Label of deformation gradient dataset. Defaults to 'F'.
        t : {'V', 'U'}, optional
            Type of the polar decomposition, 'V' for left stretch tensor and 'U' for right stretch tensor.
            Defaults to 'V'.
        m : float, optional
            Order of the strain calculation. Defaults to 0.0.

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
            Label of deformation gradient dataset. Defaults to 'F'.
        t : {'V', 'U'}, optional
            Type of the polar decomposition, 'V' for left stretch tensor and 'U' for right stretch tensor.
            Defaults to 'V'.

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
            Details of the datasets to be used:
            {arg (name to which the data is passed in func): label (in HDF5 file)}.
        args : dictionary, optional
            Arguments parsed to func.

        """
        chunk_size = 1024**2//8
        pool = mp.Pool(int(os.environ.get('OMP_NUM_THREADS',1)))
        lock = mp.Manager().Lock()

        groups = []
        with h5py.File(self.fname,'r') as f:
            for inc in self.visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            group = '/'.join([inc,ty,label,field])
                            if set(datasets.values()).issubset(f[group].keys()): groups.append(group)

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
                    dataset.attrs['creator'] = f'damask.Result.{creator} v{damask.version}' if h5py3 else \
                                               f'damask.Result.{creator} v{damask.version}'.encode()

                except (OSError,RuntimeError) as err:
                    print(f'Could not add dataset: {err}.')
            lock.release()

        pool.close()
        pool.join()


    def save_XDMF(self,output='*'):
        """
        Write XDMF file to directly visualize data in DADF5 file.

        Parameters
        ----------
        output : (list of) str
            Labels of the datasets to read.
            Defaults to '*', in which case all datasets are considered.

        """
        u = 'Unit' if self.version_minor < 12 else 'unit'                                           # compatibility hack
        if self.N_constituents != 1 or len(self.phases) != 1 or not self.structured:
            raise TypeError('XDMF output requires structured grid with single phase and single constituent.')


        attribute_type_map = defaultdict(lambda:'Matrix', ( ((),'Scalar'), ((3,),'Vector'), ((3,3),'Tensor')) )

        def number_type_map(dtype):
            if dtype in np.sctypes['int']:   return 'Int'
            if dtype in np.sctypes['uint']:  return 'UInt'
            if dtype in np.sctypes['float']: return 'Float'


        xdmf = ET.Element('Xdmf')
        xdmf.attrib={'Version':  '2.0',
                     'xmlns:xi': 'http://www.w3.org/2001/XInclude'}

        domain = ET.SubElement(xdmf, 'Domain')

        collection = ET.SubElement(domain, 'Grid')
        collection.attrib={'GridType':       'Collection',
                           'CollectionType': 'Temporal',
                           'Name':           'Increments'}

        time = ET.SubElement(collection, 'Time')
        time.attrib={'TimeType': 'List'}

        time_data = ET.SubElement(time, 'DataItem')
        times = [self.times[self.increments.index(i)] for i in self.visible['increments']]
        time_data.attrib={'Format':     'XML',
                          'NumberType': 'Float',
                          'Dimensions': f'{len(times)}'}
        time_data.text = ' '.join(map(str,times))

        attributes = []
        data_items = []

        with h5py.File(self.fname,'r') as f:
            for inc in self.visible['increments']:

                grid = ET.SubElement(collection,'Grid')
                grid.attrib = {'GridType': 'Uniform',
                               'Name':      inc}

                topology = ET.SubElement(grid, 'Topology')
                topology.attrib = {'TopologyType': '3DCoRectMesh',
                                   'Dimensions':   '{} {} {}'.format(*(self.cells+1))}

                geometry = ET.SubElement(grid, 'Geometry')
                geometry.attrib = {'GeometryType':'Origin_DxDyDz'}

                origin = ET.SubElement(geometry, 'DataItem')
                origin.attrib = {'Format':     'XML',
                                 'NumberType': 'Float',
                                 'Dimensions': '3'}
                origin.text = "{} {} {}".format(*self.origin)

                delta = ET.SubElement(geometry, 'DataItem')
                delta.attrib = {'Format':     'XML',
                                'NumberType': 'Float',
                                'Dimensions': '3'}
                delta.text="{} {} {}".format(*(self.size/self.cells))

                attributes.append(ET.SubElement(grid, 'Attribute'))
                attributes[-1].attrib = {'Name':          'u / m',
                                         'Center':        'Node',
                                         'AttributeType': 'Vector'}
                data_items.append(ET.SubElement(attributes[-1], 'DataItem'))
                data_items[-1].attrib = {'Format':     'HDF',
                                         'Precision':  '8',
                                         'Dimensions': '{} {} {} 3'.format(*(self.cells+1))}
                data_items[-1].text = f'{os.path.split(self.fname)[1]}:/{inc}/geometry/u_n'

                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                name = '/'.join([inc,ty,label,field,out])
                                shape = f[name].shape[1:]
                                dtype = f[name].dtype

                                unit = f[name].attrs[u] if h5py3 else f[name].attrs[u].decode()

                                attributes.append(ET.SubElement(grid, 'Attribute'))
                                attributes[-1].attrib = {'Name':          '/'.join([ty,field,out])+f' / {unit}',
                                                         'Center':       'Cell',
                                                         'AttributeType': attribute_type_map[shape]}
                                data_items.append(ET.SubElement(attributes[-1], 'DataItem'))
                                data_items[-1].attrib = {'Format':     'HDF',
                                                         'NumberType': number_type_map(dtype),
                                                         'Precision':  f'{dtype.itemsize}',
                                                         'Dimensions': '{} {} {} {}'.format(*self.cells,1 if shape == () else
                                                                                                        np.prod(shape))}
                                data_items[-1].text = f'{os.path.split(self.fname)[1]}:{name}'

        with open(self.fname.with_suffix('.xdmf').name,'w',newline='\n') as f:
            f.write(xml.dom.minidom.parseString(ET.tostring(xdmf).decode()).toprettyxml())


    def _mappings(self):
        grp    = 'mapping' if self.version_minor < 12 else 'cell_to'                              # compatibility hack
        name   = 'Name'    if self.version_minor < 12 else 'label'                                # compatibility hack
        member = 'member'  if self.version_minor < 12 else 'entry'                                # compatibility hack

        with h5py.File(self.fname,'r') as f:

            at_cell_ph = []
            in_data_ph = []
            for c in range(self.N_constituents):
                at_cell_ph.append({label: np.where(f['/'.join([grp,'phase'])][:,c][name] == label.encode())[0] \
                                          for label in self.visible['phases']})
                in_data_ph.append({label: f['/'.join([grp,'phase'])][member][at_cell_ph[c][label]][:,c] \
                                          for label in self.visible['phases']})

            at_cell_ho = {label: np.where(f['/'.join([grp,'homogenization'])][:][name] == label.encode())[0] \
                                 for label in self.visible['homogenizations']}
            in_data_ho = {label: f['/'.join([grp,'homogenization'])][member][at_cell_ho[label]] \
                                 for label in self.visible['homogenizations']}

        return at_cell_ph,in_data_ph,at_cell_ho,in_data_ho


    def save_VTK(self,output='*',mode='cell',constituents=None,fill_float=np.nan,fill_int=0,parallel=True):
        """
        Export to VTK cell/point data.

        Parameters
        ----------
        output : (list of) str, optional
            Labels of the datasets to place.
            Defaults to '*', in which case all datasets are exported.
        mode : {'cell', 'point'}
            Export in cell format or point format.
            Defaults to 'cell'.
        constituents : (list of) int, optional
            Constituents to consider.
            Defaults to None, in which case all constituents are considered.
        fill_float : float
            Fill value for non-existent entries of floating point type.
            Defaults to NaN.
        fill_int : int
            Fill value for non-existent entries of integer type.
            Defaults to 0.
        parallel : bool
            Write out VTK files in parallel in a separate background process.
            Defaults to True.

        """
        if mode.lower()=='cell':
            v = self.geometry0
        elif mode.lower()=='point':
            v = VTK.from_poly_data(self.coordinates0_point)

        ln = 3 if self.version_minor < 12 else 10                                                   # compatibility hack
        N_digits = int(np.floor(np.log10(max(1,int(self.increments[-1][ln:])))))+1

        constituents_ = constituents if isinstance(constituents,Iterable) else \
                      (range(self.N_constituents) if constituents is None else [constituents])

        suffixes = [''] if self.N_constituents == 1 or isinstance(constituents,int) else \
                   [f'#{c}' for c in constituents_]

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        with h5py.File(self.fname,'r') as f:

            for inc in util.show_progress(self.visible['increments']):

                u = _read(f['/'.join([inc,'geometry','u_n' if mode.lower() == 'cell' else 'u_p'])])
                v.add(u,'u')

                for ty in ['phase','homogenization']:
                    for field in self.visible['fields']:
                        outs = {}
                        for label in self.visible[ty+'s']:
                            if field not in f['/'.join([inc,ty,label])].keys(): continue

                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                data = ma.array(_read(f['/'.join([inc,ty,label,field,out])]))

                                if ty == 'phase':
                                    if out+suffixes[0] not in outs.keys():
                                        for c,suffix in zip(constituents_,suffixes):
                                            outs[out+suffix] = \
                                                _empty_like(data,self.N_materialpoints,fill_float,fill_int)

                                    for c,suffix in zip(constituents_,suffixes):
                                        outs[out+suffix][at_cell_ph[c][label]] = data[in_data_ph[c][label]]

                                if ty == 'homogenization':
                                    if out not in outs.keys():
                                        outs[out] = _empty_like(data,self.N_materialpoints,fill_float,fill_int)

                                    outs[out][at_cell_ho[label]] = data[in_data_ho[label]]

                        for label,dataset in outs.items():
                            v.add(dataset,' / '.join(['/'.join([ty,field,label]),dataset.dtype.metadata['unit']]))

                v.save(f'{self.fname.stem}_inc{inc[ln:].zfill(N_digits)}',parallel=parallel)


    def get(self,output='*',flatten=True,prune=True):
        """
        Collect data per phase/homogenization reflecting the group/folder structure in the DADF5 file.

        Parameters
        ----------
        output : (list of) str
            Labels of the datasets to read.
            Defaults to '*', in which case all datasets are read.
        flatten : bool
            Remove singular levels of the folder hierarchy.
            This might be beneficial in case of single increment,
            phase/homogenization, or field. Defaults to True.
        prune : bool
            Remove branches with no data. Defaults to True.

        Returns
        -------
        data : dict of numpy.ndarray
            Datasets structured by phase/homogenization and according to selected view.

        """
        r = {}

        with h5py.File(self.fname,'r') as f:
            for inc in util.show_progress(self.visible['increments']):
                r[inc] = {'phase':{},'homogenization':{},'geometry':{}}

                for out in _match(output,f['/'.join([inc,'geometry'])].keys()):
                    r[inc]['geometry'][out] = _read(f['/'.join([inc,'geometry',out])])

                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        r[inc][ty][label] = {}
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            r[inc][ty][label][field] = {}
                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                r[inc][ty][label][field][out] = _read(f['/'.join([inc,ty,label,field,out])])

        if prune:   r = util.dict_prune(r)
        if flatten: r = util.dict_flatten(r)

        return None if (type(r) == dict and r == {}) else r


    def place(self,output='*',flatten=True,prune=True,constituents=None,fill_float=np.nan,fill_int=0):
        """
        Merge data into spatial order that is compatible with the damask.VTK geometry representation.

        The returned data structure reflects the group/folder structure
        in the DADF5 file.

        Multi-phase data is fused into a single output.
        `place` is equivalent to `read` if only one phase/homogenization
        and one constituent is present.

        Parameters
        ----------
        output : (list of) str, optional
            Labels of the datasets to place.
            Defaults to '*', in which case all datasets are placed.
        flatten : bool
            Remove singular levels of the folder hierarchy.
            This might be beneficial in case of single increment or field.
            Defaults to True.
        prune : bool
            Remove branches with no data. Defaults to True.
        constituents : (list of) int, optional
            Constituents to consider.
            Defaults to 'None', in which case all constituents are considered.
        fill_float : float
            Fill value for non-existent entries of floating point type.
            Defaults to NaN.
        fill_int : int
            Fill value for non-existent entries of integer type.
            Defaults to 0.

        Returns
        -------
        data : dict of numpy.ma.MaskedArray
            Datasets structured by spatial position and according to selected view.

        """
        r = {}

        constituents_ = constituents if isinstance(constituents,Iterable) else \
                      (range(self.N_constituents) if constituents is None else [constituents])

        suffixes = [''] if self.N_constituents == 1 or isinstance(constituents,int) else \
                   [f'#{c}' for c in constituents_]

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        with h5py.File(self.fname,'r') as f:

            for inc in util.show_progress(self.visible['increments']):
                r[inc] = {'phase':{},'homogenization':{},'geometry':{}}

                for out in _match(output,f['/'.join([inc,'geometry'])].keys()):
                    r[inc]['geometry'][out] = ma.array(_read(f['/'.join([inc,'geometry',out])]),fill_value = fill_float)

                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            if field not in r[inc][ty].keys():
                                r[inc][ty][field] = {}

                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                data = ma.array(_read(f['/'.join([inc,ty,label,field,out])]))

                                if ty == 'phase':
                                    if out+suffixes[0] not in r[inc][ty][field].keys():
                                        for c,suffix in zip(constituents_,suffixes):
                                            r[inc][ty][field][out+suffix] = \
                                                _empty_like(data,self.N_materialpoints,fill_float,fill_int)

                                    for c,suffix in zip(constituents_,suffixes):
                                        r[inc][ty][field][out+suffix][at_cell_ph[c][label]] = data[in_data_ph[c][label]]

                                if ty == 'homogenization':
                                    if out not in r[inc][ty][field].keys():
                                        r[inc][ty][field][out] = \
                                            _empty_like(data,self.N_materialpoints,fill_float,fill_int)

                                    r[inc][ty][field][out][at_cell_ho[label]] = data[in_data_ho[label]]

        if prune:   r = util.dict_prune(r)
        if flatten: r = util.dict_flatten(r)

        return None if (type(r) == dict and r == {}) else r
