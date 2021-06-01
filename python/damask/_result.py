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

import damask
from . import VTK
from . import Orientation
from . import grid_filters
from . import mechanics
from . import tensor
from . import util

h5py3 = h5py.__version__[0] == '3'

chunk_size = 1024**2//8                                                                             # for compression in HDF5


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
    Add data to and export data from a DADF5 file.

    A DADF5 (DAMASK HDF5) file contains DAMASK results.
    Its group/folder structure reflects the layout in material.yaml.

    This class provides a customizable view on the DADF5 file.
    Upon initialization, all attributes are visible.
    Derived quantities are added to the file and existing data is
    exported based on the current view.

    Examples
    --------
    Open 'my_file.hdf5', which is assumed to contain deformation gradient 'F'
    and first Piola-Kirchhoff stress 'P', add the Mises equivalent of the
    Cauchy stress, and export it to VTK (file) and numpy.ndarray (memory).

    >>> import damask
    >>> r = damask.Result('my_file.hdf5')
    >>> r.add_Cauchy()
    >>> r.add_equivalent_Mises('sigma')
    >>> r.save_VTK()
    >>> r_last = r.view('increments',-1)
    >>> sigma_vM_last = r_last.get('sigma_vM')

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

            if self.version_major != 0 or not 12 <= self.version_minor <= 13:
                raise TypeError(f'Unsupported DADF5 version {self.version_major}.{self.version_minor}')

            self.structured = 'cells' in f['geometry'].attrs.keys()

            if self.structured:
                self.cells  = f['geometry'].attrs['cells']
                self.size   = f['geometry'].attrs['size']
                self.origin = f['geometry'].attrs['origin']
            else:
                self.add_curl = self.add_divergence = self.add_gradient = None

            r=re.compile('increment_[0-9]+')
            self.increments = sorted([i for i in f.keys() if r.match(i)],key=util.natural_sort)
            self.times      = [round(f[i].attrs['t/s'],12) for i in self.increments]
            if len(self.increments) == 0:
                raise ValueError('incomplete DADF5 file')

            self.N_materialpoints, self.N_constituents = np.shape(f['cell_to/phase'])

            self.homogenization  = f['cell_to/homogenization']['label'].astype('str')
            self.homogenizations = sorted(np.unique(self.homogenization),key=util.natural_sort)
            self.phase           = f['cell_to/phase']['label'].astype('str')
            self.phases          = sorted(np.unique(self.phase),key=util.natural_sort)

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

        Returns
        -------
        view : damask.Result
            Modified or new view on the DADF5 file.

        """
        # allow True/False and string arguments
        if  datasets is True:
            datasets = '*'
        elif datasets is False or datasets is None:
            datasets = []
        choice = list(datasets).copy() if hasattr(datasets,'__iter__') and not isinstance(datasets,str) else \
                [datasets]

        what_ = what if what.endswith('s') else what+'s'

        if   what_ == 'increments':
            choice = [c if isinstance(c,str) and c.startswith('increment_') else
                      self.increments[c] if isinstance(c,int) and c<0 else
                      f'increment_{c}' for c in choice]
        elif what_ == 'times':
            what_ = 'increments'
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

        valid = _match(choice,getattr(self,what_))
        existing = set(self.visible[what_])

        dup = self.copy()
        if   action == 'set':
            dup.visible[what_] = sorted(set(valid), key=util.natural_sort)
        elif action == 'add':
            add = existing.union(valid)
            dup.visible[what_] = sorted(add, key=util.natural_sort)
        elif action == 'del':
            diff = existing.difference(valid)
            dup.visible[what_] = sorted(diff, key=util.natural_sort)

        return dup


    def modification_enable(self):
        """
        Allow modification of existing data.

        Returns
        -------
        modified_view : damask.Result
            View without write-protection of existing data.

        """
        print(util.warn('Warning: Modification of existing datasets allowed!'))
        dup = self.copy()
        dup._allow_modification = True
        return dup

    def modification_disable(self):
        """
        Prevent modification of existing data (default case).

        Returns
        -------
        modified_view : damask.Result
            View with write-protection of existing data.

        """
        dup = self.copy()
        dup._allow_modification = False
        return dup


    def increments_in_range(self,start,end):
        """
        Get all increments within a given range.

        Parameters
        ----------
        start : int or str
            Start increment.
        end : int or str
            End increment.

        Returns
        -------
        increments : list of ints
            Increment number of all increments within the given bounds.

        """
        selected = []
        for i,inc in enumerate([int(i[10:]) for i in self.increments]):
            s,e = map(lambda x: int(x[10:] if isinstance(x,str) and x.startswith('inc') else x), (start,end))
            if s <= inc <= e:
                selected.append(self.increments[i])
        return selected


    def times_in_range(self,start,end):
        """
        Get all increments within a given time range.

        Parameters
        ----------
        start : float
            Time of start increment.
        end : float
            Time of end increment.

        Returns
        -------
        times : list of float
            Simulation time of all increments within the given bounds.

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

        Returns
        -------
        view : damask.Result
            View with only the selected attributes being visible.

        Examples
        --------
        Get a view that shows only results from the initial configuration:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r_first = r.view('increment',0)

        Get a view that shows all results between simulation times of 10 to 40:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r_t10to40 = r.view('times',r.times_in_range(10.0,40.0))

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

        Returns
        -------
        modified_view : damask.Result
            View with additional visible attributes.

        Examples
        --------
        Get a view that shows only results from first and last increment:

        >>> import damask
        >>> r_empty = damask.Result('my_file.hdf5').view('increments',False)
        >>> r_first = r_empty.view_more('increments',0)
        >>> r_first_and_last = r.first.view_more('increments',-1)

        """
        return self._manage_view('add',what,datasets)


    def view_less(self,what,datasets):
        """
        Remove from view.

        Parameters
        ----------
        what : {'increments', 'times', 'phases', 'homogenizations', 'fields'}
            Attribute to change.
        datasets : (list of) int (for increments), (list of) float (for times), (list of) str, or bool
            Name of datasets; supports '?' and '*' wildcards.
            True is equivalent to '*', False is equivalent to [].

        Returns
        -------
        modified_view : damask.Result
            View with fewer visible attributes.

        Examples
        --------
        Get a view that omits the undeformed configuration:

        >>> import damask
        >>> r_all = damask.Result('my_file.hdf5')
        >>> r_deformed = r_all.view_less('increments',0)

        """
        return self._manage_view('del',what,datasets)


    def rename(self,name_src,name_dst):
        """
        Rename/move datasets (within the same group/folder).

        This operation is discouraged because the history of the
        data becomes untraceable and data integrity is not ensured.

        Parameters
        ----------
        name_src : str
            Name of the datasets to be renamed.
        name_dst : str
            New name of the datasets.

        Examples
        --------
        Rename datasets containing the deformation gradient from 'F' to 'def_grad':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r_unprotected = r.modification_enable()
        >>> r_unprotected.rename('F','def_grad')

        """
        if not self._allow_modification:
            raise PermissionError('Renaming datasets not permitted')

        with h5py.File(self.fname,'a') as f:
            for inc in self.visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            path_src = '/'.join([inc,ty,label,field,name_src])
                            path_dst = '/'.join([inc,ty,label,field,name_dst])
                            if path_src in f.keys():
                                f[path_dst] = f[path_src]
                                f[path_dst].attrs['renamed'] = f'original name: {name_src}' if h5py3 else \
                                                               f'original name: {name_src}'.encode()
                                del f[path_src]


    def remove(self,name):
        """
        Remove/delete datasets.

        This operation is discouraged because the history of the
        data becomes untraceable and data integrity is not ensured.

        Parameters
        ----------
        name : str
            Name of the datasets to be deleted.

        Examples
        --------
        Delete the deformation gradient 'F':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r_unprotected = r.modification_enable()
        >>> r_unprotected.remove('F')

        """
        if not self._allow_modification:
            raise PermissionError('Removing datasets not permitted')

        with h5py.File(self.fname,'a') as f:
            for inc in self.visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            path = '/'.join([inc,ty,label,field,name])
                            if path in f.keys(): del f[path]


    def list_data(self):
        """Return information on all active datasets in the file."""
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
                                unit = f' / {dataset.attrs["unit"]}' if h5py3 else \
                                       f' / {dataset.attrs["unit"].decode()}'
                                description = dataset.attrs['description'] if h5py3 else \
                                              dataset.attrs['description'].decode()
                                msg = '        '.join([msg,f'{d}{unit}: {description}\n'])

        return msg


    def enable_user_function(self,func):
        globals()[func.__name__]=func
        print(f'Function {func.__name__} enabled in add_calculation.')


    @property
    def coordinates0_point(self):
        """Initial/undeformed cell center coordinates."""
        if self.structured:
            return grid_filters.coordinates0_point(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_p'][()]

    @property
    def coordinates0_node(self):
        """Initial/undeformed nodal coordinates."""
        if self.structured:
            return grid_filters.coordinates0_node(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_n'][()]

    @property
    def geometry0(self):
        """Initial/undeformed geometry."""
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
            Name of scalar, vector, or tensor dataset to take absolute value of.

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
    def add_calculation(self,formula,name,unit='n/a',description=None):
        """
        Add result of a general formula.

        Parameters
        ----------
        formula : str
            Formula to calculate resulting dataset.
            Existing datasets are referenced by '#TheirName#'.
        name : str
            Name of resulting dataset.
        unit : str, optional
            Physical unit of the result.
        description : str, optional
            Human-readable description of the result.

        Examples
        --------
        Add total dislocation density, i.e. the sum of mobile dislocation
        density 'rho_mob' and dislocation dipole density 'rho_dip' over
        all slip systems:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_calculation('np.sum(#rho_mob#,axis=1)','rho_mob_total',
        ...                    '1/m²','total mobile dislocation density')
        >>> r.add_calculation('np.sum(#rho_dip#,axis=1)','rho_dip_total',
        ...                    '1/m²','total dislocation dipole density')
        >>> r.add_calculation('#rho_dip_total#+#rho_mob_total','rho_total',
        ...                    '1/m²','total dislocation density')

        Add Mises equivalent of the Cauchy stress without storage of
        intermediate results. Define a user function for better readability:

        >>> import damask
        >>> def equivalent_stress(F,P):
        ...     sigma = damask.mechanics.stress_Cauchy(F=F,P=P)
        ...     return damask.mechanics.equivalent_stress_Mises(sigma)
        >>> r = damask.Result('my_file.hdf5')
        >>> r.enable_user_function(equivalent_stress)
        >>> r.add_calculation('equivalent_stress(#F#,#P#)','sigma_vM','Pa',
        ...                   'Mises equivalent of the Cauchy stress')

        """
        dataset_mapping  = {d:d for d in set(re.findall(r'#(.*?)#',formula))}                       # datasets used in the formula
        args             = {'formula':formula,'label':name,'unit':unit,'description':description}
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
            Name of the dataset containing the first Piola-Kirchhoff stress. Defaults to 'P'.
        F : str, optional
            Name of the dataset containing the deformation gradient. Defaults to 'F'.

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
            Name of tensor dataset.

        Examples
        --------
        Add the determinant of plastic deformation gradient 'F_p':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_determinant('F_p')

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
            Name of tensor dataset.

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
            Name of symmetric tensor dataset.
        eigenvalue : {'max', 'mid', 'min'}
            Eigenvalue. Defaults to 'max'.

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
            Name of symmetric tensor dataset.
        eigenvalue : {'max', 'mid', 'min'}
            Eigenvalue to which the eigenvector corresponds.
            Defaults to 'max'.

        """
        self._add_generic_pointwise(self._add_eigenvector,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


    @staticmethod
    def _add_IPF_color(l,q):
        m = util.scale_to_coprime(np.array(l))
        lattice =  q['meta']['lattice']
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
            Name of the dataset containing the crystallographic orientation as quaternions.
            Defaults to 'O'.

        Examples
        --------
        Add the IPF color along [0,1,1] for orientation 'O':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_IPF_color(np.array([0,1,1]))

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
            Name of symmetric tensor dataset.

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
            Name of symmetric tensorial stress or strain dataset.
        kind : {'stress', 'strain', None}, optional
            Kind of the von Mises equivalent. Defaults to None, in which case
            it is selected based on the unit of the dataset ('1' -> strain, 'Pa' -> stress).

        Examples
        --------
        Add the Mises equivalent of the Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_equivalent_Mises('sigma')

        Add the Mises equivalent of the spatial logarithmic strain 'epsilon_V^0.0(F)':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_equivalent_Mises('epsilon_V^0.0(F)')

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
            raise ValueError(f'invalid norm order {ord}')

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
            Name of vector or tensor dataset.
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
            Name of first Piola-Kirchhoff stress dataset. Defaults to 'P'.
        F : str, optional
            Name of deformation gradient dataset. Defaults to 'F'.

        Notes
        -----
        The definition of the second Piola-Kirchhoff stress (S = [F^-1 P]_sym)
        follows the standard definition in nonlinear continuum mechanics.
        As such, no intermediate configuration, for instance that reached by F_p,
        is taken into account.

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
    #         Name of the dataset containing the crystallographic orientation as quaternions.
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
        F : str
            Name of deformation gradient dataset.

        Examples
        --------
        Add the rotational part of deformation gradient 'F':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_rotation('F')

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
            Name of tensor dataset.

        Examples
        --------
        Add the hydrostatic part of the Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_spherical('sigma')

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
            Name of deformation gradient dataset. Defaults to 'F'.
        t : {'V', 'U'}, optional
            Type of the polar decomposition, 'V' for left stretch tensor and 'U' for right stretch tensor.
            Defaults to 'V'.
        m : float, optional
            Order of the strain calculation. Defaults to 0.0.

        Examples
        --------
        Add the Biot strain based on the deformation gradient 'F':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.strain(t='U',m=0.5)

        Add the plastic Euler-Almansi strain based on the
        plastic deformation gradient 'F_p':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.strain('F_p','V',-1)

        """
        self._add_generic_pointwise(self._add_strain,{'F':F},{'t':t,'m':m})


    @staticmethod
    def _add_stretch_tensor(F,t):
        return {
                'data':  (mechanics.stretch_left if t.upper() == 'V' else mechanics.stretch_right)(F['data']),
                'label': f"{t}({F['label']})",
                'meta':  {
                          'unit':        F['meta']['unit'],
                          'description': f"{'left' if t.upper() == 'V' else 'right'} stretch tensor "\
                                        +f"of {F['label']} ({F['meta']['description']})",           # noqa
                          'creator':     'add_stretch_tensor'
                          }
                 }
    def add_stretch_tensor(self,F='F',t='V'):
        """
        Add stretch tensor of a deformation gradient.

        Parameters
        ----------
        F : str, optional
            Name of deformation gradient dataset. Defaults to 'F'.
        t : {'V', 'U'}, optional
            Type of the polar decomposition, 'V' for left stretch tensor and 'U' for right stretch tensor.
            Defaults to 'V'.

        """
        self._add_generic_pointwise(self._add_stretch_tensor,{'F':F},{'t':t})


    @staticmethod
    def _add_curl(f,size):
        return {
                'data':  grid_filters.curl(size,f['data']),
                'label': f"curl({f['label']})",
                'meta':  {
                          'unit':        f['meta']['unit']+'/m',
                          'description': f"curl of {f['label']} ({f['meta']['description']})",
                          'creator':     'add_curl'
                          }
                 }
    def add_curl(self,f):
        """
        Add curl of a field.

        Parameters
        ----------
        f : str
            Name of vector or tensor field dataset.

        Notes
        -----
        This function is only available for structured grids,
        i.e. results from the grid solver.

        """
        self._add_generic_grid(self._add_curl,{'f':f},{'size':self.size})


    @staticmethod
    def _add_divergence(f,size):
        return {
                'data':  grid_filters.divergence(size,f['data']),
                'label': f"divergence({f['label']})",
                'meta':  {
                          'unit':        f['meta']['unit']+'/m',
                          'description': f"divergence of {f['label']} ({f['meta']['description']})",
                          'creator':     'add_divergence'
                          }
                 }
    def add_divergence(self,f):
        """
        Add divergence of a field.

        Parameters
        ----------
        f : str
            Name of vector or tensor field dataset.

        Notes
        -----
        This function is only available for structured grids,
        i.e. results from the grid solver.

        """
        self._add_generic_grid(self._add_divergence,{'f':f},{'size':self.size})


    @staticmethod
    def _add_gradient(f,size):
        return {
                'data':  grid_filters.gradient(size,f['data'] if len(f['data'].shape) == 4 else \
                                                    f['data'].reshape(f['data'].shape+(1,))),
                'label': f"gradient({f['label']})",
                'meta':  {
                          'unit':        f['meta']['unit']+'/m',
                          'description': f"gradient of {f['label']} ({f['meta']['description']})",
                          'creator':     'add_gradient'
                          }
                 }
    def add_gradient(self,f):
        """
        Add gradient of a field.

        Parameters
        ----------
        f : str
            Name of scalar or vector field dataset.

        Notes
        -----
        This function is only available for structured grids,
        i.e. results from the grid solver.

        """
        self._add_generic_grid(self._add_gradient,{'f':f},{'size':self.size})


    def _add_generic_grid(self,func,datasets,args={},constituents=None):
        """
        General function to add data on a regular grid.

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
        if len(datasets) != 1 or self.N_constituents !=1:
            raise NotImplementedError

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        with h5py.File(self.fname, 'a') as f:
            for increment in self.place(datasets.values(),False).items():
                for ty in increment[1].items():
                    for field in ty[1].items():
                        d = list(field[1].values())[0]
                        if np.any(d.mask): continue
                        dataset = {'f':{'data':np.reshape(d.data,tuple(self.cells)+d.data.shape[1:]),
                                        'label':list(datasets.values())[0],
                                        'meta':d.data.dtype.metadata}}
                        r = func(**dataset,**args)
                        result = r['data'].reshape((-1,)+r['data'].shape[3:])
                        for x in self.visible[ty[0]+'s']:
                            if ty[0] == 'phase':
                                result1 = result[at_cell_ph[0][x]]
                            if ty[0] == 'homogenization':
                                result1 = result[at_cell_ho[x]]

                            path = '/'.join(['/',increment[0],ty[0],x,field[0]])
                            dataset = f[path].create_dataset(r['label'],data=result1)

                            now = datetime.datetime.now().astimezone()
                            dataset.attrs['created'] = now.strftime('%Y-%m-%d %H:%M:%S%z') if h5py3 else \
                                                       now.strftime('%Y-%m-%d %H:%M:%S%z').encode()

                            for l,v in r['meta'].items():
                                dataset.attrs[l.lower()]=v if h5py3 else v.encode()
                            creator = dataset.attrs['creator'] if h5py3 else \
                                      dataset.attrs['creator'].decode()
                            dataset.attrs['creator'] = f'damask.Result.{creator} v{damask.version}' if h5py3 else \
                                                       f'damask.Result.{creator} v{damask.version}'.encode()


    def _job_pointwise(self,group,func,datasets,args,lock):
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
            return [None,None]

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
        pool = mp.Pool(int(os.environ.get('OMP_NUM_THREADS',4)))
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

        default_arg = partial(self._job_pointwise,func=func,datasets=datasets,args=args,lock=lock)

        for group,result in util.show_progress(pool.imap_unordered(default_arg,groups),len(groups)):
            if not result:
                continue
            lock.acquire()
            with h5py.File(self.fname, 'a') as f:
                try:
                    if self._allow_modification and '/'.join([group,result['label']]) in f:
                        dataset = f['/'.join([group,result['label']])]
                        dataset[...] = result['data']
                        dataset.attrs['overwritten'] = True
                    else:
                        if result['data'].size >= chunk_size*2:
                            shape  = result['data'].shape
                            chunks = (chunk_size//np.prod(shape[1:]),)+shape[1:]
                            dataset = f[group].create_dataset(result['label'],data=result['data'],
                                                              maxshape=shape, chunks=chunks,
                                                              compression='gzip', compression_opts=6,
                                                              shuffle=True,fletcher32=True)
                        else:
                            dataset = f[group].create_dataset(result['label'],data=result['data'])

                    now = datetime.datetime.now().astimezone()
                    dataset.attrs['created'] = now.strftime('%Y-%m-%d %H:%M:%S%z') if h5py3 else \
                                               now.strftime('%Y-%m-%d %H:%M:%S%z').encode()

                    for l,v in result['meta'].items():
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

        The XDMF format is only supported for structured grids
        with single phase and single constituent.
        For other cases use `save_VTK`.

        Parameters
        ----------
        output : (list of) str
            Names of the datasets included in the XDMF file.
            Defaults to '*', in which case all datasets are considered.

        """
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

                                unit = f[name].attrs['unit'] if h5py3 else \
                                       f[name].attrs['unit'].decode()

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
        """Mappings to place data spatially."""
        with h5py.File(self.fname,'r') as f:

            at_cell_ph = []
            in_data_ph = []
            for c in range(self.N_constituents):
                at_cell_ph.append({label: np.where(self.phase[:,c] == label)[0] \
                                          for label in self.visible['phases']})
                in_data_ph.append({label: f['/'.join(['cell_to','phase'])]['entry'][at_cell_ph[c][label]][:,c] \
                                          for label in self.visible['phases']})

            at_cell_ho = {label: np.where(self.homogenization[:] == label)[0] \
                                 for label in self.visible['homogenizations']}
            in_data_ho = {label: f['/'.join(['cell_to','homogenization'])]['entry'][at_cell_ho[label]] \
                                 for label in self.visible['homogenizations']}

        return at_cell_ph,in_data_ph,at_cell_ho,in_data_ho


    def save_VTK(self,output='*',mode='cell',constituents=None,fill_float=np.nan,fill_int=0,parallel=True):
        """
        Export to VTK cell/point data.

        One VTK file per visible increment is created.
        For cell data, the VTK format is a rectilinear grid (.vtr) for
        grid-based simulations and an unstructured grid (.vtu) for
        mesh-baed simulations. For point data, the VTK format is poly
        data (.vtp).

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets to export to the VTK file.
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
        else:
            raise ValueError(f'invalid mode {mode}')

        v.set_comments(util.execution_stamp('Result','save_VTK'))

        N_digits = int(np.floor(np.log10(max(1,int(self.increments[-1][10:])))))+1

        constituents_ = constituents if isinstance(constituents,Iterable) else \
                      (range(self.N_constituents) if constituents is None else [constituents])

        suffixes = [''] if self.N_constituents == 1 or isinstance(constituents,int) else \
                   [f'#{c}' for c in constituents_]

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        with h5py.File(self.fname,'r') as f:
            if self.version_minor >= 13:
                creator = f.attrs['creator'] if h5py3 else f.attrs['creator'].decode()
                created = f.attrs['created'] if h5py3 else f.attrs['created'].decode()
                v.add_comments(f'{creator} ({created})')

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

                v.save(f'{self.fname.stem}_inc{inc[10:].zfill(N_digits)}',parallel=parallel)


    def get(self,output='*',flatten=True,prune=True):
        """
        Collect data per phase/homogenization reflecting the group/folder structure in the DADF5 file.

        Parameters
        ----------
        output : (list of) str
            Names of the datasets to read.
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
        `place` is equivalent to `get` if only one phase/homogenization
        and one constituent is present.

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets to read.
            Defaults to '*', in which case all datasets are placed.
        flatten : bool
            Remove singular levels of the folder hierarchy.
            This might be beneficial in case of single increment or field.
            Defaults to True.
        prune : bool
            Remove branches with no data. Defaults to True.
        constituents : (list of) int, optional
            Constituents to consider.
            Defaults to None, in which case all constituents are considered.
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
