import multiprocessing as mp
from multiprocessing.synchronize import Lock
import re
import fnmatch
import os
import copy
import datetime
import xml.etree.ElementTree as ET                                                                  # noqa
import xml.dom.minidom
from pathlib import Path
from functools import partial
from collections import defaultdict
from collections.abc import Iterable
from typing import Optional, Union, Callable, Any, Sequence, Literal, Dict, List, Tuple

import h5py
import numpy as np
from numpy import ma

import damask
from . import VTK
from . import Orientation
from . import grid_filters
from . import mechanics
from . import tensor
from . import util
from ._typehints import FloatSequence, IntSequence

h5py3 = h5py.__version__[0] == '3'

chunk_size = 1024**2//8                                                                             # for compression in HDF5
prefix_inc = 'increment_'

def _read(dataset: h5py._hl.dataset.Dataset) -> np.ndarray:
    """Read a dataset and its metadata into a numpy.ndarray."""
    metadata = {k:(v.decode() if not h5py3 and type(v) is bytes else v) for k,v in dataset.attrs.items()}
    dtype = np.dtype(dataset.dtype,metadata=metadata)                                               # type: ignore
    return np.array(dataset,dtype=dtype)

def _match(requested,
           existing: h5py._hl.base.KeysViewHDF5) -> List[Any]:
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

def _empty_like(dataset: np.ma.core.MaskedArray,
                N_materialpoints: int,
                fill_float: float,
                fill_int: int) -> np.ma.core.MaskedArray:
    """Create empty numpy.ma.MaskedArray."""
    return ma.array(np.empty((N_materialpoints,)+dataset.shape[1:],dataset.dtype),
                    fill_value = fill_float if dataset.dtype in np.sctypes['float'] else fill_int,
                    mask = True)

class Result:
    """
    Add data to and export data from a DADF5 (DAMASK HDF5) file.

    A DADF5 file contains DAMASK results.
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
    >>> r.add_stress_Cauchy()
    >>> r.add_equivalent_Mises('sigma')
    >>> r.export_VTK()
    >>> r_last = r.view(increments=-1)
    >>> sigma_vM_last = r_last.get('sigma_vM')

    """

    def __init__(self, fname: Union[str, Path]):
        """
        New result view bound to a DADF5 file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Name of the DADF5 file to be opened.

        """
        with h5py.File(fname,'r') as f:

            self.version_major = f.attrs['DADF5_version_major']
            self.version_minor = f.attrs['DADF5_version_minor']

            if (self.version_major != 0 or not 14 <= self.version_minor <= 14) and self.version_major != 1:
                raise TypeError(f'unsupported DADF5 version "{self.version_major}.{self.version_minor}"')

            self.structured = 'cells' in f['geometry'].attrs.keys()

            if self.structured:
                self.cells  = f['geometry'].attrs['cells']
                self.size   = f['geometry'].attrs['size']
                self.origin = f['geometry'].attrs['origin']
            else:
                self.add_curl = self.add_divergence = self.add_gradient = None                      # type: ignore

            r = re.compile(rf'{prefix_inc}([0-9]+)')
            self.increments = sorted([i for i in f.keys() if r.match(i)],key=util.natural_sort)
            self.times = np.around([f[i].attrs['t/s'] for i in self.increments],12)
            if len(self.increments) == 0:
                raise ValueError('incomplete DADF5 file')

            self.N_materialpoints, self.N_constituents = np.shape(f['cell_to/phase'])

            self.homogenization  = f['cell_to/homogenization']['label'].astype('str')
            self.homogenizations = sorted(np.unique(self.homogenization),key=util.natural_sort)
            self.phase           = f['cell_to/phase']['label'].astype('str')
            self.phases          = sorted(np.unique(self.phase),key=util.natural_sort)

            self.fields: List[str] = []
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

        self.fname = Path(fname).expanduser().absolute()

        self._protected = True


    def __copy__(self) -> "Result":
        """
        Return deepcopy(self).

        Create deep copy.

        """
        return copy.deepcopy(self)

    copy = __copy__


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.

        """
        with h5py.File(self.fname,'r') as f:
            header = [f'Created by {f.attrs["creator"]}',
                      f'        on {f.attrs["created"]}',
                      f' executing "{f.attrs["call"]}"']
        visible_increments = self.visible['increments']

        first = self.view(increments=visible_increments[0:1]).list_data()

        last  = [] if len(visible_increments) < 2 else \
                self.view(increments=visible_increments[-1:]).list_data()

        in_between = [] if len(visible_increments) < 3 else \
                     [f'\n{inc}\n  ...' for inc in visible_increments[1:-1]]

        return util.srepr([util.deemph(header)] + first + in_between + last)


    def _manage_view(self,
                     action: Literal['set', 'add', 'del'],
                     increments: Union[None, int, Sequence[int], str, Sequence[str], bool] = None,
                     times: Union[None, float, Sequence[float], str, Sequence[str], bool] = None,
                     phases: Union[None, str, Sequence[str], bool] = None,
                     homogenizations: Union[None, str, Sequence[str], bool] = None,
                     fields: Union[None, str, Sequence[str], bool] = None) -> "Result":
        """
        Manage the visibility of the groups.

        Parameters
        ----------
        action : str
            Select from 'set', 'add', and 'del'.

        Returns
        -------
        view : damask.Result
            Modified or new view on the DADF5 file.

        """
        if increments is not None and times is not None:
            raise ValueError('"increments" and "times" are mutually exclusive')

        dup = self.copy()
        for what,datasets in zip(['increments','times','phases','homogenizations','fields'],
                                 [ increments,  times,  phases,  homogenizations,  fields ]):
            if  datasets is None:
                continue
            # allow True/False and string arguments
            elif datasets is True:
                datasets = '*'
            elif datasets is False:
                datasets = []
            choice = [datasets] if not hasattr(datasets,'__iter__') or isinstance(datasets,str) else list(datasets) # type: ignore

            if   what == 'increments':
                choice = [c if isinstance(c,str) and c.startswith(prefix_inc) else
                          self.increments[c] if isinstance(c,int) and c<0 else
                          f'{prefix_inc}{c}' for c in choice]
            elif what == 'times':
                atol = 1e-2 * np.min(np.diff(self.times))
                what = 'increments'
                if choice == ['*']:
                    choice = self.increments
                else:
                    iterator = np.array(choice).astype(float)
                    choice = []
                    for c in iterator:
                        idx = np.searchsorted(self.times,c,side='left')
                        if  idx<len(self.times) and np.isclose(c,self.times[idx],rtol=0,atol=atol):
                            choice.append(self.increments[idx])
                        elif idx>0 and np.isclose(c,self.times[idx-1],rtol=0,atol=atol):
                            choice.append(self.increments[idx-1])

            valid = _match(choice,getattr(self,what))
            existing = set(self.visible[what])

            if   action == 'set':
                dup.visible[what] = sorted(set(valid), key=util.natural_sort)
            elif action == 'add':
                dup.visible[what] = sorted(existing.union(valid), key=util.natural_sort)
            elif action == 'del':
                dup.visible[what] = sorted(existing.difference(valid), key=util.natural_sort)

        return dup


    def increments_in_range(self,
                            start: Union[None, str, int] = None,
                            end: Union[None, str, int] = None) -> Sequence[int]:
        """
        Get all increments within a given range.

        Parameters
        ----------
        start : int or str, optional
            Start increment. Defaults to first.
        end : int or str, optional
            End increment. Defaults to last.

        Returns
        -------
        increments : list of ints
            Increment number of all increments within the given bounds.

        """
        s,e = map(lambda x: int(x.split(prefix_inc)[-1] if isinstance(x,str) and x.startswith(prefix_inc) else x),
                  (self.incs[ 0] if start is None else start,
                   self.incs[-1] if  end  is None else  end))
        return [i for i in self.incs if s <= i <= e]

    def times_in_range(self,
                       start: Optional[float] = None,
                       end: Optional[float] = None) -> Sequence[float]:
        """
        Get times of all increments within a given time range.

        Parameters
        ----------
        start : float, optional
            Time of start increment. Defaults to time of first.
        end : float, optional
            Time of end increment. Defaults to time of last.

        Returns
        -------
        times : list of float
            Time of each increment within the given bounds.

        """
        s,e = (self.times[ 0] if start is None else start,
               self.times[-1] if  end  is None else  end)
        return [t for t in self.times if s <= t <= e]


    def view(self,*,
             increments: Union[None, int, Sequence[int], str, Sequence[str], bool] = None,
             times: Union[None, float, Sequence[float], str, Sequence[str], bool] = None,
             phases: Union[None, str, Sequence[str], bool] = None,
             homogenizations: Union[None, str, Sequence[str], bool] = None,
             fields: Union[None, str, Sequence[str], bool] = None,
             protected: Optional[bool] = None) -> "Result":
        """
        Set view.

        Wildcard matching with '?' and '*' is supported.
        True is equivalent to '*', False is equivalent to [].

        Parameters
        ----------
        increments: (list of) int, (list of) str, or bool, optional.
            Numbers of increments to select.
        times: (list of) float, (list of) str, or bool, optional.
            Simulation times of increments to select.
        phases: (list of) str, or bool, optional.
            Names of phases to select.
        homogenizations: (list of) str, or bool, optional.
            Names of homogenizations to select.
        fields: (list of) str, or bool, optional.
            Names of fields to select.
        protected: bool, optional.
            Protection status of existing data.

        Returns
        -------
        view : damask.Result
            View with only the selected attributes being visible.

        Examples
        --------
        Get a view that shows only results from the initial configuration:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r_first = r.view(increments=0)

        Get a view that shows all results between simulation times of 10 to 40:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r_t10to40 = r.view(times=r.times_in_range(10.0,40.0))

        """
        dup = self._manage_view('set',increments,times,phases,homogenizations,fields)
        if protected is not None:
            if not protected:
                print(util.warn('Warning: Modification of existing datasets allowed!'))
            dup._protected = protected

        return dup


    def view_more(self,*,
                  increments: Union[None, int, Sequence[int], str, Sequence[str], bool] = None,
                  times: Union[None, float, Sequence[float], str, Sequence[str], bool] = None,
                  phases: Union[None, str, Sequence[str], bool] = None,
                  homogenizations: Union[None, str, Sequence[str], bool] = None,
                  fields: Union[None, str, Sequence[str], bool] = None) -> "Result":
        """
        Add to view.

        Wildcard matching with '?' and '*' is supported.
        True is equivalent to '*', False is equivalent to [].

        Parameters
        ----------
        increments: (list of) int, (list of) str, or bool, optional.
            Numbers of increments to select.
        times: (list of) float, (list of) str, or bool, optional.
            Simulation times of increments to select.
        phases: (list of) str, or bool, optional.
            Names of phases to select.
        homogenizations: (list of) str, or bool, optional.
            Names of homogenizations to select.
        fields: (list of) str, or bool, optional.
            Names of fields to select.

        Returns
        -------
        modified_view : damask.Result
            View with additional visible attributes.

        Examples
        --------
        Get a view that shows only results from first and last increment:

        >>> import damask
        >>> r_empty = damask.Result('my_file.hdf5').view(increments=False)
        >>> r_first = r_empty.view_more(increments=0)
        >>> r_first_and_last = r.first.view_more(increments=-1)

        """
        return self._manage_view('add',increments,times,phases,homogenizations,fields)


    def view_less(self,*,
                  increments: Union[None, int, Sequence[int], str, Sequence[str], bool] = None,
                  times: Union[None, float, Sequence[float], str, Sequence[str], bool] = None,
                  phases: Union[None, str, Sequence[str], bool] = None,
                  homogenizations: Union[None, str, Sequence[str], bool] = None,
                  fields: Union[None, str, Sequence[str], bool] = None) -> "Result":
        """
        Remove from view.

        Wildcard matching with '?' and '*' is supported.
        True is equivalent to '*', False is equivalent to [].

        Parameters
        ----------
        increments: (list of) int, (list of) str, or bool, optional.
            Numbers of increments to select.
        times: (list of) float, (list of) str, or bool, optional.
            Simulation times of increments to select.
        phases: (list of) str, or bool, optional.
            Names of phases to select.
        homogenizations: (list of) str, or bool, optional.
            Names of homogenizations to select.
        fields: (list of) str, or bool, optional.
            Names of fields to select.

        Returns
        -------
        modified_view : damask.Result
            View with fewer visible attributes.

        Examples
        --------
        Get a view that omits the undeformed configuration:

        >>> import damask
        >>> r_all = damask.Result('my_file.hdf5')
        >>> r_deformed = r_all.view_less(increments=0)

        """
        return self._manage_view('del',increments,times,phases,homogenizations,fields)


    def rename(self,
               name_src: str,
               name_dst: str):
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
        >>> r_unprotected = r.view(protected=False)
        >>> r_unprotected.rename('F','def_grad')

        """
        if self._protected:
            raise PermissionError('rename datasets')

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


    def remove(self, name: str):
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
        >>> r_unprotected = r.view(protected=False)
        >>> r_unprotected.remove('F')

        """
        if self._protected:
            raise PermissionError('delete datasets')

        with h5py.File(self.fname,'a') as f:
            for inc in self.visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            path = '/'.join([inc,ty,label,field,name])
                            if path in f.keys(): del f[path]


    def list_data(self) -> List[str]:
        """
        Collect information on all active datasets in the file.

        Returns
        -------
        data : list of str
            Line-formatted information about active datasets.

        """
        msg = []
        with h5py.File(self.fname,'r') as f:
            for inc in self.visible['increments']:
                msg += [f'\n{inc} ({self.times[self.increments.index(inc)]} s)']
                for ty in ['phase','homogenization']:
                    msg += [f'  {ty}']
                    for label in self.visible[ty+'s']:
                        msg += [f'    {label}']
                        for field in _match(self.visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            msg += [f'      {field}']
                            for d in f['/'.join([inc,ty,label,field])].keys():
                                dataset = f['/'.join([inc,ty,label,field,d])]
                                unit = dataset.attrs["unit"] if h5py3 else \
                                       dataset.attrs["unit"].decode()
                                description = dataset.attrs['description'] if h5py3 else \
                                              dataset.attrs['description'].decode()
                                msg += [f'        {d} / {unit}: {description}']

        return msg


    def enable_user_function(self,
                             func: Callable):
        globals()[func.__name__]=func
        print(f'Function {func.__name__} enabled in add_calculation.')


    @property
    def simulation_setup_files(self):
        """Simulation setup files used to generate the Result object."""
        files = []
        with h5py.File(self.fname,'r') as f_in:
            f_in['setup'].visititems(lambda name,obj: files.append(name) if isinstance(obj,h5py.Dataset) else None)
        return files

    @property
    def incs(self):
        return [int(i.split(prefix_inc)[-1]) for i in self.increments]


    @property
    def coordinates0_point(self) -> np.ndarray:
        """Initial/undeformed cell center coordinates."""
        if self.structured:
            return grid_filters.coordinates0_point(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_p'][()]

    @property
    def coordinates0_node(self) -> np.ndarray:
        """Initial/undeformed nodal coordinates."""
        if self.structured:
            return grid_filters.coordinates0_node(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        else:
            with h5py.File(self.fname,'r') as f:
                return f['geometry/x_n'][()]

    @property
    def geometry0(self) -> VTK:
        """Initial/undeformed geometry."""
        if self.structured:
            return VTK.from_image_data(self.cells,self.size,self.origin)
        else:
            with h5py.File(self.fname,'r') as f:
                return VTK.from_unstructured_grid(f['/geometry/x_n'][()],
                                                  f['/geometry/T_c'][()]-1,
                                                  f['/geometry/T_c'].attrs['VTK_TYPE'] if h5py3 else \
                                                  f['/geometry/T_c'].attrs['VTK_TYPE'].decode())


    @staticmethod
    def _add_absolute(x: Dict[str, Any]) -> Dict[str, Any]:
        return {
                'data':  np.abs(x['data']),
                'label': f'|{x["label"]}|',
                'meta':  {
                          'unit':        x['meta']['unit'],
                          'description': f"absolute value of {x['label']} ({x['meta']['description']})",
                          'creator':     'add_absolute'
                          }
                 }
    def add_absolute(self, x: str):
        """
        Add absolute value.

        Parameters
        ----------
        x : str
            Name of scalar, vector, or tensor dataset to take absolute value of.

        """
        self._add_generic_pointwise(self._add_absolute,{'x':x})


    @staticmethod
    def _add_calculation(**kwargs) -> Dict[str, Any]:
        formula = kwargs['formula']
        for d in re.findall(r'#(.*?)#',formula):
            formula = formula.replace(f'#{d}#',f"kwargs['{d}']['data']")
        data = eval(formula)

        if not hasattr(data,'shape') or data.shape[0] != kwargs[d]['data'].shape[0]:
            raise ValueError('"{}" results in invalid shape'.format(kwargs['formula']))

        return {
                'data':  data,
                'label': kwargs['label'],
                'meta':  {
                          'unit':        kwargs['unit'],
                          'description': f"{kwargs['description']} (formula: {kwargs['formula']})",
                          'creator':     'add_calculation'
                          }
                 }
    def add_calculation(self,
                        formula: str,
                        name: str,
                        unit: str = 'n/a',
                        description: Optional[str] = None):
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
        >>> r.add_calculation('#rho_dip_total#+#rho_mob_total#','rho_total',
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
        dataset_mapping = {d:d for d in set(re.findall(r'#(.*?)#',formula))}                        # datasets used in the formula
        args             = {'formula':formula,'label':name,'unit':unit,'description':description}
        self._add_generic_pointwise(self._add_calculation,dataset_mapping,args)


    @staticmethod
    def _add_stress_Cauchy(P: Dict[str, Any], F: Dict[str, Any]) -> Dict[str, Any]:
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
    def add_stress_Cauchy(self,
                          P: str = 'P',
                          F: str = 'F'):
        """
        Add Cauchy stress calculated from first Piola-Kirchhoff stress and deformation gradient.

        Parameters
        ----------
        P : str, optional
            Name of the dataset containing the first Piola-Kirchhoff stress.
            Defaults to 'P'.
        F : str, optional
            Name of the dataset containing the deformation gradient.
            Defaults to 'F'.

        """
        self._add_generic_pointwise(self._add_stress_Cauchy,{'P':P,'F':F})


    @staticmethod
    def _add_determinant(T: Dict[str, Any]) -> Dict[str, Any]:
        return {
                'data':  np.linalg.det(T['data']),
                'label': f"det({T['label']})",
                'meta':  {
                          'unit':        T['meta']['unit'],
                          'description': f"determinant of tensor {T['label']} ({T['meta']['description']})",
                          'creator':     'add_determinant'
                          }
                }
    def add_determinant(self, T: str):
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
    def _add_deviator(T: Dict[str, Any]) -> Dict[str, Any]:
        return {
                'data':  tensor.deviatoric(T['data']),
                'label': f"s_{T['label']}",
                'meta':  {
                          'unit':        T['meta']['unit'],
                          'description': f"deviator of tensor {T['label']} ({T['meta']['description']})",
                          'creator':     'add_deviator'
                          }
                 }
    def add_deviator(self, T: str):
        """
        Add the deviatoric part of a tensor.

        Parameters
        ----------
        T : str
            Name of tensor dataset.

        Examples
        --------
        Add the deviatoric part of Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_deviator('sigma')

        """
        self._add_generic_pointwise(self._add_deviator,{'T':T})


    @staticmethod
    def _add_eigenvalue(T_sym: Dict[str, Any], eigenvalue: Literal['max, mid, min']) -> Dict[str, Any]:
        if   eigenvalue == 'max':
            label,p = 'maximum',2
        elif eigenvalue == 'mid':
            label,p = 'intermediate',1
        elif eigenvalue == 'min':
            label,p = 'minimum',0
        else:
            raise ValueError(f'invalid eigenvalue: {eigenvalue}')

        return {
                'data': tensor.eigenvalues(T_sym['data'])[:,p],
                'label': f"lambda_{eigenvalue}({T_sym['label']})",
                'meta' : {
                          'unit':        T_sym['meta']['unit'],
                          'description': f"{label} eigenvalue of {T_sym['label']} ({T_sym['meta']['description']})",
                          'creator':     'add_eigenvalue'
                         }
                }
    def add_eigenvalue(self,
                       T_sym: str,
                       eigenvalue: Literal['max', 'mid', 'min'] = 'max'):
        """
        Add eigenvalues of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Name of symmetric tensor dataset.
        eigenvalue : {'max', 'mid', 'min'}, optional
            Eigenvalue. Defaults to 'max'.

        Examples
        --------
        Add the minimum eigenvalue of Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_eigenvalue('sigma','min')

        """
        self._add_generic_pointwise(self._add_eigenvalue,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


    @staticmethod
    def _add_eigenvector(T_sym: Dict[str, Any], eigenvalue: Literal['max', 'mid', 'min']) -> Dict[str, Any]:
        if   eigenvalue == 'max':
            label,p = 'maximum',2
        elif eigenvalue == 'mid':
            label,p = 'intermediate',1
        elif eigenvalue == 'min':
            label,p = 'minimum',0
        else:
            raise ValueError(f'invalid eigenvalue: {eigenvalue}')

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
    def add_eigenvector(self,
                        T_sym: str,
                        eigenvalue: Literal['max', 'mid', 'min'] = 'max'):
        """
        Add eigenvector of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Name of symmetric tensor dataset.
        eigenvalue : {'max', 'mid', 'min'}, optional
            Eigenvalue to which the eigenvector corresponds.
            Defaults to 'max'.

        """
        self._add_generic_pointwise(self._add_eigenvector,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


    @staticmethod
    def _add_IPF_color(l: FloatSequence, q: Dict[str, Any]) -> Dict[str, Any]:
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
    def add_IPF_color(self,
                      l: FloatSequence,
                      q: str = 'O'):
        """
        Add RGB color tuple of inverse pole figure (IPF) color.

        Parameters
        ----------
        l : numpy.array of shape (3)
            Lab frame direction for inverse pole figure.
        q : str, optional
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
    def _add_maximum_shear(T_sym: Dict[str, Any]) -> Dict[str, Any]:
        return {
                'data':  mechanics.maximum_shear(T_sym['data']),
                'label': f"max_shear({T_sym['label']})",
                'meta':  {
                          'unit':        T_sym['meta']['unit'],
                          'description': f"maximum shear component of {T_sym['label']} ({T_sym['meta']['description']})",
                          'creator':     'add_maximum_shear'
                          }
                 }
    def add_maximum_shear(self, T_sym: str):
        """
        Add maximum shear components of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Name of symmetric tensor dataset.

        """
        self._add_generic_pointwise(self._add_maximum_shear,{'T_sym':T_sym})


    @staticmethod
    def _add_equivalent_Mises(T_sym: Dict[str, Any], kind: str) -> Dict[str, Any]:
        k = kind
        if k is None:
            if T_sym['meta']['unit'] == '1':
                k = 'strain'
            elif T_sym['meta']['unit'] == 'Pa':
                k = 'stress'
        if k not in ['stress', 'strain']:
            raise ValueError(f'invalid von Mises kind "{kind}"')

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
    def add_equivalent_Mises(self,
                             T_sym: str,
                             kind: Optional[str] = None):
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
    def _add_norm(x: Dict[str, Any], ord: Union[int, float, Literal['fro', 'nuc']]) ->  Dict[str, Any]:
        o = ord
        if len(x['data'].shape) == 2:
            axis: Union[int, Tuple[int, int]] = 1
            t = 'vector'
            if o is None: o = 2
        elif len(x['data'].shape) == 3:
            axis = (1,2)
            t = 'tensor'
            if o is None: o = 'fro'
        else:
            raise ValueError(f'invalid shape of {x["label"]}')

        return {
                'data':  np.linalg.norm(x['data'],ord=o,axis=axis,keepdims=True),
                'label': f"|{x['label']}|_{o}",
                'meta':  {
                          'unit':        x['meta']['unit'],
                          'description': f"{o}-norm of {t} {x['label']} ({x['meta']['description']})",
                          'creator':     'add_norm'
                          }
                 }
    def add_norm(self,
                 x: str,
                 ord: Union[None, int, float, Literal['fro', 'nuc']] = None):
        """
        Add the norm of a vector or tensor.

        Parameters
        ----------
        x : str
            Name of vector or tensor dataset.
        ord : {non-zero int, inf, -inf, 'fro', 'nuc'}, optional
            Order of the norm. inf means NumPy's inf object. For details refer to numpy.linalg.norm.

        """
        self._add_generic_pointwise(self._add_norm,{'x':x},{'ord':ord})


    @staticmethod
    def _add_stress_second_Piola_Kirchhoff(P: Dict[str, Any], F: Dict[str, Any]) -> Dict[str, Any]:
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
    def add_stress_second_Piola_Kirchhoff(self,
                                          P: str = 'P',
                                          F: str = 'F'):
        r"""
        Add second Piola-Kirchhoff stress calculated from first Piola-Kirchhoff stress and deformation gradient.

        Parameters
        ----------
        P : str, optional
            Name of first Piola-Kirchhoff stress dataset. Defaults to 'P'.
        F : str, optional
            Name of deformation gradient dataset. Defaults to 'F'.

        Notes
        -----
        The definition of the second Piola-Kirchhoff stress
        :math:`\vb{S} = \left(\vb{F}^{-1} \vb{P}\right)_\text{sym}`
        follows the standard definition in nonlinear continuum mechanics.
        As such, no intermediate configuration, for instance that reached by :math:`\vb{F}_\text{p}`,
        is taken into account.

        """
        self._add_generic_pointwise(self._add_stress_second_Piola_Kirchhoff,{'P':P,'F':F})



    @staticmethod
    def _add_pole(q: Dict[str, Any],
                  uvw: FloatSequence,
                  hkl: FloatSequence,
                  with_symmetry: bool,
                  normalize: bool) -> Dict[str, Any]:
        c = q['meta']['c/a'] if 'c/a' in q['meta'] else 1
        brackets = ['[]','()','⟨⟩','{}'][(uvw is None)*1+with_symmetry*2]
        label = 'p^' + '{}{} {} {}{}'.format(brackets[0],
                                              *(uvw if uvw else hkl),
                                              brackets[-1],)
        ori = Orientation(q['data'],lattice=q['meta']['lattice'],a=1,c=c)

        return {
                'data': ori.to_pole(uvw=uvw,hkl=hkl,with_symmetry=with_symmetry,normalize=normalize),
                'label': label,
                'meta' : {
                          'unit':        '1',
                          'description': f'{"normalized " if normalize else ""}lab frame vector along lattice ' \
                                         + ('direction' if uvw is not None else 'plane') \
                                         + ('s' if with_symmetry else ''),
                          'creator':     'add_pole'
                          }
                }
    def add_pole(self,
                 q: str = 'O',
                 *,
                 uvw: Optional[FloatSequence] = None,
                 hkl: Optional[FloatSequence] = None,
                 with_symmetry: bool = False,
                 normalize: bool = True):
        """
        Add lab frame vector along lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        q : str, optional
            Name of the dataset containing the crystallographic orientation as quaternions.
            Defaults to 'O'.
        uvw|hkl : numpy.ndarray of shape (3)
            Miller indices of crystallographic direction or plane normal.
        with_symmetry : bool, optional
            Calculate all N symmetrically equivalent vectors.
            Defaults to True.
        normalize : bool, optional
            Normalize output vector.
            Defaults to True.

        """
        self._add_generic_pointwise(self._add_pole,
                                    {'q':q},
                                    {'uvw':uvw,'hkl':hkl,'with_symmetry':with_symmetry,'normalize':normalize})


    @staticmethod
    def _add_rotation(F: Dict[str, Any]) -> Dict[str, Any]:
        return {
                'data':  mechanics.rotation(F['data']).as_matrix(),
                'label': f"R({F['label']})",
                'meta':  {
                          'unit':        F['meta']['unit'],
                          'description': f"rotational part of {F['label']} ({F['meta']['description']})",
                          'creator':     'add_rotation'
                          }
                 }
    def add_rotation(self, F: str):
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
    def _add_spherical(T: Dict[str, Any]) -> Dict[str, Any]:
        return {
                'data':  tensor.spherical(T['data'],False),
                'label': f"p_{T['label']}",
                'meta':  {
                          'unit':        T['meta']['unit'],
                          'description': f"spherical component of tensor {T['label']} ({T['meta']['description']})",
                          'creator':     'add_spherical'
                          }
                 }
    def add_spherical(self, T: str):
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
    def _add_strain(F: Dict[str, Any], t: Literal['V', 'U'], m: float) -> Dict[str, Any]:
        side = 'left' if t == 'V' else 'right'
        return {
                'data':  mechanics.strain(F['data'],t,m),
                'label': f"epsilon_{t}^{m}({F['label']})",
                'meta':  {
                          'unit':        F['meta']['unit'],
                          'description': f'Seth-Hill strain tensor of order {m} based on {side} stretch tensor '+\
                                         f"of {F['label']} ({F['meta']['description']})",
                          'creator':     'add_strain'
                          }
                 }
    def add_strain(self,
                   F: str = 'F',
                   t: Literal['V', 'U'] = 'V',
                   m: float = 0.0):
        r"""
        Add strain tensor (Seth-Hill family) of a deformation gradient.

        By default, the logarithmic strain based on the
        left stretch tensor is added.

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
        Add the Euler-Almansi strain:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_strain(t='V',m=-1.0)

        Add the plastic Biot strain:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_strain('F_p','U',0.5)

        Notes
        -----
        The presence of rotational parts in the elastic and plastic deformation gradient
        calls for the use of
        material/Lagragian strain measures (based on 'U') for plastic strains and
        spatial/Eulerian strain measures (based on 'V') for elastic strains
        when calculating averages.

        The strain is defined as:

        .. math::

            m = 0 \\\\
            \vb*{\epsilon}_V^{(0)} = \ln (\vb{V}) \\\\
            \vb*{\epsilon}_U^{(0)} = \ln (\vb{U}) \\\\

            m \neq 0 \\\\
            \vb*{\epsilon}_V^{(m)} = \frac{1}{2m} (\vb{V}^{2m} - \vb{I}) \\\\
            \vb*{\epsilon}_U^{(m)} = \frac{1}{2m} (\vb{U}^{2m} - \vb{I})

        References
        ----------
        | https://en.wikipedia.org/wiki/Finite_strain_theory
        | https://de.wikipedia.org/wiki/Verzerrungstensor

        """
        self._add_generic_pointwise(self._add_strain,{'F':F},{'t':t,'m':m})


    @staticmethod
    def _add_stretch_tensor(F: Dict[str, Any], t: str) -> Dict[str, Any]:
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
    def add_stretch_tensor(self,
                           F: str = 'F',
                           t: Literal['V', 'U'] = 'V'):
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
    def _add_curl(f: Dict[str, Any], size: np.ndarray) -> Dict[str, Any]:
        return {
                'data':  grid_filters.curl(size,f['data']),
                'label': f"curl({f['label']})",
                'meta':  {
                          'unit':        f['meta']['unit']+'/m',
                          'description': f"curl of {f['label']} ({f['meta']['description']})",
                          'creator':     'add_curl'
                          }
                 }
    def add_curl(self, f: str):
        """
        Add curl of a field.

        Parameters
        ----------
        f : str
            Name of vector or tensor field dataset.

        Notes
        -----
        This function is only available for structured grids,
        i.e. fields resulting from the grid solver.

        """
        self._add_generic_grid(self._add_curl,{'f':f},{'size':self.size})


    @staticmethod
    def _add_divergence(f: Dict[str, Any], size: np.ndarray) -> Dict[str, Any]:
        return {
                'data':  grid_filters.divergence(size,f['data']),
                'label': f"divergence({f['label']})",
                'meta':  {
                          'unit':        f['meta']['unit']+'/m',
                          'description': f"divergence of {f['label']} ({f['meta']['description']})",
                          'creator':     'add_divergence'
                          }
                 }
    def add_divergence(self, f: str):
        """
        Add divergence of a field.

        Parameters
        ----------
        f : str
            Name of vector or tensor field dataset.

        Notes
        -----
        This function is only available for structured grids,
        i.e. fields resulting from the grid solver.

        """
        self._add_generic_grid(self._add_divergence,{'f':f},{'size':self.size})


    @staticmethod
    def _add_gradient(f: Dict[str, Any], size: np.ndarray) -> Dict[str, Any]:
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
    def add_gradient(self, f: str):
        """
        Add gradient of a field.

        Parameters
        ----------
        f : str
            Name of scalar or vector field dataset.

        Notes
        -----
        This function is only available for structured grids,
        i.e. fields resulting from the grid solver.

        """
        self._add_generic_grid(self._add_gradient,{'f':f},{'size':self.size})


    def _add_generic_grid(self,
                          func: Callable,
                          datasets: Dict[str, str],
                          args: Dict[str, str] = {},
                          constituents = None):
        """
        General function to add data on a regular grid.

        Parameters
        ----------
        func : function
            Callback function that calculates a new dataset from one or
            more datasets per DADF5 group.
        datasets : dictionary
            Details of the datasets to be used:
            {arg (name to which the data is passed in func): label (in DADF5 file)}.
        args : dictionary, optional
            Arguments parsed to func.

        """
        if len(datasets) != 1 or self.N_constituents != 1:
            raise NotImplementedError

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        increments = self.place(list(datasets.values()),False)
        if not increments: raise RuntimeError("received invalid dataset")
        with h5py.File(self.fname, 'a') as f:
            for increment in increments.items():
                for ty in increment[1].items():
                    for field in ty[1].items():
                        d: np.ma.MaskedArray = list(field[1].values())[0]
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
                            h5_dataset = f[path].create_dataset(r['label'],data=result1)

                            now = datetime.datetime.now().astimezone()
                            h5_dataset.attrs['created'] = now.strftime('%Y-%m-%d %H:%M:%S%z') if h5py3 else \
                                                          now.strftime('%Y-%m-%d %H:%M:%S%z').encode()

                            for l,v in r['meta'].items():
                                h5_dataset.attrs[l.lower()]=v if h5py3 else v.encode()
                            creator = h5_dataset.attrs['creator'] if h5py3 else \
                                      h5_dataset.attrs['creator'].decode()
                            h5_dataset.attrs['creator'] = f'damask.Result.{creator} v{damask.version}' if h5py3 else \
                                                          f'damask.Result.{creator} v{damask.version}'.encode()


    def _job_pointwise(self,
                       group: str,
                       callback: Callable,
                       datasets: Dict[str, str],
                       args: Dict[str, str],
                       lock: Lock) -> List[Union[None, Any]]:
        """Execute job for _add_generic_pointwise."""
        try:
            datasets_in = {}
            lock.acquire()
            with h5py.File(self.fname,'r') as f:
                for arg,label in datasets.items():
                    loc  = f[group+'/'+label]
                    datasets_in[arg]={'data' :loc[()],
                                      'label':label,
                                      'meta': {k:(v.decode() if not h5py3 and type(v) is bytes else v) \
                                               for k,v in loc.attrs.items()}}
            lock.release()
            r = callback(**datasets_in,**args)
            return [group,r]
        except Exception as err:
            print(f'Error during calculation: {err}.')
            return [None,None]


    def _add_generic_pointwise(self,
                               func: Callable,
                               datasets: Dict[str, Any],
                               args: Dict[str, Any] = {}):
        """
        General function to add pointwise data.

        Parameters
        ----------
        callback : function
            Callback function that calculates a new dataset from one or
            more datasets per DADF5 group.
        datasets : dictionary
            Details of the datasets to be used:
            {arg (name to which the data is passed in func): label (in DADF5 file)}.
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

        default_arg = partial(self._job_pointwise,callback=func,datasets=datasets,args=args,lock=lock)

        for group,result in util.show_progress(pool.imap_unordered(default_arg,groups),len(groups)):# type: ignore
            if not result:
                continue
            lock.acquire()
            with h5py.File(self.fname, 'a') as f:
                try:
                    if not self._protected and '/'.join([group,result['label']]) in f:
                        dataset = f['/'.join([group,result['label']])]
                        dataset[...] = result['data']
                        dataset.attrs['overwritten'] = True
                    else:
                        shape = result['data'].shape
                        if compress := result['data'].size >= chunk_size*2:
                            chunks = (chunk_size//np.prod(shape[1:]),)+shape[1:]
                        else:
                            chunks = shape
                        dataset = f[group].create_dataset(result['label'],data=result['data'],
                                                          maxshape=shape, chunks=chunks,
                                                          compression = 'gzip' if compress else None,
                                                          compression_opts = 6 if compress else None,
                                                          shuffle=True,fletcher32=True)

                    now = datetime.datetime.now().astimezone()
                    dataset.attrs['created'] = now.strftime('%Y-%m-%d %H:%M:%S%z') if h5py3 else \
                                               now.strftime('%Y-%m-%d %H:%M:%S%z').encode()

                    for l,v in result['meta'].items():
                        dataset.attrs[l.lower()]=v.encode() if not h5py3 and type(v) is str else v
                    creator = dataset.attrs['creator'] if h5py3 else \
                              dataset.attrs['creator'].decode()
                    dataset.attrs['creator'] = f'damask.Result.{creator} v{damask.version}' if h5py3 else \
                                               f'damask.Result.{creator} v{damask.version}'.encode()

                except (OSError,RuntimeError) as err:
                    print(f'Could not add dataset: {err}.')
            lock.release()

        pool.close()
        pool.join()


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


    def get(self,
            output: Union[str, List[str]] = '*',
            flatten: bool = True,
            prune: bool = True) -> Optional[Dict[str,Any]]:
        """
        Collect data per phase/homogenization reflecting the group/folder structure in the DADF5 file.

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets to read.
            Defaults to '*', in which case all datasets are read.
        flatten : bool, optional
            Remove singular levels of the folder hierarchy.
            This might be beneficial in case of single increment,
            phase/homogenization, or field. Defaults to True.
        prune : bool, optional
            Remove branches with no data. Defaults to True.

        Returns
        -------
        data : dict of numpy.ndarray
            Datasets structured by phase/homogenization and according to selected view.

        """
        r: Dict[str,Any] = {}

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


    def place(self,
              output: Union[str, List[str]] = '*',
              flatten: bool = True,
              prune: bool = True,
              constituents: Optional[IntSequence] = None,
              fill_float: float = np.nan,
              fill_int: int = 0) -> Optional[Dict[str,Any]]:
        """
        Merge data into spatial order that is compatible with the damask.VTK geometry representation.

        The returned data structure reflects the group/folder structure in the DADF5 file.

        Multi-phase data is fused into a single output.
        `place` is equivalent to `get` if only one phase/homogenization
        and one constituent is present.

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets to read.
            Defaults to '*', in which case all visible datasets are placed.
        flatten : bool, optional
            Remove singular levels of the folder hierarchy.
            This might be beneficial in case of single increment or field.
            Defaults to True.
        prune : bool, optional
            Remove branches with no data. Defaults to True.
        constituents : (list of) int, optional
            Constituents to consider.
            Defaults to None, in which case all constituents are considered.
        fill_float : float, optional
            Fill value for non-existent entries of floating point type.
            Defaults to NaN.
        fill_int : int, optional
            Fill value for non-existent entries of integer type.
            Defaults to 0.

        Returns
        -------
        data : dict of numpy.ma.MaskedArray
            Datasets structured by spatial position and according to selected view.

        """
        r: Dict[str,Any] = {}

        constituents_ = map(int,constituents) if isinstance(constituents,Iterable) else \
                      (range(self.N_constituents) if constituents is None else [constituents])      # type: ignore

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


    def export_XDMF(self,
                    output: Union[str, List[str]] = '*',
                    target_dir: Union[None, str, Path] = None,
                    absolute_path: bool = False):
        """
        Write XDMF file to directly visualize data from DADF5 file.

        The XDMF format is only supported for structured grids
        with single phase and single constituent.
        For other cases use `export_VTK`.

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets included in the XDMF file.
            Defaults to '*', in which case all datasets are considered.
        target_dir : str or pathlib.Path, optional
            Directory to save XDMF file. Will be created if non-existent.
        absolute_path : bool, optional
            Store absolute (instead of relative) path to DADF5 file.
            Defaults to False, i.e. the XDMF file expects the
            DADF5 file at a stable relative path.

        """
        if self.N_constituents != 1 or len(self.phases) != 1 or not self.structured:
            raise TypeError('XDMF output requires structured grid with single phase and single constituent.')

        attribute_type_map = defaultdict(lambda:'Matrix', ( ((),'Scalar'), ((3,),'Vector'), ((3,3),'Tensor')) )

        def number_type_map(dtype):
            if dtype in np.sctypes['int']:   return 'Int'
            if dtype in np.sctypes['uint']:  return 'UInt'
            if dtype in np.sctypes['float']: return 'Float'


        xdmf = ET.Element('Xdmf')
        xdmf.attrib = {'Version':  '2.0',
                       'xmlns:xi': 'http://www.w3.org/2001/XInclude'}

        domain = ET.SubElement(xdmf, 'Domain')

        collection = ET.SubElement(domain, 'Grid')
        collection.attrib = {'GridType':       'Collection',
                             'CollectionType': 'Temporal',
                             'Name':           'Increments'}

        time = ET.SubElement(collection, 'Time')
        time.attrib = {'TimeType': 'List'}

        time_data = ET.SubElement(time, 'DataItem')
        times = [self.times[self.increments.index(i)] for i in self.visible['increments']]
        time_data.attrib = {'Format':     'XML',
                            'NumberType': 'Float',
                            'Dimensions': f'{len(times)}'}
        time_data.text = ' '.join(map(str,times))

        attributes = []
        data_items = []

        hdf5_name = self.fname.name
        hdf5_dir  = self.fname.parent
        xdmf_dir  = Path.cwd() if target_dir is None else Path(target_dir)
        hdf5_link = (hdf5_dir if absolute_path else Path(os.path.relpath(hdf5_dir,xdmf_dir.resolve())))/hdf5_name

        with h5py.File(self.fname,'r') as f:
            for inc in self.visible['increments']:

                grid = ET.SubElement(collection,'Grid')
                grid.attrib = {'GridType': 'Uniform',
                               'Name':      inc}

                topology = ET.SubElement(grid, 'Topology')
                topology.attrib = {'TopologyType': '3DCoRectMesh',
                                   'Dimensions':   '{} {} {}'.format(*(self.cells[::-1]+1))}

                geometry = ET.SubElement(grid, 'Geometry')
                geometry.attrib = {'GeometryType':'Origin_DxDyDz'}

                origin = ET.SubElement(geometry, 'DataItem')
                origin.attrib = {'Format':     'XML',
                                 'NumberType': 'Float',
                                 'Dimensions': '3'}
                origin.text = "{} {} {}".format(*self.origin[::-1])

                delta = ET.SubElement(geometry, 'DataItem')
                delta.attrib = {'Format':     'XML',
                                'NumberType': 'Float',
                                'Dimensions': '3'}
                delta.text="{} {} {}".format(*(self.size/self.cells)[::-1])

                attributes.append(ET.SubElement(grid, 'Attribute'))
                attributes[-1].attrib = {'Name':          'u / m',
                                         'Center':        'Node',
                                         'AttributeType': 'Vector'}
                data_items.append(ET.SubElement(attributes[-1], 'DataItem'))
                data_items[-1].attrib = {'Format':     'HDF',
                                         'Precision':  '8',
                                         'Dimensions': '{} {} {} 3'.format(*(self.cells[::-1]+1))}
                data_items[-1].text = f'{hdf5_link}:/{inc}/geometry/u_n'
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
                                                         'Dimensions': '{} {} {} {}'.format(*self.cells[::-1],1 if shape == () else
                                                                                                        np.prod(shape))}
                                data_items[-1].text = f'{hdf5_link}:{name}'

        xdmf_dir.mkdir(parents=True,exist_ok=True)
        with util.open_text((xdmf_dir/hdf5_name).with_suffix('.xdmf'),'w') as f:
            f.write(xml.dom.minidom.parseString(ET.tostring(xdmf).decode()).toprettyxml())


    def export_VTK(self,
                   output: Union[str,List[str]] = '*',
                   mode: str = 'cell',
                   constituents: Optional[IntSequence] = None,
                   target_dir: Union[None, str, Path] = None,
                   fill_float: float = np.nan,
                   fill_int: int = 0,
                   parallel: bool = True):
        """
        Export to VTK cell/point data.

        One VTK file per visible increment is created.
        For point data, the VTK format is PolyData (.vtp).
        For cell data, the file format is either ImageData (.vti)
        or UnstructuredGrid (.vtu) for grid-based or mesh-based simulations,
        respectively.

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets to export to the VTK file.
            Defaults to '*', in which case all visible datasets are exported.
        mode : {'cell', 'point'}, optional
            Export in cell format or point format.
            Defaults to 'cell'.
        constituents : (list of) int, optional
            Constituents to consider.
            Defaults to None, in which case all constituents are considered.
        target_dir : str or pathlib.Path, optional
            Directory to save VTK files. Will be created if non-existent.
        fill_float : float, optional
            Fill value for non-existent entries of floating point type.
            Defaults to NaN.
        fill_int : int, optional
            Fill value for non-existent entries of integer type.
            Defaults to 0.
        parallel : bool, optional
            Write VTK files in parallel in a separate background process.
            Defaults to True.

        """
        if mode.lower()=='cell':
            v = self.geometry0
        elif mode.lower()=='point':
            v = VTK.from_poly_data(self.coordinates0_point)
        else:
            raise ValueError(f'invalid mode "{mode}"')

        v.comments = [util.execution_stamp('Result','export_VTK')]

        N_digits = int(np.floor(np.log10(max(1,self.incs[-1]))))+1

        constituents_ = constituents if isinstance(constituents,Iterable) else \
                        (range(self.N_constituents) if constituents is None else [constituents])    # type: ignore

        suffixes = [''] if self.N_constituents == 1 or isinstance(constituents,int) else \
                   [f'#{c}' for c in constituents_]

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        vtk_dir = Path.cwd() if target_dir is None else Path(target_dir)
        vtk_dir.mkdir(parents=True,exist_ok=True)

        with h5py.File(self.fname,'r') as f:
            if self.version_minor >= 13:
                creator = f.attrs['creator'] if h5py3 else f.attrs['creator'].decode()
                created = f.attrs['created'] if h5py3 else f.attrs['created'].decode()
                v.comments += [f'{creator} ({created})']

            for inc in util.show_progress(self.visible['increments']):

                u = _read(f['/'.join([inc,'geometry','u_n' if mode.lower() == 'cell' else 'u_p'])])
                v = v.set('u',u)

                for ty in ['phase','homogenization']:
                    for field in self.visible['fields']:
                        outs: Dict[str, np.ma.core.MaskedArray] = {}
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
                            v = v.set(' / '.join(['/'.join([ty,field,label]),dataset.dtype.metadata['unit']]),dataset)


                v.save(vtk_dir/f'{self.fname.stem}_inc{inc.split(prefix_inc)[-1].zfill(N_digits)}',
                       parallel=parallel)

    def export_DADF5(self,
                     fname,
                     output: Union[str, List[str]] = '*',
                     mapping = None):
        """
        Export visible components into a new DADF5 file.

        A DADF5 (DAMASK HDF5) file contains DAMASK results.
        Its group/folder structure reflects the layout in material.yaml.

        Parameters
        ----------
        fname : str or pathlib.Path
            Name of the DADF5 file to be created.
        output : (list of) str, optional
            Names of the datasets to export.
            Defaults to '*', in which case all visible datasets are exported.
        mapping : numpy.ndarray of int, shape (:,:,:), optional
            Indices for regridding.

        """
        if Path(fname).expanduser().absolute() == self.fname:
            raise PermissionError(f'cannot overwrite {self.fname}')

        def cp(path_in,path_out,label,mapping):
            if mapping is None:
                path_in.copy(label,path_out)
            else:
                path_out.create_dataset(label,data=path_in[label][()][mapping])
                path_out[label].attrs.update(path_in[label].attrs)


        with h5py.File(self.fname,'r') as f_in, h5py.File(fname,'w') as f_out:
            f_out.attrs.update(f_in.attrs)
            for g in ['setup','geometry'] + (['cell_to'] if mapping is None else []):
                f_in.copy(g,f_out)

            if mapping is not None:
                cells = mapping.shape
                mapping_flat = mapping.flatten(order='F')
                f_out['geometry'].attrs['cells'] = cells
                f_out.create_group('cell_to')                                                       # ToDo: attribute missing
                mappings = {'phase':{},'homogenization':{}}                                         # type: ignore

                mapping_phase = f_in['cell_to']['phase'][()][mapping_flat]
                for p in np.unique(mapping_phase['label']):
                    m = mapping_phase['label'] == p
                    mappings['phase'][p] = mapping_phase[m]['entry']
                    c = np.count_nonzero(m)
                    mapping_phase[m] = list(zip((p,)*c,tuple(np.arange(c))))
                f_out['cell_to'].create_dataset('phase',data=mapping_phase.reshape(np.prod(mapping_flat.shape),-1))

                mapping_homog = f_in['cell_to']['homogenization'][()][mapping]
                for h in np.unique(mapping_homog['label']):
                    m = mapping_homog['label'] == h
                    mappings['homogenization'][h] = mapping_homog[m]['entry']
                    c = np.count_nonzero(m)
                    mapping_homog[mapping_homog['label'] == h] = list(zip((h,)*c,tuple(np.arange(c))))
                f_out['cell_to'].create_dataset('homogenization',data=mapping_homog.flatten())


            for inc in util.show_progress(self.visible['increments']):
                f_in.copy(inc,f_out,shallow=True)
                if mapping is None:
                    for label in ['u_p','u_n']:
                        f_in[inc]['geometry'].copy(label,f_out[inc]['geometry'])
                else:
                    u_p = f_in[inc]['geometry']['u_p'][()][mapping_flat]
                    f_out[inc]['geometry'].create_dataset('u_p',data=u_p)
                    u_n = np.zeros((len(mapping_flat),3))                                           # ToDo: needs implementation
                    f_out[inc]['geometry'].create_dataset('u_n',data=u_n)


                for label in self.homogenizations:
                    f_in[inc]['homogenization'].copy(label,f_out[inc]['homogenization'],shallow=True)
                for label in self.phases:
                    f_in[inc]['phase'].copy(label,f_out[inc]['phase'],shallow=True)

                for ty in ['phase','homogenization']:
                    for label in self.visible[ty+'s']:
                        for field in _match(self.visible['fields'],f_in['/'.join([inc,ty,label])].keys()):
                            p = '/'.join([inc,ty,label,field])
                            for out in _match(output,f_in[p].keys()):
                                cp(f_in[p],f_out[p],out,None if mapping is None else mappings[ty][label.encode()])


    def export_simulation_setup(self,
                     output: Union[str, List[str]] = '*',
                     target_dir: Union[None, str, Path] = None,
                     overwrite: bool = False,
                     ):
        """
        Export original simulation setup of the Result object.

        Parameters
        ----------
        output : (list of) str, optional
            Names of the datasets to export to the file.
            Defaults to '*', in which case all setup files are exported.
        target_dir : str or pathlib.Path, optional
            Directory to save setup files. Will be created if non-existent.
        overwrite : bool, optional
            Overwrite any existing setup files.
            Defaults to False.

        """
        def export(name: str,
                   obj: Union[h5py.Dataset,h5py.Group],
                   output: Union[str,List[str]],
                   cfg_dir: Path,
                   overwrite: bool):

            cfg = cfg_dir/name

            if type(obj) == h5py.Dataset and _match(output,[name]):
                if cfg.exists() and not overwrite:
                    raise PermissionError(f'"{cfg}" exists')
                else:
                    cfg.parent.mkdir(parents=True,exist_ok=True)
                    with util.open_text(cfg,'w') as f_out: f_out.write(obj[0].decode())

        cfg_dir = (Path.cwd() if target_dir is None else Path(target_dir))
        with h5py.File(self.fname,'r') as f_in:
            f_in['setup'].visititems(partial(export,
                                             output=output,
                                             cfg_dir=cfg_dir,
                                             overwrite=overwrite))
