import re
import fnmatch
import os
import copy
import xml.etree.ElementTree as ET                                                                  # noqa
import xml.dom.minidom
import functools
from pathlib import Path
from collections import defaultdict
from collections.abc import Iterable
from typing import Optional, Union, Callable, Any, Sequence, Literal, NamedTuple

import h5py
import numpy as np
from numpy import ma
from scipy import interpolate
import logging

import damask
from . import VTK
from . import Orientation
from . import Rotation
from . import grid_filters
from . import mechanics
from . import tensor
from . import util
from ._typehints import FloatSequence, IntSequence, DADF5Dataset, BravaisLattice


logger = logging.getLogger(__name__)

chunk_size = 1024**2//8                                                                             # for compression in HDF5

prefix_inc = 'increment_'

class MappingsTuple(NamedTuple):
    at_cell_ph: list[dict[str, np.ndarray]]
    in_data_ph: list[dict[str, np.ndarray]]
    at_cell_ho: dict[str, np.ndarray]
    in_data_ho: dict[str, np.ndarray]


def _read(dataset: h5py._hl.dataset.Dataset) -> np.ndarray:
    """Read a dataset and its metadata into a numpy.ndarray."""
    dtype = np.dtype(dataset.dtype, metadata=dict(dataset.attrs.items()))                           # type: ignore
    return np.array(dataset, dtype=dtype)

def _read_dt(dataset: h5py._hl.dataset.Dataset) -> np.dtype:
    """Only read the metadata of an item without loading the full array."""
    return np.dtype(dataset.dtype, metadata=dict(dataset.attrs.items()))

def _get_common_metadata(dtypes: list[np.dtype]) -> np.dtype:
    metadata_list = [dtype.metadata for dtype in dtypes if dtype.metadata]
    common_metadata = {}
    for key in set.intersection(*map(set, metadata_list)):
        value = metadata_list[0][key]
        if all(np.array_equal(meta[key], value) for meta in metadata_list):
            common_metadata[key] = value
    return np.dtype(dtypes[0].base.type, metadata=common_metadata)

def _match(requested,
           existing: h5py._hl.base.KeysViewHDF5) -> list[str]:
    """Find matches among two sets of labels."""
    def flatten_list(list_of_lists):
        return [e for e_ in list_of_lists for e in e_]

    if requested is True:
        requested = '*'
    elif requested is False or requested is None:
        requested = []

    requested_ = util.to_list(requested)

    return sorted(set(flatten_list([fnmatch.filter(existing,r) for r in requested_])),
                  key=util.natural_sort)

def _empty_like(dataset_shape: tuple[int],
                dtype: np.dtype,
                N_materialpoints: int,
                fill_float: float,
                fill_int: int) -> np.ma.core.MaskedArray:
    """Create empty numpy.ma.MaskedArray."""
    return ma.array(np.empty((N_materialpoints,) + dataset_shape[1:], dtype=dtype),
                    fill_value=fill_float if np.issubdtype(dtype, np.floating) else fill_int,
                    mask=True)


class Result:
    r"""
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
    and first Piola-Kirchhoff stress 'P', add the von Mises equivalent of the
    Cauchy stress, and export it to VTK (file) and numpy.ndarray (memory).

    >>> import damask
    >>> r = damask.Result('my_file.hdf5')
    >>> r
    \x1b[2mCreated by DAMASK_grid ...
            on ...
     executing "..."\x1b[0m
    ...
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

            r = re.compile(rf'{prefix_inc}([0-9]+)')
            self._increments = sorted([i for i in f.keys() if r.match(i)],key=util.natural_sort)
            self._times = {int(i.split('_')[1]):np.around(f[i].attrs['t/s'],12) for i in self._increments}
            if len(self._increments) == 0:
                raise ValueError('incomplete DADF5 file')

            self.N_materialpoints, self.N_constituents = np.shape(f['cell_to/phase'])

            self.homogenization   = f['cell_to/homogenization']['label'].astype('str')
            self._homogenizations = sorted(np.unique(self.homogenization),key=util.natural_sort)
            self.phase            = f['cell_to/phase']['label'].astype('str')
            self._phases          = sorted(np.unique(self.phase),key=util.natural_sort)

            fields: list[str] = []
            for c in self._phases:
                fields += f['/'.join([self._increments[0],'phase',c])].keys()
            for m in self._homogenizations:
                fields += f['/'.join([self._increments[0],'homogenization',m])].keys()
            self._fields = sorted(set(fields),key=util.natural_sort)                                # make unique

        self._visible = {'increments':      self._increments,
                         'phases':          self._phases,
                         'homogenizations': self._homogenizations,
                         'fields':          self._fields,
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
        visible_increments = self._visible['increments']

        first = self.view(increments=visible_increments[0:1]).list_data()

        last  = [] if len(visible_increments) < 2 else \
                self.view(increments=visible_increments[-1:]).list_data()

        in_between = [] if len(visible_increments) < 3 else \
                     [f'\n{inc}\n  ...' for inc in visible_increments[1:-1]]

        return util.srepr([util.deemph(header)] + first + in_between + last)


    def __bool__(self) -> bool:
        """
        Return bool(self).

        Report whether file is still a valid HDF5 file.
        """
        return h5py.is_hdf5(self.fname)


    def _manage_view(self,
                     action: Literal['set', 'add', 'del'],
                     increments: Optional[Union[int, Sequence[int], str, Sequence[str], bool]] = None,
                     times: Optional[Union[float, Sequence[float], Literal['*'], bool]] = None,
                     phases: Optional[Union[str, Sequence[str], bool]] = None,
                     homogenizations: Optional[Union[str, Sequence[str], bool]] = None,
                     fields: Optional[Union[str, Sequence[str], bool]] = None) -> "Result":
        """
        Manage the visibility of the groups.

        Parameters
        ----------
        action : str
            Select from 'set', 'add', and 'del'.
        increments : (list of) int, (list of) str, or bool, optional
            Numbers of increments to select. Mutually exclusive with 'times'.
        times : (list of) float, '*', or bool, optional
            Simulation times of increments to select. Mutually exclusive with
            'increments'.
        phases : (list of) str, or bool, optional
            Names of phases to select.
        homogenizations : (list of) str, or bool, optional
            Names of homogenizations to select.
        fields : (list of) str, or bool, optional
            Names of fields to select.
        protected : bool, optional
            Protection status of existing data.

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
            choice = util.to_list(datasets)
            N_expected = len(choice)

            if   what == 'increments':
                choice = [c if isinstance(c,str) and c.startswith(prefix_inc) else
                          self._increments[c] if isinstance(c,int) and c<0 else
                          f'{prefix_inc}{c}' for c in choice]
            elif what == 'times':
                times = list(self._times.values())
                atol = 1.e-2 * np.min(np.diff(times))
                what = 'increments'
                if choice == ['*']:
                    choice = self._increments
                else:
                    iterator = np.array(choice).astype(float)
                    choice = []
                    for c in iterator:
                        idx = np.searchsorted(times,c,side='left')
                        if  idx<len(self._times) and np.isclose(c,times[idx],rtol=0,atol=atol):
                            choice.append(self._increments[idx])
                        elif idx>0 and np.isclose(c,times[idx-1],rtol=0,atol=atol):
                            choice.append(self._increments[idx-1])

            valid = _match(choice,getattr(self,'_'+what))
            if len(valid) < N_expected:
                w = what if times is None else 'times'
                logger.warning(f'Found only "{list(map(str,valid))}" when requesting "{datasets}" for "{w}".')

            existing = set(self._visible[what])
            if   action == 'set':
                dup._visible[what] = sorted(set(valid), key=util.natural_sort)
            elif action == 'add':
                dup._visible[what] = sorted(existing.union(valid), key=util.natural_sort)
            elif action == 'del':
                dup._visible[what] = sorted(existing.difference(valid), key=util.natural_sort)

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
                  (self._incs[ 0] if start is None else start,
                   self._incs[-1] if  end  is None else  end))
        return [i for i in self._incs if s <= i <= e]

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
               self.times[-1] if end   is None else end)
        return [t for t in self.times if s <= t <= e]


    def view(self,*,
             increments: Optional[Union[int, Sequence[int], str, Sequence[str], bool]] = None,
             times: Optional[Union[float, Sequence[float], Literal['*'], bool]] = None,
             phases: Optional[Union[str, Sequence[str], bool]] = None,
             homogenizations: Optional[Union[str, Sequence[str], bool]] = None,
             fields: Optional[Union[str, Sequence[str], bool]] = None,
             protected: Optional[bool] = None) -> "Result":
        """
        Set view.

        Wildcard matching with '?' and '*' is supported.
        True is equivalent to '*', False is equivalent to [].

        Parameters
        ----------
        increments : (list of) int, (list of) str, or bool, optional
            Numbers of increments to select. Mutually exclusive with 'times'.
        times : (list of) float, '*', or bool, optional
            Simulation times of increments to select. Mutually exclusive with
            'increments'.
        phases : (list of) str, or bool, optional
            Names of phases to select.
        homogenizations : (list of) str, or bool, optional
            Names of homogenizations to select.
        fields : (list of) str, or bool, optional
            Names of fields to select.
        protected : bool, optional
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
                logger.warning(util.warn('Modification of existing datasets allowed!'))
            dup._protected = protected

        return dup


    def view_more(self,*,
                  increments: Optional[Union[int, Sequence[int], str, Sequence[str], bool]] = None,
                  times: Optional[Union[float, Sequence[float], Literal['*'], bool]] = None,
                  phases: Optional[Union[str, Sequence[str], bool]] = None,
                  homogenizations: Optional[Union[str, Sequence[str], bool]] = None,
                  fields: Optional[Union[str, Sequence[str], bool]] = None) -> "Result":
        """
        Add to view.

        Wildcard matching with '?' and '*' is supported.
        True is equivalent to '*', False is equivalent to [].

        Parameters
        ----------
        increments : (list of) int, (list of) str, or bool, optional
            Numbers of increments to select. Mutually exclusive with 'times'.
        times : (list of) float, '*', or bool, optional
            Simulation times of increments to select. Mutually exclusive with
            'increments'.
        phases : (list of) str, or bool, optional
            Names of phases to select.
        homogenizations : (list of) str, or bool, optional
            Names of homogenizations to select.
        fields : (list of) str, or bool, optional
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
        >>> r_first_and_last = r_first.view_more(increments=-1)
        """
        return self._manage_view('add',increments,times,phases,homogenizations,fields)


    def view_less(self,*,
                  increments: Optional[Union[int, Sequence[int], str, Sequence[str], bool]] = None,
                  times: Optional[Union[float, Sequence[float], Literal['*'], bool]] = None,
                  phases: Optional[Union[str, Sequence[str], bool]] = None,
                  homogenizations: Optional[Union[str, Sequence[str], bool]] = None,
                  fields: Optional[Union[str, Sequence[str], bool]] = None) -> "Result":
        """
        Remove from view.

        Wildcard matching with '?' and '*' is supported.
        True is equivalent to '*', False is equivalent to [].

        Parameters
        ----------
        increments : (list of) int, (list of) str, or bool, optional
            Numbers of increments to select. Mutually exclusive with 'times'.
        times : (list of) float, '*', or bool, optional
            Simulation times of increments to select. Mutually exclusive with
            'increments'.
        phases : (list of) str, or bool, optional
            Names of phases to select.
        homogenizations : (list of) str, or bool, optional
            Names of homogenizations to select.
        fields : (list of) str, or bool, optional
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


    def view_all(self):
        """
        Make all attributes visible.

        Returns
        -------
        modified_view : damask.Result
            View with all attributes visible.
        """
        return self.view(increments='*',phases='*',homogenizations='*',fields='*')


    def rename(self,
               name_src: str,
               name_dst: str):
        r"""
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
        \x1b[93m\x1b[1mWarning: Modification of existing datasets allowed!\x1b[0m\x1b[0m
        >>> r_unprotected.rename('F','def_grad')
        """
        if self._protected:
            raise PermissionError('rename datasets')

        with h5py.File(self.fname,'a') as f:
            for inc in self._visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self._visible[ty+'s']:
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            path_src = '/'.join([inc,ty,label,field,name_src])
                            path_dst = '/'.join([inc,ty,label,field,name_dst])
                            if path_src in f.keys():
                                f[path_dst] = f[path_src]
                                f[path_dst].attrs['renamed'] = f'original name: {name_src}'
                                del f[path_src]


    def remove(self, name: str):
        r"""
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
        \x1b[93m\x1b[1mWarning: Modification of existing datasets allowed!\x1b[0m\x1b[0m
        >>> r_unprotected.remove('F')
        """
        if self._protected:
            raise PermissionError('delete datasets')

        with h5py.File(self.fname,'a') as f:
            for inc in self._visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self._visible[ty+'s']:
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            path = '/'.join([inc,ty,label,field,name])
                            if path in f.keys(): del f[path]


    def list_data(self) -> list[str]:
        """
        Collect information on all active datasets in the file.

        Returns
        -------
        data : list of str
            Line-formatted information about active datasets.
        """
        msg = []
        with h5py.File(self.fname,'r') as f:
            for inc in self._visible['increments']:
                msg.append(f'\n{inc} ({self._times[int(inc.split("_")[1])]} s)')
                for ty in ['phase','homogenization']:
                    msg.append(f'  {ty}')
                    for label in self._visible[ty+'s']:
                        msg.append(f'    {label}')
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            msg.append(f'      {field}')
                            for d in f['/'.join([inc,ty,label,field])].keys():
                                dataset = f['/'.join([inc,ty,label,field,d])]
                                unit = dataset.attrs["unit"]
                                description = dataset.attrs['description']
                                msg.append(f'        {d} / {unit}: {description}')

        return msg


    def enable_user_function(self,
                             func: Callable):
        globals()[func.__name__]=func
        logger.info(f'Function {func.__name__} enabled in add_calculation.')


    @property
    def simulation_setup_files(self):
        """Simulation setup files used to generate the Result object."""
        files = []
        with h5py.File(self.fname,'r') as f_in:
            f_in['setup'].visititems(lambda name,obj: files.append(name) if isinstance(obj,h5py.Dataset) else None)
        return files

    @property
    def _incs(self):
        return [int(i.split(prefix_inc)[-1]) for i in self._increments]

    @property
    def increments(self):
        return [i for i in self._visible['increments']]

    @property
    def times(self):
        return [self._times[int(i.split(prefix_inc)[-1])] for i in self.increments]

    @property
    def phases(self):
        return self._visible['phases']

    @property
    def homogenizations(self):
        return self._visible['homogenizations']

    @property
    def fields(self):
        return self._visible['fields']


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
                                                  f['/geometry/T_c'].attrs['VTK_TYPE'])


    def add_absolute(self, x: str):
        """
        Add absolute value.

        Parameters
        ----------
        x : str
            Name of scalar, vector, or tensor dataset to take absolute value of.

        Notes
        -----
        For details refer to ``numpy.abs``.
        """
        def absolute(x: DADF5Dataset) -> DADF5Dataset:
            return {
                    'data':  np.abs(x['data']),
                    'label': f'|{x["label"]}|',
                    'meta':  {
                              'unit':        x['meta']['unit'],
                              'description': f"absolute value of {x['label']} ({x['meta']['description']})",
                              'creator':     'add_absolute'
                              }
                     }

        self._add_generic_pointwise(absolute,{'x':x})


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
        [...]
        >>> r.add_calculation('np.sum(#rho_dip#,axis=1)','rho_dip_total',
        ...                    '1/m²','total dislocation dipole density')
        [...]
        >>> r.add_calculation('#rho_dip_total#+#rho_mob_total#','rho_total',
        ...                    '1/m²','total dislocation density')
        [...]

        Add von Mises equivalent of the Cauchy stress without storage of
        intermediate results. Define a user function for better readability:

        >>> import damask
        >>> def equivalent_stress(F,P):
        ...     sigma = damask.mechanics.stress_Cauchy(F=F,P=P)
        ...     return damask.mechanics.equivalent_stress_Mises(sigma)
        >>> r = damask.Result('my_file.hdf5')
        >>> r.enable_user_function(equivalent_stress)
        Function equivalent_stress enabled in add_calculation.
        >>> r.add_calculation('equivalent_stress(#F#,#P#)','sigma_vM','Pa',
        ...                   'von Mises equivalent of the Cauchy stress')
        [...]
        """
        def calculation(**kwargs) -> DADF5Dataset:
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

        dataset_mapping = {d:d for d in set(re.findall(r'#(.*?)#',formula))}                        # datasets used in the formula
        args             = {'formula':formula,'label':name,'unit':unit,'description':description}
        self._add_generic_pointwise(calculation,dataset_mapping,args)


    def add_determinant(self, T: str):
        """
        Add the determinant of a tensor.

        Parameters
        ----------
        T : str
            Name of tensor dataset.

        Notes
        -----
        For details refer to ``numpy.linalg.det``.

        Examples
        --------
        Add the determinant of plastic deformation gradient 'F_p':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_determinant('F_p')
        [...]
        """
        def determinant(T: DADF5Dataset) -> DADF5Dataset:
            return {
                    'data':  np.linalg.det(T['data']),
                    'label': f"det({T['label']})",
                    'meta':  {
                              'unit':        T['meta']['unit'],
                              'description': f"determinant of tensor {T['label']} ({T['meta']['description']})",
                              'creator':     'add_determinant'
                              }
                    }

        self._add_generic_pointwise(determinant,{'T':T})


    def add_deviator(self, T: str):
        """
        Add the deviatoric part of a tensor.

        Parameters
        ----------
        T : str
            Name of tensor dataset.

        Notes
        -----
        For details refer to :func:`damask.tensor.deviatoric`.

        Examples
        --------
        Add the deviatoric part of Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_deviator('sigma')
        [...]
        """
        def deviator(T: DADF5Dataset) -> DADF5Dataset:
            return {
                    'data':  tensor.deviatoric(T['data']),
                    'label': f"s_{T['label']}",
                    'meta':  {
                              'unit':        T['meta']['unit'],
                              'description': f"deviator of tensor {T['label']} ({T['meta']['description']})",
                              'creator':     'add_deviator'
                              }
                     }

        self._add_generic_pointwise(deviator,{'T':T})


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

        Notes
        -----
        For details refer to :func:`damask.tensor.eigenvalues`.

        Examples
        --------
        Add the minimum eigenvalue of Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_eigenvalue('sigma','min')
        [...]
        """
        def eigenval(T_sym: DADF5Dataset, eigenvalue: Literal['max', 'mid', 'min']) -> DADF5Dataset:
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

        self._add_generic_pointwise(eigenval,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


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

        Notes
        -----
        For details refer to :func:`damask.tensor.eigenvectors`.
        """
        def eigenvector(T_sym: DADF5Dataset, eigenvalue: Literal['max', 'mid', 'min']) -> DADF5Dataset:
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

        self._add_generic_pointwise(eigenvector,{'T_sym':T_sym},{'eigenvalue':eigenvalue})


    def add_equivalent_Mises(self,
                             T_sym: str,
                             kind: Optional[str] = None):
        """
        Add the equivalent von Mises stress or strain of a symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Name of symmetric tensorial stress or strain dataset.
        kind : {'stress', 'strain', None}, optional
            Kind of the von Mises equivalent. Defaults to None, in which case
            it is selected based on the unit of the dataset ('1' -> strain, 'Pa' -> stress).

        Notes
        -----
        For details refer to :func:`damask.mechanics.equivalent_stress_Mises` and
        :func:`damask.mechanics.equivalent_strain_Mises`.

        Examples
        --------
        Add the von Mises equivalent of the Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_equivalent_Mises('sigma')
        [...]

        Add the von Mises equivalent of the spatial logarithmic strain 'epsilon_V^0.0(F)':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_equivalent_Mises('epsilon_V^0.0(F)')
        [...]
        """
        def equivalent_Mises(T_sym: DADF5Dataset, kind: str) -> DADF5Dataset:
            k = kind
            if k is None:
                if T_sym['meta']['unit'] == '1':
                    k = 'strain'
                elif T_sym['meta']['unit'] == 'Pa':
                    k = 'stress'
            if k not in ['stress', 'strain']:
                raise ValueError(f'invalid von Mises kind "{kind}"')

            return {
                    'data':  mechanics._equivalent_Mises(T_sym['data'],
                                                        2./3. if k == 'strain' else 3./2.),
                    'label': f"{T_sym['label']}_vM",
                    'meta':  {
                              'unit':        T_sym['meta']['unit'],
                              'description': f"Mises equivalent {k} of {T_sym['label']} ({T_sym['meta']['description']})",
                              'creator':     'add_Mises'
                              }
                    }

        self._add_generic_pointwise(equivalent_Mises,{'T_sym':T_sym},{'kind':kind})


    def add_IPF_color(self,
                      l: FloatSequence,
                      q: str = 'O'):
        """
        Add RGB color tuple of inverse pole figure (IPF) color.

        Parameters
        ----------
        l : numpy.array of shape (3) or compatible
            Lab frame direction for inverse pole figure.
        q : str, optional
            Name of the dataset containing the crystallographic orientation as quaternions.
            Defaults to 'O'.

        Notes
        -----
        For details refer to :func:`damask.Orientation.IPF_color`.

        Examples
        --------
        Add the IPF color along x-direction for orientation 'O':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_IPF_color(l = [1,0,0], q = 'O')
        [...]
        """
        def IPF_color(l: FloatSequence, q: DADF5Dataset) -> DADF5Dataset:
            m = util.scale_to_coprime(np.array(l))
            lattice: BravaisLattice = q['meta']['lattice']                                          # type: ignore[assignment]
            o = Orientation(rotation = q['data'],lattice=lattice)

            return {
                    'data': (o.IPF_color(l)*255).astype(np.uint8),
                    'label': 'IPFcolor_({} {} {})'.format(*m),
                    'meta' : {
                              'unit':        '8-bit RGB',
                              'lattice':     lattice,
                              'description': 'Inverse Pole Figure (IPF) colors along sample direction ({} {} {})'.format(*m),
                              'creator':     'add_IPF_color'
                             }
                   }

        self._add_generic_pointwise(IPF_color,{'q':q},{'l':l})


    def add_maximum_shear(self, T_sym: str):
        """
        Add maximum shear components of symmetric tensor.

        Parameters
        ----------
        T_sym : str
            Name of symmetric tensor dataset.

        Notes
        -----
        For details refer to :func:`damask.mechanics.maximum_shear`.
        """
        def maximum_shear(T_sym: DADF5Dataset) -> DADF5Dataset:
            return {
                    'data':  mechanics.maximum_shear(T_sym['data']),
                    'label': f"max_shear({T_sym['label']})",
                    'meta':  {
                              'unit':        T_sym['meta']['unit'],
                              'description': f"maximum shear component of {T_sym['label']} ({T_sym['meta']['description']})",
                              'creator':     'add_maximum_shear'
                              }
                     }

        self._add_generic_pointwise(maximum_shear,{'T_sym':T_sym})


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
            Order of the norm. inf means NumPy's inf object.

        Notes
        -----
        For details refer to numpy.linalg.norm.
        """
        def norm(x: DADF5Dataset, ord: Union[int, float, Literal['fro', 'nuc']]) -> DADF5Dataset:
            o = ord
            if len(x['data'].shape) == 2:
                axis: Union[int, tuple[int, int]] = 1
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

        self._add_generic_pointwise(norm,{'x':x},{'ord':ord})


    def add_pole(self,
                 q: str = 'O',
                 *,
                 uvw: Optional[IntSequence] = None,
                 hkl: Optional[IntSequence] = None,
                 with_symmetry: bool = False,
                 normalize: bool = True):
        """
        Add lab frame vector along lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        q : str, optional
            Name of the dataset containing the crystallographic orientation as quaternions.
            Defaults to 'O'.
        uvw|hkl : numpy.ndarray of shape (3) or compatible
            Miller indices of crystallographic direction or plane normal.
        with_symmetry : bool, optional
            Calculate all N symmetrically equivalent vectors.
            Defaults to False.
        normalize : bool, optional
            Normalize output vector.
            Defaults to True.

        Notes
        -----
        For details refer to :func:`damask.Orientation.to_frame`.
        """
        def pole(q: DADF5Dataset,
                 uvw: IntSequence, hkl: IntSequence,
                 with_symmetry: bool,
                 normalize: bool) -> DADF5Dataset:
            c = q['meta']['c/a'] if 'c/a' in q['meta'] else 1.0
            brackets = ['[]','()','⟨⟩','{}'][(uvw is None)*1+with_symmetry*2]
            label = 'p^' + '{}{} {} {}{}'.format(brackets[0],
                                                 *(uvw if uvw else hkl),
                                                 brackets[-1],)
            lattice: BravaisLattice = q['meta']['lattice']                                          # type: ignore[assignment]
            ori = Orientation(q['data'],lattice=lattice,a=1,c=c)

            return {
                    'data': np.moveaxis(ori.to_frame(uvw=uvw,hkl=hkl,
                                                     with_symmetry=with_symmetry,
                                                     normalize=normalize),0,-2 if with_symmetry else 0),
                    'label': label,
                    'meta' : {
                              'unit':        '1',
                              'description': f'{"normalized " if normalize else ""}lab frame vector along lattice ' \
                                             + ('plane' if uvw is None else 'direction') \
                                             + ('s' if with_symmetry else ''),
                              'creator':     'add_pole'
                              }
                    }

        self._add_generic_pointwise(pole,{'q':q},{'uvw':uvw,'hkl':hkl,'with_symmetry':with_symmetry,'normalize':normalize})


    def add_rotation(self, F: str):
        """
        Add rotational part of a deformation gradient.

        Parameters
        ----------
        F : str
            Name of deformation gradient dataset.

        Notes
        -----
        For details refer to :func:`damask.mechanics.rotation`.

        Examples
        --------
        Add the rotational part of deformation gradient 'F':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_rotation('F')
        [...]
        """
        def rotation(F: DADF5Dataset) -> DADF5Dataset:
            return {
                    'data':  mechanics.rotation(F['data']).as_matrix(),
                    'label': f"R({F['label']})",
                    'meta':  {
                              'unit':        F['meta']['unit'],
                              'description': f"rotational part of {F['label']} ({F['meta']['description']})",
                              'creator':     'add_rotation'
                              }
                     }

        self._add_generic_pointwise(rotation,{'F':F})


    def add_spherical(self, T: str):
        """
        Add the spherical (hydrostatic) part of a tensor.

        Parameters
        ----------
        T : str
            Name of tensor dataset.

        Notes
        -----
        For details refer to :func:`damask.tensor.spherical`.

        Examples
        --------
        Add the hydrostatic part of the Cauchy stress 'sigma':

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_spherical('sigma')
        [...]
        """
        def spherical(T: DADF5Dataset) -> DADF5Dataset:
            return {
                    'data':  tensor.spherical(T['data'],False),
                    'label': f"p_{T['label']}",
                    'meta':  {
                              'unit':        T['meta']['unit'],
                              'description': f"spherical component of tensor {T['label']} ({T['meta']['description']})",
                              'creator':     'add_spherical'
                              }
                     }

        self._add_generic_pointwise(spherical,{'T':T})


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

        Notes
        -----
        For details refer to :func:`damask.mechanics.strain`.

        Examples
        --------
        Add the Euler-Almansi strain:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_strain(t='V',m=-1.0)
        [...]

        Add the plastic Biot strain:

        >>> import damask
        >>> r = damask.Result('my_file.hdf5')
        >>> r.add_strain('F_p','U',0.5)
        [...]
        """
        def strain(F: DADF5Dataset, t: Literal['V', 'U'], m: float) -> DADF5Dataset:
            side = 'left' if t == 'V' else 'right'
            return {
                    'data':  mechanics.strain(F['data'],t,m),
                    'label': f"epsilon_{t}^{m}({F['label']})",
                    'meta':  {
                              'unit':        F['meta']['unit'],
                              'description': f'Seth-Hill strain tensor of order {m} based on {side} stretch tensor '
                                             f"of {F['label']} ({F['meta']['description']})",
                              'creator':     'add_strain'
                              }
                     }

        self._add_generic_pointwise(strain,{'F':F},{'t':t,'m':m})


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

        Notes
        -----
        For details refer to :func:`damask.mechanics.stress_Cauchy`.
        """
        def stress_Cauchy(P: DADF5Dataset, F: DADF5Dataset) -> DADF5Dataset:
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

        self._add_generic_pointwise(stress_Cauchy,{'P':P,'F':F})


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
        For details refer to :func:`damask.mechanics.stress_second_Piola_Kirchhoff`.
        """
        def stress_second_Piola_Kirchhoff(P: DADF5Dataset, F: DADF5Dataset) -> DADF5Dataset:
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

        self._add_generic_pointwise(stress_second_Piola_Kirchhoff,{'P':P,'F':F})


    def add_stretch_tensor(self,
                           F: str = 'F',
                           t: Literal['V', 'U'] = 'V'):
        """
        Add stretch tensor of a deformation gradient.

        Notes
        -----
        For details refer to :func:`damask.mechanics.stretch_left` and
        :func:`damask.mechanics.stretch_right`.

        Parameters
        ----------
        F : str, optional
            Name of deformation gradient dataset. Defaults to 'F'.
        t : {'V', 'U'}, optional
            Type of the polar decomposition, 'V' for left stretch tensor and 'U' for right stretch tensor.
            Defaults to 'V'.
        """
        def stretch_tensor(F: DADF5Dataset, t: str) -> DADF5Dataset:
            return {
                    'data':  mechanics._polar_decomposition(F['data'],t)[0],
                    'label': f"{t}({F['label']})",
                    'meta':  {
                              'unit':        F['meta']['unit'],
                              'description': f"{'left' if t.upper() == 'V' else 'right'} stretch tensor "
                                             f"of {F['label']} ({F['meta']['description']})",           # noqa
                              'creator':     'add_stretch_tensor'
                              }
                     }

        self._add_generic_pointwise(stretch_tensor,{'F':F},{'t':t})


    def add_curl(self, f: str):
        """
        Add curl of a field.

        Parameters
        ----------
        f : str
            Name of vector or tensor field dataset.

        Notes
        -----
        For details refer to :func:`damask.grid_filters.curl`.

        This function is implemented only for structured grids
        with one constituent and a single phase.
        """
        def curl(f: DADF5Dataset, size: np.ndarray) -> DADF5Dataset:
            return {
                    'data':  grid_filters.curl(size,f['data']),
                    'label': f"curl({f['label']})",
                    'meta':  {
                              'unit':        f['meta']['unit']+'/m',
                              'description': f"curl of {f['label']} ({f['meta']['description']})",
                              'creator':     'add_curl'
                              }
                     }

        self._add_generic_grid(curl,{'f':f},{'size':self.size})


    def add_divergence(self, f: str):
        """
        Add divergence of a field.

        Parameters
        ----------
        f : str
            Name of vector or tensor field dataset.

        Notes
        -----
        For details refer to :func: `damask.grid_filters.divergence`.

        This function is implemented only for structured grids
        with one constituent and a single phase.
        """
        def divergence(f: DADF5Dataset, size: np.ndarray) -> DADF5Dataset:
            return {
                    'data':  grid_filters.divergence(size,f['data']),
                    'label': f"divergence({f['label']})",
                    'meta':  {
                              'unit':        f['meta']['unit']+'/m',
                              'description': f"divergence of {f['label']} ({f['meta']['description']})",
                              'creator':     'add_divergence'
                              }
                     }

        self._add_generic_grid(divergence,{'f':f},{'size':self.size})


    def add_gradient(self, f: str):
        """
        Add gradient of a field.

        Parameters
        ----------
        f : str
            Name of scalar or vector field dataset.

        Notes
        -----
        For details refer to :func:`damask.grid_filters.gradient`.

        This function is implemented only for structured grids
        with one constituent and a single phase.
        """
        def gradient(f: DADF5Dataset, size: np.ndarray) -> DADF5Dataset:
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

        self._add_generic_grid(gradient,{'f':f},{'size':self.size})


    def _add_generic_grid(self,
                          func: Callable[..., DADF5Dataset],
                          datasets: dict[str, str],
                          args: dict[str,Any]):
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
        args : dictionary
            Arguments parsed to func.
        """
        if self.N_constituents != 1 or len(datasets) != 1 or not self.structured:
            raise NotImplementedError('not a structured grid with one constituent and a single phase')

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        increments = self.place(list(datasets.values()),False)
        if not increments: raise RuntimeError('received invalid dataset')
        with h5py.File(self.fname, 'a') as f:
            for increment in increments.items():
                for ty in increment[1].items():
                    for field in ty[1].items():
                        d: np.ma.MaskedArray = list(field[1].values())[0]
                        if np.any(d.mask): continue

                        dataset = {'f':{'data':grid_filters.unravel(d.data,self.cells),
                                        'label':list(datasets.values())[0],
                                        'meta':d.data.dtype.metadata}}
                        r = func(**dataset,**args)
                        result = grid_filters.ravel(r['data'])
                        for x in self._visible[ty[0]+'s']:
                            path = '/'.join(['/',increment[0],ty[0],x,field[0]])
                            if ty[0] == 'phase':
                                result1 = result[at_cell_ph[0][x]]
                            if ty[0] == 'homogenization':
                                result1 = result[at_cell_ho[x]]
                            if not self._protected and '/'.join([path,r['label']]) in f:
                                h5_dataset = f['/'.join([path,r['label']])]
                                h5_dataset[...] = result1
                                h5_dataset.attrs['overwritten'] = True
                            else:
                                h5_dataset = f[path].create_dataset(r['label'],data=result1)

                            h5_dataset.attrs['created'] = util.time_stamp()

                            for l,v in r['meta'].items():
                                h5_dataset.attrs[l.lower()] = v
                            creator = h5_dataset.attrs['creator']
                            h5_dataset.attrs['creator'] = f'damask.Result.{creator} v{damask.version}'


    def _add_generic_pointwise(self,
                               func: Callable[..., DADF5Dataset],
                               datasets: dict[str, str],
                               args: Optional[dict[str, Any]] = None):
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
        args = args if args else {}

        def job_pointwise(group: str,
                          callback: Callable[..., DADF5Dataset],
                          datasets: dict[str, str],
                          args: dict[str, str]) -> Union[None, DADF5Dataset]:
            try:
                datasets_in = {}
                with h5py.File(self.fname,'r') as f:
                    for arg,label in datasets.items():
                        loc  = f[group+'/'+label]
                        datasets_in[arg]={'data' :loc[()],
                                          'label':label,
                                          'meta': {k: v for k,v in loc.attrs.items()}}
                return callback(**datasets_in,**args)
            except Exception as err:
                logger.error(f'Error during pointwise calculation: {err}.')
                return None

        groups = []
        with h5py.File(self.fname,'r') as f:
            for inc in self._visible['increments']:
                for ty in ['phase','homogenization']:
                    for label in self._visible[ty+'s']:
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            group = '/'.join([inc,ty,label,field])
                            if set(datasets.values()).issubset(f[group].keys()): groups.append(group)

        if len(groups) == 0:
            logger.warning('No matching dataset found, no data was added.')
            return


        for group in util.show_progress(groups):
            if not (result := job_pointwise(group, callback=func, datasets=datasets, args=args)):   # type: ignore
                continue
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

                    dataset.attrs['created'] = util.time_stamp()

                    for l,v in result['meta'].items():
                        dataset.attrs[l.lower()] = v
                    creator = dataset.attrs['creator']
                    dataset.attrs['creator'] = f'damask.Result.{creator} v{damask.version}'

                except (OSError,RuntimeError) as err:
                    logger.error(f'Could not add dataset: {err}.')


    def _mappings(self) -> MappingsTuple:
        """Mappings to place data spatially."""
        with h5py.File(self.fname,'r') as f:

            at_cell_ph = []
            in_data_ph = []
            for c in range(self.N_constituents):
                at_cell_ph.append({label: np.where(self.phase[:,c] == label)[0] \
                                          for label in self._visible['phases']})
                in_data_ph.append({label: f['/'.join(['cell_to','phase'])]['entry'][at_cell_ph[c][label]][:,c] \
                                          for label in self._visible['phases']})

            at_cell_ho = {label: np.where(self.homogenization[:] == label)[0] \
                                 for label in self._visible['homogenizations']}
            in_data_ho = {label: f['/'.join(['cell_to','homogenization'])]['entry'][at_cell_ho[label]] \
                                 for label in self._visible['homogenizations']}

        return MappingsTuple(at_cell_ph, in_data_ph, at_cell_ho, in_data_ho)


    def get(self,
            output: Union[str, list[str]] = '*',
            flatten: bool = True,
            prune: bool = True) -> Union[None,dict[str,Any]]:
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
        r: dict[str,Any] = {}

        with h5py.File(self.fname,'r') as f:
            for inc in util.show_progress(self._visible['increments']):
                r[inc] = {'phase':{},'homogenization':{},'geometry':{}}

                for out in _match(output,f['/'.join([inc,'geometry'])].keys()):
                    r[inc]['geometry'][out] = _read(f['/'.join([inc,'geometry',out])])

                for ty in ['phase','homogenization']:
                    for label in self._visible[ty+'s']:
                        r[inc][ty][label] = {}
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            r[inc][ty][label][field] = {}
                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                r[inc][ty][label][field][out] = _read(f['/'.join([inc,ty,label,field,out])])

        if prune:   r = util.dict_prune(r)
        if flatten: r = util.dict_flatten(r)

        return None if (type(r) is dict and r == {}) else r


    def place(self,
              output: Union[str, list[str]] = '*',
              flatten: bool = True,
              prune: bool = True,
              constituents: Optional[IntSequence] = None,
              fill_float: float = np.nan,
              fill_int: int = 0) -> Optional[dict[str,Any]]:
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
        r: dict[str,Any] = {}

        constituents_ = map(int,constituents) if isinstance(constituents,Iterable) else \
                      (range(self.N_constituents) if constituents is None else [constituents])      # type: ignore

        suffixes = [''] if self.N_constituents == 1 or isinstance(constituents,int) else \
                   [f'#{c}' for c in constituents_]                                                 # type: ignore

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        with h5py.File(self.fname,'r') as f:

            for inc in util.show_progress(self._visible['increments']):
                r[inc] = {'phase':{},'homogenization':{},'geometry':{}}

                for out in _match(output,f['/'.join([inc,'geometry'])].keys()):
                    r[inc]['geometry'][out] = ma.array(_read(f['/'.join([inc,'geometry',out])]),fill_value = fill_float)

                for ty in ['phase','homogenization']:

                    dtypes_by_out: dict[str, Any] = {}
                    for label in self._visible[ty + 's']:
                        for field in _match(self._visible['fields'], f['/'.join([inc, ty, label])].keys()):
                            for out in _match(output, f['/'.join([inc, ty, label, field])].keys()):
                                dtypes_by_out.setdefault(out, []).append(_read_dt(f['/'.join([inc, ty, label, field, out])]))
                    dtype_by_out = {out: _get_common_metadata(dtypes) for out, dtypes in dtypes_by_out.items()}

                    for label in self._visible[ty+'s']:
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            if field not in r[inc][ty].keys():
                                r[inc][ty][field] = {}

                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                data = ma.array(_read(f['/'.join([inc,ty,label,field,out])]))

                                if ty == 'phase':
                                    if out+suffixes[0] not in r[inc][ty][field].keys():
                                        for c,suffix in zip(constituents_,suffixes):
                                            r[inc][ty][field][out+suffix] = \
                                                _empty_like(data.shape,dtype_by_out[out],
                                                            self.N_materialpoints,fill_float,fill_int)

                                    for c,suffix in zip(constituents_,suffixes):
                                        r[inc][ty][field][out+suffix][at_cell_ph[c][label]] = data[in_data_ph[c][label]] # type: ignore

                                if ty == 'homogenization':
                                    if out not in r[inc][ty][field].keys():
                                        r[inc][ty][field][out] = \
                                            _empty_like(data.shape,dtype_by_out[out],
                                                        self.N_materialpoints,fill_float,fill_int)

                                    r[inc][ty][field][out][at_cell_ho[label]] = data[in_data_ho[label]]

        if prune:   r = util.dict_prune(r)
        if flatten: r = util.dict_flatten(r)

        return None if (type(r) is dict and r == {}) else r


    def export_XDMF(self,
                    output: Union[str, list[str]] = '*',
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

        Notes
        -----
        This function is implemented only for structured grids with
        one constituent and a single phase.
        """
        if self.N_constituents != 1 or len(self.phases) != 1 or not self.structured:
            raise NotImplementedError('not a structured grid with one constituent and a single phase')

        attribute_type_map = defaultdict(lambda:'Matrix', ( ((),'Scalar'), ((3,),'Vector'), ((3,3),'Tensor')) )

        def number_type_map(dtype):
            if np.issubdtype(dtype,np.signedinteger):   return 'Int'
            if np.issubdtype(dtype,np.unsignedinteger): return 'UInt'
            if np.issubdtype(dtype,np.floating):        return 'Float'
            raise TypeError(f'invalid type "{dtype}"')


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
        time_data.attrib = {'Format':     'XML',
                            'NumberType': 'Float',
                            'Dimensions': f'{len(self.times)}'}
        time_data.text = ' '.join(map(str,self.times))

        attributes = []
        data_items = []

        hdf5_name = self.fname.name
        hdf5_dir  = self.fname.parent
        out_dir   = Path.cwd() if target_dir is None else Path(target_dir)
        hdf5_link = (hdf5_dir if absolute_path else Path(os.path.relpath(hdf5_dir,out_dir.resolve())))/hdf5_name

        with h5py.File(self.fname,'r') as f:
            for inc in self._visible['increments']:

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
                    for label in self._visible[ty+'s']:
                        for field in _match(self._visible['fields'],f['/'.join([inc,ty,label])].keys()):
                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                name = '/'.join([inc,ty,label,field,out])
                                shape = f[name].shape[1:]
                                dtype = f[name].dtype

                                unit = f[name].attrs['unit']

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

        out_dir.mkdir(parents=True,exist_ok=True)
        with util.open_text((out_dir/hdf5_name).with_suffix('.xdmf'),'w') as f:
            f.write(xml.dom.minidom.parseString(ET.tostring(xdmf).decode()).toprettyxml())


    def export_VTK(self,
                   output: Union[str,list[str]] = '*',
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
        For cell data, the file format is either ImageData (.vti) or
        UnstructuredGrid (.vtu) for grid-based or mesh-based simulations,
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

        N_digits = int(np.floor(np.log10(max(1,self._incs[-1]))))+1

        constituents_ = constituents if isinstance(constituents,Iterable) else \
                        (range(self.N_constituents) if constituents is None else [constituents])    # type: ignore

        suffixes = [''] if self.N_constituents == 1 or isinstance(constituents,int) else \
                   [f'#{c}' for c in constituents_]

        at_cell_ph,in_data_ph,at_cell_ho,in_data_ho = self._mappings()

        out_dir = Path.cwd() if target_dir is None else Path(target_dir)
        out_dir.mkdir(parents=True,exist_ok=True)

        with h5py.File(self.fname,'r') as f:
            creator = f.attrs['creator']
            created = f.attrs['created']
            v.comments.append(f'{creator} ({created})')

            for inc in util.show_progress(self._visible['increments']):

                u = _read(f['/'.join([inc,'geometry','u_n' if mode.lower() == 'cell' else 'u_p'])])
                v = v.set('u',u)

                for ty in ['phase','homogenization']:

                    dtypes_by_out: dict[str, Any] = {}
                    for label in self._visible[ty + 's']:
                        for field in _match(self._visible['fields'], f['/'.join([inc, ty, label])].keys()):
                            for out in _match(output, f['/'.join([inc, ty, label, field])].keys()):
                                dtypes_by_out.setdefault(out, []).append(_read_dt(f['/'.join([inc, ty, label, field, out])]))
                    dtype_by_out = {out: _get_common_metadata(dtypes) for out, dtypes in dtypes_by_out.items()}

                    for field in self._visible['fields']:
                        outs: dict[str, np.ma.core.MaskedArray] = {}
                        for label in self._visible[ty+'s']:
                            if field not in f['/'.join([inc,ty,label])].keys(): continue

                            for out in _match(output,f['/'.join([inc,ty,label,field])].keys()):
                                data = ma.array(_read(f['/'.join([inc,ty,label,field,out])]))

                                if ty == 'phase':
                                    if out+suffixes[0] not in outs.keys():
                                        for c,suffix in zip(constituents_,suffixes):
                                            outs[out+suffix] = \
                                                _empty_like(data.shape,dtype_by_out[out],
                                                            self.N_materialpoints,fill_float,fill_int)

                                    for c,suffix in zip(constituents_,suffixes):
                                        outs[out+suffix][at_cell_ph[c][label]] = data[in_data_ph[c][label]]

                                if ty == 'homogenization':
                                    if out not in outs.keys():
                                        outs[out] = _empty_like(data.shape,dtype_by_out[out],
                                                                self.N_materialpoints,fill_float,fill_int)
                                    outs[out][at_cell_ho[label]] = data[in_data_ho[label]]

                        for label,dataset in outs.items():
                            v = v.set(' / '.join(['/'.join([ty,field,label]),dataset.dtype.metadata.get('unit')]),dataset)

                v.save(out_dir/f'{self.fname.stem}_inc{inc.split(prefix_inc)[-1].zfill(N_digits)}',
                       parallel=parallel)


    def export_DREAM3D(self,
                       q: str = 'O',
                       target_dir: Union[None, str, Path] = None):
        """
        Export the visible components to DREAM3D compatible files.

        One DREAM3D file per visible increment is created.
        The geometry is based on the undeformed configuration.

        Parameters
        ----------
        q : str, optional
            Name of the dataset containing the crystallographic orientation as quaternions.
            Defaults to 'O'.

        target_dir : str or pathlib.Path, optional
            Directory to save DREAM3D files. Will be created if non-existent.

        Notes
        -----
        This function is implemented only for structured grids with
        one constituent.
        """
        def add_attribute(obj,name,data):
            """DREAM.3D requires fixed length string."""
            if isinstance(data,str):
                tid = h5py.h5t.C_S1.copy()
                tid.set_size(len(data)+1)
                obj.attrs.create(name,data,dtype=h5py.Datatype(tid))
            else:
                obj.attrs.create(name,data)

        def create_and_open(obj,name):
            obj.create_group(name)
            return obj[name]

        if self.N_constituents != 1 or not self.structured:
            raise NotImplementedError('not a structured grid with one constituent')

        N_digits = int(np.floor(np.log10(max(1,self._incs[-1]))))+1

        at_cell_ph,in_data_ph,_,_ = self._mappings()

        out_dir = Path.cwd() if target_dir is None else Path(target_dir)
        out_dir.mkdir(parents=True,exist_ok=True)

        with h5py.File(self.fname,'r') as f:
            for inc in util.show_progress(self._visible['increments']):
                for c in range(self.N_constituents):
                    crystal_structure = [999]
                    phase_name = ['Unknown Phase Type']
                    cell_orientation = np.zeros((np.prod(self.cells),3),np.float32)
                    phase_ID = np.zeros((np.prod(self.cells)),dtype=np.int32)
                    count = 1
                    for label in self._visible['phases']:
                        try:
                            data = _read(f['/'.join([inc,'phase',label,'mechanical',q])])
                            lattice = data.dtype.metadata['lattice']
                            # Map to DREAM.3D IDs
                            if lattice == 'hP':
                                crystal_structure.append(0)
                            elif lattice in ['cI','cF']:
                                crystal_structure.append(1)
                            elif lattice == 'tI':
                                crystal_structure.append(8)

                            cell_orientation[at_cell_ph[c][label],:] = \
                                Rotation(data[in_data_ph[c][label],:]).as_Euler_angles().astype(np.float32)
                            phase_ID[at_cell_ph[c][label]] = count
                            phase_name.append(label)
                            count +=1
                        except KeyError:
                            pass


                with h5py.File(f'{out_dir}/{self.fname.stem}_inc{inc.split(prefix_inc)[-1].zfill(N_digits)}.dream3d','w') as f_out:
                    add_attribute(f_out,'FileVersion','7.0')

                    for g in ['DataContainerBundles','Pipeline']:                                   # empty groups (needed)
                        f_out.create_group(g)

                    data_container = create_and_open(f_out,'DataContainers/SyntheticVolumeDataContainer')

                    cell = create_and_open(data_container,'CellData')
                    add_attribute(cell,'AttributeMatrixType',np.array([3],np.uint32))
                    add_attribute(cell,'TupleDimensions', np.array(self.cells,np.uint64))

                    cell['Phases'] = np.reshape(phase_ID,tuple(np.flip(self.cells))+(1,))
                    cell['EulerAngles'] = cell_orientation.reshape(tuple(np.flip(self.cells))+(3,))
                    for dataset in ['Phases','EulerAngles']:
                        add_attribute(cell[dataset],'DataArrayVersion',np.array([2],np.int32))
                        add_attribute(cell[dataset],'Tuple Axis Dimensions','x={},y={},z={}'.format(*np.array(self.cells)))
                        add_attribute(cell[dataset],'TupleDimensions', np.array(self.cells,np.uint64))
                    add_attribute(cell['Phases'], 'ComponentDimensions', np.array([1],np.uint64))
                    add_attribute(cell['Phases'], 'ObjectType', 'DataArray<int32_t>')
                    add_attribute(cell['EulerAngles'], 'ComponentDimensions', np.array([3],np.uint64))
                    add_attribute(cell['EulerAngles'], 'ObjectType', 'DataArray<float>')

                    cell_ensemble =  create_and_open(data_container,'CellEnsembleData')

                    cell_ensemble['CrystalStructures'] = np.array(crystal_structure,np.uint32).reshape(-1,1)
                    cell_ensemble['PhaseTypes'] = np.array([999] + [0]*(len(crystal_structure)-1),np.uint32).reshape(-1,1)
                    tid = h5py.h5t.C_S1.copy()
                    tid.set_size(h5py.h5t.VARIABLE)
                    tid.set_cset(h5py.h5t.CSET_ASCII)
                    cell_ensemble.create_dataset(name='PhaseName',data = phase_name, dtype=h5py.Datatype(tid))

                    cell_ensemble.attrs['AttributeMatrixType'] = np.array([11],np.uint32)
                    cell_ensemble.attrs['TupleDimensions']     = np.array([len(self._phases) + 1], np.uint64)
                    for group in ['CrystalStructures','PhaseTypes','PhaseName']:
                        add_attribute(cell_ensemble[group], 'ComponentDimensions', np.array([1],np.uint64))
                        add_attribute(cell_ensemble[group], 'Tuple Axis Dimensions', f'x={len(self._phases)+1}')
                        add_attribute(cell_ensemble[group], 'DataArrayVersion', np.array([2],np.int32))
                        add_attribute(cell_ensemble[group], 'TupleDimensions', np.array([len(self._phases) + 1],np.uint64))
                    for group in ['CrystalStructures','PhaseTypes']:
                        add_attribute(cell_ensemble[group], 'ObjectType', 'DataArray<uint32_t>')
                    add_attribute(cell_ensemble['PhaseName'], 'ObjectType', 'StringDataArray')

                    geom = create_and_open(data_container,'_SIMPL_GEOMETRY')
                    geom['DIMENSIONS'] = np.array(self.cells,np.int64)
                    geom['ORIGIN']     = np.array(self.origin,np.float32)
                    geom['SPACING']    = np.float32(self.size/self.cells)
                    names = ['GeometryName',  'GeometryTypeName','GeometryType','SpatialDimensionality','UnitDimensionality']
                    values = ['ImageGeometry','ImageGeometry', np.array([0],np.uint32)] + [np.array([3],np.uint32)]*2
                    for name,value in zip(names,values):
                        add_attribute(geom,name,value)


    def export_DADF5(self,
                     fname,
                     output: Union[str, list[str]] = '*',
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
            Indices for regridding. Only applicable for grid
            solver results.
        """
        if Path(fname).expanduser().absolute() == self.fname:
            raise PermissionError(f'cannot overwrite "{self.fname}"')

        if mapping is not None and not self.structured:
            raise PermissionError('cannot regrid unstructured mesh')


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
                if self.version_major == 1 and self.version_minor > 0:
                    f_out['geometry']['cells'][...] = cells

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


            for inc in util.show_progress(self._visible['increments']):
                f_in.copy(inc,f_out,shallow=True)
                if mapping is None:
                    for label in ['u_p','u_n']:
                        f_in[inc]['geometry'].copy(label,f_out[inc]['geometry'])
                else:
                    u_p = f_in[inc]['geometry']['u_p'][()][mapping_flat]
                    f_out[inc]['geometry'].create_dataset('u_p',data=u_p)
                    delta = self.size/np.array(mapping.shape)
                    c_0_p = tuple([np.linspace(delta[i]/2,self.size[i]-delta[i]/2,mapping.shape[i]) for i in [0,1,2]])
                    interpolator = interpolate.RegularGridInterpolator(c_0_p,np.reshape(u_p,tuple(cells)+(3,)),
                                                                       fill_value=None,bounds_error=False,method='linear')
                    c_0_n = grid_filters.coordinates0_node(mapping.shape,self.size).reshape(-1,3)
                    f_out[inc]['geometry'].create_dataset('u_n',data=interpolator(c_0_n))
                    f_out[inc]['geometry/u_n'].attrs.update(f_in[inc]['geometry/u_n'].attrs)


                for label in self._homogenizations:
                    f_in[inc]['homogenization'].copy(label,f_out[inc]['homogenization'],shallow=True)
                for label in self._phases:
                    f_in[inc]['phase'].copy(label,f_out[inc]['phase'],shallow=True)

                for ty in ['phase','homogenization']:
                    for label in self._visible[ty+'s']:
                        for field in _match(self._visible['fields'],f_in['/'.join([inc,ty,label])].keys()):
                            p = '/'.join([inc,ty,label,field])
                            for out in _match(output,f_in[p].keys()):
                                cp(f_in[p],f_out[p],out,None if mapping is None else mappings[ty][label.encode()])


    def export_simulation_setup(self,
                     output: Union[str, list[str]] = '*',
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
                   output: Union[str,list[str]],
                   cfg_dir: Path,
                   overwrite: bool):

            cfg = cfg_dir/name

            if type(obj) is h5py.Dataset and _match(output,[name]):
                if cfg.exists() and not overwrite:
                    raise PermissionError(f'"{cfg}" exists')
                else:
                    cfg.parent.mkdir(parents=True,exist_ok=True)
                    with util.open_text(cfg,'w') as f_out: f_out.write(obj[0].decode())

        cfg_dir = (Path.cwd() if target_dir is None else Path(target_dir))
        with h5py.File(self.fname,'r') as f_in:
            f_in['setup'].visititems(functools.partial(export,
                                                       output=output,
                                                       cfg_dir=cfg_dir,
                                                       overwrite=overwrite))
