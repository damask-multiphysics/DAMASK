# SPDX-License-Identifier: AGPL-3.0-or-later
from pathlib import Path
import os
import glob
import shutil
import re

import h5py
import numpy as np
import pytest

import damask

def pytest_addoption(parser):
    parser.addoption('--update', action='store_true', default=False,
                     help='Update reference results.')
    parser.addoption('--damask-root',default=os.environ.get('DAMASK_ROOT'),
                     help='DAMASK root directory.')
    parser.addoption('--rng-entropy',
                     help='Entropy for random seed generator.')
    parser.addoption('--mpi-launcher',default='mpiexec',
                     help='MPI launcher executable.')

@pytest.fixture
def update(request):
    """Store current results as new reference data."""
    return request.config.getoption('--update')

@pytest.fixture
def damask_root(request):
    """DAMASK root directory."""
    if (damask_root := request.config.getoption('--damask-root')) is not None:
        return Path(damask_root)
    else:
        return damask_root

@pytest.fixture
def np_rng(request):
    """Instance of numpy.random.Generator."""
    e = request.config.getoption('--rng-entropy')
    print('\nrng entropy:',sq := np.random.SeedSequence(e if e is None else int(e)).entropy)
    return np.random.default_rng(seed=sq)

@pytest.fixture
def mpi_launcher(request):
    """Name of the MPI launcher."""
    return request.config.getoption('--mpi-launcher')


@pytest.fixture
def res_path_base():
    """Directory containing testing resources."""
    return Path(__file__).parent/'resources'


@pytest.fixture
def copy_files():
    """Copy files from source to destination directory."""
    def _copy_files(src,dst,files=None):
        for f in (files if files else glob.glob(os.path.join(src,'*'))):
            shutil.copy(os.path.join(src,f),dst)
    return _copy_files


@pytest.fixture
def h5py_dataset_iterator():
    """Iterate over all datasets in an HDF5 file."""
    def _h5py_dataset_iterator(g, prefix=''):
        for key,item in g.items():
            path = os.path.join(prefix, key)
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from _h5py_dataset_iterator(item, path)
    return _h5py_dataset_iterator


@pytest.fixture
def assert_allclose():
    """
    Asserts the element-wise equality of two arrays within relative+absolute tolerance.

    Parameters
    ----------
    a : np.ndarray or h5py Dataset
        Array to compare.
    b : np.ndarray or h5py Dataset
        Array to compare.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.
    N : int
        Maximum number of reported deviating values.
        Defaults to 16. N=None outputs all values.
    msg : str
        Informative message.

    """
    def _assert_allclose(a,b,
                         rtol: float = 1e-05,
                         atol: float = 1e-08,
                         N = 16,
                         msg = '',
                         ):
        assert np.logical_or(np.isclose(a,b, rtol=rtol,atol=atol),
                             np.isclose(b,a, rtol=rtol,atol=atol)).all(), \
               report_nonclose(a,b, rtol=rtol,atol=atol, N=N, msg=msg)

    return _assert_allclose

def report_nonclose(a, b, rtol, atol, N, msg):
    """
    Report where values of two arrays deviate from each other.

    Output is sorted from large to small magnitude of absolute difference.

    Parameters
    ----------
    a : np.ndarray or h5py Dataset
        Array to compare.
    b : np.ndarray or h5py Dataset
        Array to compare.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.
    N : int
        Maximum number of reported deviating values.
        Defaults to all.
    msg : str
        Informative message added to output.
        Defaults to None.

    Returns
    -------
    description : str
        Description of differences.

    """
    a_ = np.atleast_1d(np.asarray(a))
    b_ = np.atleast_1d(np.asarray(b))
    if a_.shape != b_.shape: b_ = np.broadcast_to(b_,a.shape)

    absdiff = np.abs(a_ - b_)
    absmax  = np.max(np.abs(np.stack((a_, b_))), axis=0)

    idx = (absdiff > atol + rtol * absmax).nonzero()
    ids = np.argsort(absdiff[idx])[::-1]
    n = len(idx[0])
    split = N is not None and N<n
    head,tail = (slice(0,(N+1)//2),slice(n-N//2,n)) if split else (slice(0,n),slice(n,n))
    diffs = [ f'abs / rel diff {absdiff[idx][i]:>16.10g} / {absdiff[idx][i]/absmax[idx][i]:<16.10g}'
             +f' between {a_[idx][i]:>16.10g} and {b_[idx][i]:<16.10g}'
             +f' at {np.transpose(idx)[i]}' for i in list(ids[head])+list(ids[tail])]
    if split: diffs.insert(head.stop,'...')

    return '\n'.join(['']+diffs
                     +['',f'fraction {n}/{a_.size} outside tolerance']
                     +(['',msg] if msg is not None else [])
                     )


def petsc_version_():
    """Return PETSc version."""
    if (petsc_dir := os.environ.get('PETSC_DIR')):
        with open(Path(petsc_dir)/'include/petscversion.h','r') as f:
            c = f.read()
            v = '.'.join([re.search(fr'#define PETSC_VERSION_{p}\s*([0-9]+)+',c).group(1) \
                          for p in ['MAJOR','MINOR','SUBMINOR']])
            print(f'PETSc version: {v}')
        return damask.util.version(v)
    else:
        raise EnvironmentError('environment variable "PETSC_DIR" not defined')


pytest.petsc_version = petsc_version_


@pytest.fixture
def petsc_version():
    return petsc_version_
