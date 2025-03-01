from pathlib import Path
import os
import datetime

import numpy as np
import pytest
import matplotlib as mpl
if os.name == 'posix' and 'DISPLAY' not in os.environ:
    mpl.use('Agg')
import matplotlib.pyplot as plt

import damask


patched_version = '99.99.99-9999-pytest'
@pytest.fixture
def patch_damask_version(monkeypatch):
    """Set damask.version for reproducible tests results."""
    monkeypatch.setattr(damask, 'version', patched_version)


patched_date = datetime.datetime(2019, 11, 2, 11, 58, 0)
@pytest.fixture
def patch_datetime_now(monkeypatch):
    """Set datetime.datetime.now for reproducible tests results."""
    class mydatetime:
        @classmethod
        def now(cls):
            return patched_date

    monkeypatch.setattr(datetime, 'datetime', mydatetime)


@pytest.fixture
def patch_execution_stamp(monkeypatch):
    """Set damask.util.execution_stamp for reproducible tests results."""
    def execution_stamp(class_name,function_name=None):
        _function_name = '' if function_name is None else f'.{function_name}'
        return f'damask.{class_name}{_function_name} v{patched_version} ({patched_date})'

    monkeypatch.setattr(damask.util, 'execution_stamp', execution_stamp)


@pytest.fixture
def patch_plt_show(monkeypatch):
    def _None(block=None):
        pass
    monkeypatch.setattr(plt, 'show', _None, raising=True)


def pytest_addoption(parser):
    parser.addoption('--update', action='store_true', default=False,
                     help='Update reference results.')
    parser.addoption('--damask-root',default=os.environ.get('DAMASK_ROOT'),
                     help='DAMASK root directory.')
    parser.addoption('--rng-entropy',
                     help='Entropy for random seed generator.')

@pytest.fixture
def update(request):
    """Store current results as new reference data."""
    return request.config.getoption('--update')

@pytest.fixture
def damask_root(request):
    """DAMASK root directory."""
    if (damask_root := request.config.getoption('--damask-root')) is not None:
        return Path(damask_root).expanduser()
    else:
        return damask_root

@pytest.fixture
def np_rng(request):
    """Instance of numpy.random.Generator."""
    e = request.config.getoption('--rng-entropy')
    print('\nrng entropy: ',sq := np.random.SeedSequence(e if e is None else int(e)).entropy)
    return np.random.default_rng(seed=sq)

# https://stackoverflow.com/questions/51883573
def pytest_collection_modifyitems(config, items):
    if config.getoption('--damask-root') is None:
        need_damask_root = pytest.mark.skip(reason='need --damask-root to run')
        for item in items:
            if 'need_damask_root' in item.keywords: item.add_marker(need_damask_root)

def pytest_configure(config):
    config.addinivalue_line(
        'markers', 'need_damask_root: mark test to run only if DAMASK root is given'
    )


@pytest.fixture
def res_path_base():
    """Directory containing testing resources."""
    return Path(__file__).parent/'resources'


@pytest.fixture
def set_of_quaternions(np_rng):
    """A set of n random rotations."""
    def random_quaternions(N,rng):
        r = np_rng.random((N,3))

        A = np.sqrt(r[:,2])
        B = np.sqrt(1.0-r[:,2])
        qu = np.column_stack([np.cos(2.0*np.pi*r[:,0])*A,
                              np.sin(2.0*np.pi*r[:,1])*B,
                              np.cos(2.0*np.pi*r[:,1])*B,
                              np.sin(2.0*np.pi*r[:,0])*A])
        qu[:,0]*=np.sign(qu[:,0])

        return qu

    n = 600
    scatter=1.e-2
    specials = np.array([
                         [1.0, 0.0, 0.0, 0.0],
                        #----------------------
                         [0.0, 1.0, 0.0, 0.0],
                         [0.0, 0.0, 1.0, 0.0],
                         [0.0, 0.0, 0.0, 1.0],
                         [0.0,-1.0, 0.0, 0.0],
                         [0.0, 0.0,-1.0, 0.0],
                         [0.0, 0.0, 0.0,-1.0],
                        #----------------------
                         [1.0, 1.0, 0.0, 0.0],
                         [1.0, 0.0, 1.0, 0.0],
                         [1.0, 0.0, 0.0, 1.0],
                         [0.0, 1.0, 1.0, 0.0],
                         [0.0, 1.0, 0.0, 1.0],
                         [0.0, 0.0, 1.0, 1.0],
                        #----------------------
                         [1.0,-1.0, 0.0, 0.0],
                         [1.0, 0.0,-1.0, 0.0],
                         [1.0, 0.0, 0.0,-1.0],
                         [0.0, 1.0,-1.0, 0.0],
                         [0.0, 1.0, 0.0,-1.0],
                         [0.0, 0.0, 1.0,-1.0],
                        #----------------------
                         [0.0, 1.0,-1.0, 0.0],
                         [0.0, 1.0, 0.0,-1.0],
                         [0.0, 0.0, 1.0,-1.0],
                        #----------------------
                         [0.0,-1.0,-1.0, 0.0],
                         [0.0,-1.0, 0.0,-1.0],
                         [0.0, 0.0,-1.0,-1.0],
                        #----------------------
                         [1.0, 1.0, 1.0, 0.0],
                         [1.0, 1.0, 0.0, 1.0],
                         [1.0, 0.0, 1.0, 1.0],
                         [1.0,-1.0, 1.0, 0.0],
                         [1.0,-1.0, 0.0, 1.0],
                         [1.0, 0.0,-1.0, 1.0],
                         [1.0, 1.0,-1.0, 0.0],
                         [1.0, 1.0, 0.0,-1.0],
                         [1.0, 0.0, 1.0,-1.0],
                         [1.0,-1.0,-1.0, 0.0],
                         [1.0,-1.0, 0.0,-1.0],
                         [1.0, 0.0,-1.0,-1.0],
                        #----------------------
                         [0.0, 1.0, 1.0, 1.0],
                         [0.0, 1.0,-1.0, 1.0],
                         [0.0, 1.0, 1.0,-1.0],
                         [0.0,-1.0, 1.0, 1.0],
                         [0.0,-1.0,-1.0, 1.0],
                         [0.0,-1.0, 1.0,-1.0],
                         [0.0,-1.0,-1.0,-1.0],
                        #----------------------
                         [1.0, 1.0, 1.0, 1.0],
                         [1.0,-1.0, 1.0, 1.0],
                         [1.0, 1.0,-1.0, 1.0],
                         [1.0, 1.0, 1.0,-1.0],
                         [1.0,-1.0,-1.0, 1.0],
                         [1.0,-1.0, 1.0,-1.0],
                         [1.0, 1.0,-1.0,-1.0],
                         [1.0,-1.0,-1.0,-1.0],
                        ])
    specials /= np.linalg.norm(specials,axis=1).reshape(-1,1)
    specials_scatter = specials + np.broadcast_to((np_rng.random(4)*2.-1.)*scatter,specials.shape)
    specials_scatter /= np.linalg.norm(specials_scatter,axis=1).reshape(-1,1)
    specials_scatter[specials_scatter[:,0]<0]*=-1

    return np.vstack((specials,
                      specials_scatter,
                      random_quaternions(n-2*len(specials),np_rng),
                      ))

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
    a_ = np.asarray(a)
    b_ = np.asarray(b)
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
