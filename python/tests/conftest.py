from pathlib import Path
import datetime
import os

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
    parser.addoption('--damaskroot',
                     help='DAMASK root directory.')

@pytest.fixture
def update(pytestconfig):
    """Store current results as new reference results."""
    return pytestconfig.getoption('--update')

@pytest.fixture
def damaskroot(pytestconfig):
    """Specify DAMASK root directory."""
    return pytestconfig.getoption('--damaskroot')

# https://stackoverflow.com/questions/51883573
def pytest_collection_modifyitems(config, items):
    if config.getoption('--damaskroot') is None:
        need_damaskroot = pytest.mark.skip(reason='need --damaskroot to run')
        for item in items:
            if 'need_damaskroot' in item.keywords: item.add_marker(need_damaskroot)

def pytest_configure(config):
    config.addinivalue_line(
        'markers', 'need_damaskroot: mark test to run only if DAMASK root is given'
    )

@pytest.fixture
def res_path_base():
    """Directory containing testing resources."""
    return Path(__file__).parent/'resources'


@pytest.fixture
def set_of_quaternions():
    """A set of n random rotations."""
    def random_quaternions(N):
        r = np.random.rand(N,3)

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
    specials_scatter = specials + np.broadcast_to((np.random.rand(4)*2.-1.)*scatter,specials.shape)
    specials_scatter /= np.linalg.norm(specials_scatter,axis=1).reshape(-1,1)
    specials_scatter[specials_scatter[:,0]<0]*=-1

    return np.vstack((specials,
                      specials_scatter,
                      random_quaternions(n-2*len(specials)),
                      ))
