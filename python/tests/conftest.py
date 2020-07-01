from pathlib import Path
import numpy as np

import pytest

def pytest_addoption(parser):
    parser.addoption("--update",
                     action="store_true",
                     default=False)

@pytest.fixture
def update(request):
    """Store current results as new reference results."""
    return request.config.getoption("--update")

@pytest.fixture
def reference_dir_base():
    """Directory containing reference results."""
    return Path(__file__).parent/'reference'

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

    n = 1100
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
    specials_scatter = specials + np.broadcast_to(np.random.rand(4)*scatter,specials.shape)
    specials_scatter /= np.linalg.norm(specials_scatter,axis=1).reshape(-1,1)
    specials_scatter[specials_scatter[:,0]<0]*=-1

    return np.array([s for s in specials] + \
                    [s for s in specials_scatter] + \
                    [s for s in random_quaternions(n-2*len(specials))])
