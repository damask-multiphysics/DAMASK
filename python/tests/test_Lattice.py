import random

import pytest
import numpy as np

from damask import Rotation
from damask import Symmetry

def in_FZ(system,rho):
    """Non-vectorized version of 'in_FZ'."""
    rho_abs = abs(rho)

    if system == 'cubic':
        return     np.sqrt(2.0)-1.0 >= rho_abs[0] \
               and np.sqrt(2.0)-1.0 >= rho_abs[1] \
               and np.sqrt(2.0)-1.0 >= rho_abs[2] \
               and 1.0 >= rho_abs[0] + rho_abs[1] + rho_abs[2]
    elif system == 'hexagonal':
        return     1.0 >= rho_abs[0] and 1.0 >= rho_abs[1] and 1.0 >= rho_abs[2] \
               and 2.0 >= np.sqrt(3)*rho_abs[0] + rho_abs[1] \
               and 2.0 >= np.sqrt(3)*rho_abs[1] + rho_abs[0] \
               and 2.0 >= np.sqrt(3) + rho_abs[2]
    elif system == 'tetragonal':
        return     1.0 >= rho_abs[0] and 1.0 >= rho_abs[1] \
               and np.sqrt(2.0) >= rho_abs[0] + rho_abs[1] \
               and np.sqrt(2.0) >= rho_abs[2] + 1.0
    elif system == 'orthorhombic':
        return     1.0 >= rho_abs[0] and 1.0 >= rho_abs[1] and 1.0 >= rho_abs[2]
    else:
        return np.all(np.isfinite(rho_abs))


def in_disorientation_SST(system,rho):
    """Non-vectorized version of 'in_Disorientation_SST'."""
    epsilon = 0.0
    if system == 'cubic':
        return rho[0] >= rho[1]+epsilon              and rho[1] >= rho[2]+epsilon and rho[2] >= epsilon
    elif system == 'hexagonal':
        return rho[0] >= np.sqrt(3)*(rho[1]-epsilon) and rho[1] >= epsilon        and rho[2] >= epsilon
    elif system == 'tetragonal':
        return rho[0] >= rho[1]-epsilon              and rho[1] >= epsilon        and rho[2] >= epsilon
    elif system == 'orthorhombic':
        return rho[0] >= epsilon                     and rho[1] >= epsilon        and rho[2] >= epsilon
    else:
        return True


def in_SST(system,vector,proper = False):
    """Non-vectorized version of 'in_SST'."""
    if system == 'cubic':
        basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                       [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                       [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                   'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                       [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                       [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                }
    elif system == 'hexagonal':
        basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                       [ 1.            , -np.sqrt(3.)   ,  0. ],
                                       [ 0.            ,  2.            ,  0. ] ]),
                 'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                       [-1.            ,  np.sqrt(3.)   ,  0. ],
                                       [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                }
    elif system == 'tetragonal':
        basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                       [ 1.            , -1.            ,  0. ],
                                       [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                 'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                       [-1.            ,  1.            ,  0. ],
                                       [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                }
    elif system == 'orthorhombic':
        basis = {'improper':np.array([ [ 0., 0., 1.],
                                       [ 1., 0., 0.],
                                       [ 0., 1., 0.] ]),
                   'proper':np.array([ [ 0., 0., 1.],
                                       [-1., 0., 0.],
                                       [ 0., 1., 0.] ]),
                }
    else:
        return True

    v = np.array(vector,dtype=float)
    if proper:
        theComponents = np.around(np.dot(basis['improper'],v),12)
        inSST = np.all(theComponents >= 0.0)
        if not inSST:
            theComponents = np.around(np.dot(basis['proper'],v),12)
            inSST = np.all(theComponents >= 0.0)
    else:
        v[2] = abs(v[2])
        theComponents = np.around(np.dot(basis['improper'],v),12)
        inSST = np.all(theComponents >= 0.0)

    return inSST


@pytest.fixture
def set_of_rodrigues(set_of_quaternions):
    return Rotation(set_of_quaternions).as_Rodrigues(vector=True)[:200]

class TestSymmetry:

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_in_FZ_vectorize(self,set_of_rodrigues,system):
        result = Symmetry(system).in_FZ(set_of_rodrigues.reshape(50,4,3)).reshape(200)
        for i,r in enumerate(result):
            assert r == in_FZ(system,set_of_rodrigues[i])

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_in_disorientation_SST_vectorize(self,set_of_rodrigues,system):
        result = Symmetry(system).in_disorientation_SST(set_of_rodrigues.reshape(50,4,3)).reshape(200)
        for i,r in enumerate(result):
            assert r == in_disorientation_SST(system,set_of_rodrigues[i])

    @pytest.mark.parametrize('proper',[True,False])
    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_in_SST_vectorize(self,system,proper):
        vecs = np.random.rand(20,4,3)
        result = Symmetry(system).in_SST(vecs,proper).reshape(20*4)
        for i,r in enumerate(result):
            assert r == in_SST(system,vecs.reshape(20*4,3)[i],proper)

    @pytest.mark.parametrize('invalid_symmetry',['fcc','bcc','hello'])
    def test_invalid_symmetry(self,invalid_symmetry):
        with pytest.raises(KeyError):
            s = Symmetry(invalid_symmetry)                                                          # noqa

    def test_equal(self):
        symmetry = random.choice(Symmetry.crystal_systems)
        print(symmetry)
        assert Symmetry(symmetry) == Symmetry(symmetry)

    def test_not_equal(self):
        symmetries = random.sample(Symmetry.crystal_systems,k=2)
        assert Symmetry(symmetries[0]) != Symmetry(symmetries[1])

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_in_FZ(self,system):
          assert Symmetry(system).in_FZ(np.zeros(3))

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_in_disorientation_SST(self,system):
          assert Symmetry(system).in_disorientation_SST(np.zeros(3))

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    @pytest.mark.parametrize('proper',[True,False])
    def test_in_SST(self,system,proper):
          assert Symmetry(system).in_SST(np.zeros(3),proper)

    @pytest.mark.parametrize('function',['in_FZ','in_disorientation_SST','in_SST'])
    def test_invalid_argument(self,function):
        s = Symmetry()                                                                              # noqa
        with pytest.raises(ValueError):
            eval(f's.{function}(np.ones(4))')
