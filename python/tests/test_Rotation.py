import os

import pytest
import numpy as np

from damask import Rotation

n = 1000

@pytest.fixture
def default():
    """A set of n random rotations."""
    specials =[np.array([ 1.0, 0.0, 0.0, 0.0]),
               #-----------------------------------------------
               np.array([ 1.0, 1.0, 0.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 1.0, 0.0, 1.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 1.0, 0.0, 0.0, 1.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 1.0, 1.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 1.0, 0.0, 1.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 0.0, 1.0, 1.0])*np.sqrt(2.)*.5,
               #-----------------------------------------------
               np.array([ 1.0,-1.0, 0.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 1.0, 0.0,-1.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 1.0, 0.0, 0.0,-1.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 1.0,-1.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 1.0, 0.0,-1.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 0.0, 1.0,-1.0])*np.sqrt(2.)*.5,
               #-----------------------------------------------
               np.array([ 0.0, 1.0,-1.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 1.0, 0.0,-1.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 0.0, 1.0,-1.0])*np.sqrt(2.)*.5,
               #-----------------------------------------------
               np.array([ 0.0,-1.0,-1.0, 0.0])*np.sqrt(2.)*.5,
               np.array([ 0.0,-1.0, 0.0,-1.0])*np.sqrt(2.)*.5,
               np.array([ 0.0, 0.0,-1.0,-1.0])*np.sqrt(2.)*.5,
               #-----------------------------------------------
              ]
    return [Rotation.fromQuaternion(s) for s in specials] + \
           [Rotation.fromRandom() for r in range(n-len(specials))]

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Rotation')


class TestRotation:

    def test_Eulers(self,default):
        for rot in default:
            c = Rotation.fromEulers(rot.asEulers())
            ok = np.allclose(rot.asQuaternion(),c.asQuaternion())
            if np.isclose(rot.asQuaternion()[0],0.0,atol=1.e-13,rtol=0.0):
                ok = ok or np.allclose(rot.asQuaternion(),c.asQuaternion()*-1.)
            assert ok

    def test_AxisAngle(self,default):
        for rot in default:
            c = Rotation.fromAxisAngle(rot.asAxisAngle())
            assert np.allclose(rot.asEulers(),c.asEulers())

    def test_Matrix(self,default):
        for rot in default:
            c = Rotation.fromMatrix(rot.asMatrix())
            ok = np.allclose(rot.asAxisAngle(),c.asAxisAngle())
            if np.isclose(rot.asAxisAngle()[3],np.pi):
                ok = ok or np.allclose(rot.asAxisAngle(),c.asAxisAngle()*np.array([-1.,-1.,-1.,1.]))
            assert ok

    def test_Rodriques(self,default):
        for rot in default:
            c = Rotation.fromRodrigues(rot.asRodrigues())
            assert np.allclose(rot.asMatrix(),c.asMatrix())

    def test_Homochoric(self,default):
        for rot in default:
            c = Rotation.fromHomochoric(rot.asHomochoric())
            assert np.allclose(np.clip(rot.asRodrigues(),None,1.e9),np.clip(c.asRodrigues(),None,1.e9))

    def test_Cubochoric(self,default):
        for rot in default:
            c = Rotation.fromCubochoric(rot.asCubochoric())
            assert np.allclose(rot.asHomochoric(),c.asHomochoric())

    def test_Quaternion(self,default):
        for rot in default:
            c = Rotation.fromQuaternion(rot.asQuaternion())
            assert np.allclose(rot.asCubochoric(),c.asCubochoric())
