import pytest
import numpy as np

from damask import Rotation
   
n = 1000

@pytest.fixture
def default():
    """A set of n random rotations."""
    return [Rotation.fromRandom() for r in range(n)]


class TestRotation:

    def test_Eulers(self,default):
        for rot in default:
            assert np.allclose(rot.asQuaternion(), 
                               Rotation.fromEulers(rot.asEulers()).asQuaternion())


    def test_AxisAngle(self,default):
        for rot in default:
            assert np.allclose(rot.asEulers(),
                               Rotation.fromAxisAngle(rot.asAxisAngle()).asEulers())


    def test_Matrix(self,default):
        for rot in default:
            assert np.allclose(rot.asAxisAngle(),
                               Rotation.fromMatrix(rot.asMatrix()).asAxisAngle())


    def test_Rodriques(self,default):
        for rot in default:
            assert np.allclose(rot.asMatrix(),
                               Rotation.fromRodrigues(rot.asRodrigues()).asMatrix())


    def test_Homochoric(self,default):
        for rot in default:
            assert np.allclose(rot.asRodrigues(),
                               Rotation.fromHomochoric(rot.asHomochoric()).asRodrigues())


    def test_Cubochoric(self,default):
        for rot in default:
            assert np.allclose(rot.asHomochoric(),
                               Rotation.fromCubochoric(rot.asCubochoric()).asHomochoric())


    def test_Quaternion(self,default):
        for rot in default:
            assert np.allclose(rot.asCubochoric(),
                               Rotation.fromQuaternion(rot.asQuaternion()).asCubochoric())
