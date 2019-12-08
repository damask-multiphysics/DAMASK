import os

import pytest
import numpy as np

import damask
from damask import Rotation
from damask import Orientation
   
n = 1000

@pytest.fixture
def default():
    """A set of n random rotations."""
    return [Rotation.fromRandom() for r in range(n)]

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Rotation')


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


    @pytest.mark.parametrize('model',['Bain','KS','GT','GT_prime','NW','Pitsch'])
    @pytest.mark.parametrize('lattice',['fcc','bcc'])
    def test_relationship_forward_backward(self,model,lattice):
        ori = Orientation(Rotation.fromRandom(),lattice)
        for i,r in enumerate(ori.relatedOrientations(model)):
            ori2 = r.relatedOrientations(model)[i]
            misorientation = ori.rotation.misorientation(ori2.rotation)
            assert misorientation.asAxisAngle(degrees=True)[3]<1.0e-5

    @pytest.mark.parametrize('model',['Bain','KS','GT','GT_prime','NW','Pitsch'])
    @pytest.mark.parametrize('lattice',['fcc','bcc'])
    def test_relationship_reference(self,update,reference_dir,model,lattice):
        reference = os.path.join(reference_dir,'{}_{}.txt'.format(lattice,model))
        ori = Orientation(Rotation(),lattice)
        eu = np.array([o.rotation.asEulers(degrees=True) for o in ori.relatedOrientations(model)])
        if update: 
            coords = np.array([(1,i+1) for i,x in enumerate(eu)])
            table = damask.Table(eu,{'Eulers':(3,)})
            table.add('pos',coords)
            table.to_ASCII(reference)
        assert np.allclose(eu,damask.Table.from_ASCII(reference).get('Eulers'))
