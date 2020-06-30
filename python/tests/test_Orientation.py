import os
from itertools import permutations

import pytest
import numpy as np

import damask
from damask import Rotation
from damask import Orientation
from damask import Lattice

n = 1000

def IPF_color(orientation,direction):
    """TSL color of inverse pole figure for given axis (non-vectorized)."""
    for o in orientation.equivalent:
        pole = o.rotation@direction
        inSST,color = orientation.lattice.in_SST(pole,color=True)
        if inSST: break

    return color

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Rotation')


class TestOrientation:

    @pytest.mark.parametrize('color',[{'label':'red',  'RGB':[1,0,0],'direction':[0,0,1]},
                                      {'label':'green','RGB':[0,1,0],'direction':[0,1,1]},
                                      {'label':'blue', 'RGB':[0,0,1],'direction':[1,1,1]}])
    @pytest.mark.parametrize('lattice',['fcc','bcc'])
    def test_IPF_cubic(self,color,lattice):
        cube = damask.Orientation(damask.Rotation(),lattice)
        for direction in set(permutations(np.array(color['direction']))):
            assert np.allclose(cube.IPF_color(np.array(direction)),np.array(color['RGB']))

    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_IPF_equivalent(self,set_of_quaternions,lattice):
        direction = np.random.random(3)*2.0-1
        for ori in Orientation(Rotation(set_of_quaternions),lattice)[200]:
            color = ori.IPF_color(direction)
            for equivalent in ori.equivalent:
                assert np.allclose(color,equivalent.IPF_color(direction))


    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_IPF_vectorize(self,set_of_quaternions,lattice):
        for ori in Orientation(Rotation(set_of_quaternions),lattice)[200]:
            direction = np.random.random(3)*2.0-1
            assert np.allclose(ori.IPF_color(direction),IPF_color(ori,direction))


    @pytest.mark.parametrize('model',['Bain','KS','GT','GT_prime','NW','Pitsch'])
    @pytest.mark.parametrize('lattice',['fcc','bcc'])
    def test_relationship_forward_backward(self,model,lattice):
        ori = Orientation(Rotation.from_random(),lattice)
        for i,r in enumerate(ori.relatedOrientations(model)):
            ori2 = r.relatedOrientations(model)[i]
            misorientation = ori.rotation.misorientation(ori2.rotation)
            assert misorientation.asAxisAngle(degrees=True)[3]<1.0e-5

    @pytest.mark.parametrize('model',['Bain','KS','GT','GT_prime','NW','Pitsch'])
    @pytest.mark.parametrize('lattice',['fcc','bcc'])
    def test_relationship_reference(self,update,reference_dir,model,lattice):
        reference = os.path.join(reference_dir,'{}_{}.txt'.format(lattice,model))
        ori = Orientation(Rotation(),lattice)
        eu = np.array([o.rotation.as_Eulers(degrees=True) for o in ori.relatedOrientations(model)])
        if update:
            coords = np.array([(1,i+1) for i,x in enumerate(eu)])
            table = damask.Table(eu,{'Eulers':(3,)})
            table.add('pos',coords)
            table.to_ASCII(reference)
        assert np.allclose(eu,damask.Table.from_ASCII(reference).get('Eulers'))
