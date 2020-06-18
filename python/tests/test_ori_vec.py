import pytest
import numpy as np

from damask import Rotation
from damask import Orientation
from damask import Lattice

rot0= Rotation.from_random()
rot1= Rotation.from_random()
rot2= Rotation.from_random()
rot3= Rotation.from_random()

class TestOrientation_vec:
    @pytest.mark.xfail
    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_equivalentOrientations_vec(self,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),rot2.as_quaternion(),rot3.as_quaternion()])
        rot_vec=Rotation.from_quaternion(quat)
        ori_vec=Orientation(rot_vec,lattice)

        for s in range(len(ori_vec.lattice.symmetry.symmetryOperations())):
            assert all(ori_vec.equivalent_vec()[s,0].rotation.as_Eulers() == \
                        ori0.equivalentOrientations()[s].rotation.as_Eulers())
            assert all(ori_vec.equivalent_vec()[s,1].rotation.as_quaternion() == \
                        ori1.equivalentOrientations()[s].rotation.as_quaternion())
            assert all(ori_vec.equivalent_vec()[s,2].rotation.as_Rodrigues() == \
                        ori2.equivalentOrientations()[s].rotation.as_Rodrigues())
            assert all(ori_vec.equivalent_vec()[s,3].rotation.as_cubochoric() == \
                        ori3.equivalentOrientations()[s].rotation.as_cubochoric())

    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_inFZ_vec(self,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)
        #ensure 1 of them is in FZ
        ori4=ori0.reduced()
        rot4=ori4.rotation

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),\
                        rot2.as_quaternion(),rot3.as_quaternion(), rot4.as_quaternion()])
        rot_vec=Rotation.from_quaternion(quat)
        ori_vec=Orientation(rot_vec,lattice)

        assert ori_vec.inFZ_vec()[0] == ori0.inFZ()
        assert ori_vec.inFZ_vec()[1] == ori1.inFZ()
        assert ori_vec.inFZ_vec()[2] == ori2.inFZ()
        assert ori_vec.inFZ_vec()[3] == ori3.inFZ()
        assert ori_vec.inFZ_vec()[4] == ori4.inFZ()


    @pytest.mark.parametrize('model',['Bain','KS','GT','GT_prime','NW','Pitsch'])
    @pytest.mark.parametrize('lattice',['fcc','bcc'])
    def test_relatedOrientations_vec(self,model,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),rot2.as_quaternion(),rot3.as_quaternion()])
        rot_vec=Rotation.from_quaternion(quat)
        ori_vec=Orientation(rot_vec,lattice)

        for s in range(len(ori1.lattice.relationOperations(model)['rotations'])):
            assert all(ori_vec.relatedOrientations_vec(model)[s,0].rotation.as_Eulers() == \
                        ori0.relatedOrientations(model)[s].rotation.as_Eulers())
            assert all(ori_vec.relatedOrientations_vec(model)[s,1].rotation.as_quaternion() == \
                        ori1.relatedOrientations(model)[s].rotation.as_quaternion())
            assert all(ori_vec.relatedOrientations_vec(model)[s,2].rotation.as_Rodrigues() == \
                        ori2.relatedOrientations(model)[s].rotation.as_Rodrigues())
            assert all(ori_vec.relatedOrientations_vec(model)[s,3].rotation.as_cubochoric() == \
                        ori3.relatedOrientations(model)[s].rotation.as_cubochoric())

