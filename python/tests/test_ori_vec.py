import pytest
import numpy as np

from damask import Rotation
from damask import Orientation
from damask import Lattice

rot0= Rotation.from_random()  ; 
rot1= Rotation.from_random()  ; 
rot2= Rotation.from_random()  ; 
rot3= Rotation.from_random()  ; 

#disorientation

#fromaverage
#average

class TestOrientation_vec:
    #@pytest.mark.xfail
    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_equivalent_vec(self,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),rot2.as_quaternion(),rot3.as_quaternion()])
        ori_vec=Orientation(quat,lattice)

        for s in range(len(ori_vec.lattice.symmetry.symmetryOperations())):
            assert all(ori_vec.equivalent_vec.rotation.as_Eulers()[s,0] == \
                        ori0.equivalentOrientations()[s].rotation.as_Eulers())
            assert all(ori_vec.equivalent_vec.rotation.as_quaternion()[s,1] == \
                        ori1.equivalentOrientations()[s].rotation.as_quaternion())
            assert all(ori_vec.equivalent_vec.rotation.as_Rodrigues()[s,2] == \
                        ori2.equivalentOrientations()[s].rotation.as_Rodrigues())
            assert all(ori_vec.equivalent_vec.rotation.as_cubochoric()[s,3] == \
                        ori3.equivalentOrientations()[s].rotation.as_cubochoric())

    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_inFZ_vec(self,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)
        ori4=ori0.reduced() ; rot4=ori4.rotation #ensure 1 of them is in FZ

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),\
                        rot2.as_quaternion(),rot3.as_quaternion(), rot4.as_quaternion()])
        ori_vec=Orientation(quat,lattice)

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
        ori_vec=Orientation(quat,lattice)


        for s in range(len(ori1.lattice.relationOperations(model)['rotations'])):
            assert all(ori_vec.relatedOrientations_vec(model).rotation.as_Eulers()[s,0] == \
                        ori0.relatedOrientations(model)[s].rotation.as_Eulers())
            assert all(ori_vec.relatedOrientations_vec(model).rotation.as_quaternion()[s,1] == \
                        ori1.relatedOrientations(model)[s].rotation.as_quaternion())
            assert all(ori_vec.relatedOrientations_vec(model).rotation.as_Rodrigues()[s,2] == \
                        ori2.relatedOrientations(model)[s].rotation.as_Rodrigues())
            assert all(ori_vec.relatedOrientations_vec(model).rotation.as_cubochoric()[s,3] == \
                        ori3.relatedOrientations(model)[s].rotation.as_cubochoric())

    @pytest.mark.parametrize('lattice',Lattice.lattices)
    def test_reduced_vec(self,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)
        #ensure 1 of them is in FZ
        ori4=ori0.reduced()
        rot4=ori4.rotation

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),\
                        rot2.as_quaternion(),rot3.as_quaternion(), rot4.as_quaternion()])
        ori_vec=Orientation(quat,lattice)

        assert all(ori_vec.reduced_vec.rotation.as_Eulers()[0] ==     ori0.reduced().rotation.as_Eulers()     )
        assert all(ori_vec.reduced_vec.rotation.as_quaternion()[1] == ori1.reduced().rotation.as_quaternion() )
        assert all(ori_vec.reduced_vec.rotation.as_Rodrigues()[2] ==  ori2.reduced().rotation.as_Rodrigues()  )
        assert all(ori_vec.reduced_vec.rotation.as_cubochoric()[3] == ori3.reduced().rotation.as_cubochoric() )
        assert all(ori_vec.reduced_vec.rotation.as_axis_angle()[4] == ori4.reduced().rotation.as_axis_angle() )


    @pytest.mark.parametrize('lattice',['bcc','fcc','bct'])
    def test_IPFcolor_vec(self,lattice):
        ori0=Orientation(rot0,lattice)
        ori1=Orientation(rot1,lattice)
        ori2=Orientation(rot2,lattice)
        ori3=Orientation(rot3,lattice)

        quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),\
                        rot2.as_quaternion(),rot3.as_quaternion()])
        ori_vec=Orientation(quat,lattice)

        assert np.allclose( ori_vec.IPFcolor_vec(np.array([0,0,1]))[0],ori0.IPFcolor(np.array([0,0,1])))
        assert np.allclose( ori_vec.IPFcolor_vec(np.array([0,2,1]))[1],ori1.IPFcolor(np.array([0,2,1])))
        assert np.allclose( ori_vec.IPFcolor_vec(np.array([0,3,1]))[2],ori2.IPFcolor(np.array([0,3,1])))
        assert np.allclose( ori_vec.IPFcolor_vec(np.array([4,0,1]))[3],ori3.IPFcolor(np.array([4,0,1])))
