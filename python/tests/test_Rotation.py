import os

import pytest
import numpy as np

from damask import Rotation
import rotation_conversion

n = 1100
atol=1.e-4
scatter=1.e-2

@pytest.fixture
def default():
    """A set of n random rotations."""
    specials = np.array(
              [np.array([ 1.0, 0.0, 0.0, 0.0]),
               #-----------------------------------------------
               np.array([0.0, 1.0, 0.0, 0.0]),
               np.array([0.0, 0.0, 1.0, 0.0]),
               np.array([0.0, 0.0, 0.0, 1.0]),
               np.array([0.0,-1.0, 0.0, 0.0]),
               np.array([0.0, 0.0,-1.0, 0.0]),
               np.array([0.0, 0.0, 0.0,-1.0]),
               #-----------------------------------------------
               np.array([1.0, 1.0, 0.0, 0.0])/np.sqrt(2.),
               np.array([1.0, 0.0, 1.0, 0.0])/np.sqrt(2.),
               np.array([1.0, 0.0, 0.0, 1.0])/np.sqrt(2.),
               np.array([0.0, 1.0, 1.0, 0.0])/np.sqrt(2.),
               np.array([0.0, 1.0, 0.0, 1.0])/np.sqrt(2.),
               np.array([0.0, 0.0, 1.0, 1.0])/np.sqrt(2.),
               #-----------------------------------------------
               np.array([1.0,-1.0, 0.0, 0.0])/np.sqrt(2.),
               np.array([1.0, 0.0,-1.0, 0.0])/np.sqrt(2.),
               np.array([1.0, 0.0, 0.0,-1.0])/np.sqrt(2.),
               np.array([0.0, 1.0,-1.0, 0.0])/np.sqrt(2.),
               np.array([0.0, 1.0, 0.0,-1.0])/np.sqrt(2.),
               np.array([0.0, 0.0, 1.0,-1.0])/np.sqrt(2.),
               #-----------------------------------------------
               np.array([0.0, 1.0,-1.0, 0.0])/np.sqrt(2.),
               np.array([0.0, 1.0, 0.0,-1.0])/np.sqrt(2.),
               np.array([0.0, 0.0, 1.0,-1.0])/np.sqrt(2.),
               #-----------------------------------------------
               np.array([0.0,-1.0,-1.0, 0.0])/np.sqrt(2.),
               np.array([0.0,-1.0, 0.0,-1.0])/np.sqrt(2.),
               np.array([0.0, 0.0,-1.0,-1.0])/np.sqrt(2.),
               #-----------------------------------------------
               np.array([1.0, 1.0, 1.0, 0.0])/np.sqrt(3.),
               np.array([1.0, 1.0, 0.0, 1.0])/np.sqrt(3.),
               np.array([1.0, 0.0, 1.0, 1.0])/np.sqrt(3.),
               np.array([1.0,-1.0, 1.0, 0.0])/np.sqrt(3.),
               np.array([1.0,-1.0, 0.0, 1.0])/np.sqrt(3.),
               np.array([1.0, 0.0,-1.0, 1.0])/np.sqrt(3.),
               np.array([1.0, 1.0,-1.0, 0.0])/np.sqrt(3.),
               np.array([1.0, 1.0, 0.0,-1.0])/np.sqrt(3.),
               np.array([1.0, 0.0, 1.0,-1.0])/np.sqrt(3.),
               np.array([1.0,-1.0,-1.0, 0.0])/np.sqrt(3.),
               np.array([1.0,-1.0, 0.0,-1.0])/np.sqrt(3.),
               np.array([1.0, 0.0,-1.0,-1.0])/np.sqrt(3.),
               #-----------------------------------------------
               np.array([0.0, 1.0, 1.0, 1.0])/np.sqrt(3.),
               np.array([0.0, 1.0,-1.0, 1.0])/np.sqrt(3.),
               np.array([0.0, 1.0, 1.0,-1.0])/np.sqrt(3.),
               np.array([0.0,-1.0, 1.0, 1.0])/np.sqrt(3.),
               np.array([0.0,-1.0,-1.0, 1.0])/np.sqrt(3.),
               np.array([0.0,-1.0, 1.0,-1.0])/np.sqrt(3.),
               np.array([0.0,-1.0,-1.0,-1.0])/np.sqrt(3.),
               #-----------------------------------------------
               np.array([1.0, 1.0, 1.0, 1.0])/2.,
               np.array([1.0,-1.0, 1.0, 1.0])/2.,
               np.array([1.0, 1.0,-1.0, 1.0])/2.,
               np.array([1.0, 1.0, 1.0,-1.0])/2.,
               np.array([1.0,-1.0,-1.0, 1.0])/2.,
               np.array([1.0,-1.0, 1.0,-1.0])/2.,
               np.array([1.0, 1.0,-1.0,-1.0])/2.,
               np.array([1.0,-1.0,-1.0,-1.0])/2.,
              ])
    specials_scatter = specials + np.broadcast_to(np.random.rand(4)*scatter,specials.shape)
    specials_scatter /= np.linalg.norm(specials_scatter,axis=1).reshape(-1,1)
    specials_scatter[specials_scatter[:,0]<0]*=-1

    return [Rotation.from_quaternion(s) for s in specials] + \
           [Rotation.from_quaternion(s) for s in specials_scatter] + \
           [Rotation.from_random() for _ in range(n-len(specials)-len(specials_scatter))]

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Rotation')


class TestRotation:

    @pytest.mark.parametrize('degrees',[True,False])
    def test_Eulers(self,default,degrees):
        for rot in default:
            m = rot.as_quaternion()
            o = Rotation.from_Eulers(rot.as_Eulers(degrees),degrees).as_quaternion()
            ok = np.allclose(m,o,atol=atol)
            if np.isclose(rot.as_quaternion()[0],0.0,atol=atol):
                ok = ok or np.allclose(m*-1.,o,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.isclose(np.linalg.norm(o),1.0)

    @pytest.mark.parametrize('P',[1,-1])
    @pytest.mark.parametrize('normalise',[True,False])
    @pytest.mark.parametrize('degrees',[True,False])
    def test_AxisAngle(self,default,degrees,normalise,P):
        c = np.array([P*-1,P*-1,P*-1,1.])
        for rot in default:
            m = rot.as_Eulers()
            o = Rotation.from_axis_angle(rot.as_axis_angle(degrees)*c,degrees,normalise,P).as_Eulers()
            u = np.array([np.pi*2,np.pi,np.pi*2])
            ok = np.allclose(m,o,atol=atol)
            ok = ok or np.allclose(np.where(np.isclose(m,u),m-u,m),np.where(np.isclose(o,u),o-u,o),atol=atol)
            if np.isclose(m[1],0.0,atol=atol) or np.isclose(m[1],np.pi,atol=atol):
                sum_phi = np.unwrap([m[0]+m[2],o[0]+o[2]])
                ok = ok or np.isclose(sum_phi[0],sum_phi[1],atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and (np.zeros(3)-1.e-9 <= o).all() and (o <= np.array([np.pi*2.,np.pi,np.pi*2.])+1.e-9).all()

    def test_Matrix(self,default):
        for rot in default:
            m = rot.as_axis_angle()
            o = Rotation.from_axis_angle(rot.as_axis_angle()).as_axis_angle()
            ok = np.allclose(m,o,atol=atol)
            if np.isclose(m[3],np.pi,atol=atol):
                ok = ok or np.allclose(m*np.array([-1.,-1.,-1.,1.]),o,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.isclose(np.linalg.norm(o[:3]),1.0) and o[3]<=np.pi++1.e-9

    @pytest.mark.parametrize('P',[1,-1])
    @pytest.mark.parametrize('normalise',[True,False])
    def test_Rodrigues(self,default,normalise,P):
        c = np.array([P*-1,P*-1,P*-1,1.])
        for rot in default:
            m = rot.as_matrix()
            o = Rotation.from_Rodrigues(rot.as_Rodrigues()*c,normalise,P).as_matrix()
            ok = np.allclose(m,o,atol=atol)
            print(m,o)
            assert ok and np.isclose(np.linalg.det(o),1.0)

    @pytest.mark.parametrize('P',[1,-1])
    def test_Homochoric(self,default,P):
        cutoff = np.tan(np.pi*.5*(1.-1e-4))
        for rot in default:
            m = rot.as_Rodrigues()
            o = Rotation.from_homochoric(rot.as_homochoric()*P*-1,P).as_Rodrigues()
            ok = np.allclose(np.clip(m,None,cutoff),np.clip(o,None,cutoff),atol=atol)
            ok = ok or np.isclose(m[3],0.0,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.isclose(np.linalg.norm(o[:3]),1.0)

    @pytest.mark.parametrize('P',[1,-1])
    def test_Cubochoric(self,default,P):
        for rot in default:
            m = rot.as_homochoric()
            o = Rotation.from_cubochoric(rot.as_cubochoric()*P*-1,P).as_homochoric()
            ok = np.allclose(m,o,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.linalg.norm(o) < (3.*np.pi/4.)**(1./3.) + 1.e-9

    @pytest.mark.parametrize('P',[1,-1])
    def test_Quaternion(self,default,P):
        c = np.array([1,P*-1,P*-1,P*-1])
        for rot in default:
            m = rot.as_cubochoric()
            o = Rotation.from_quaternion(rot.as_quaternion()*c,False,P).as_cubochoric()
            ok = np.allclose(m,o,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and o.max() < np.pi**(2./3.)*0.5+1.e-9

    @pytest.mark.parametrize('function',[Rotation.from_quaternion,
                                         Rotation.from_Eulers,
                                         Rotation.from_axis_angle,
                                         Rotation.from_matrix,
                                         Rotation.from_Rodrigues,
                                         Rotation.from_homochoric])
    def test_invalid_shape(self,function):
        invalid_shape = np.random.random(np.random.randint(8,32,(3)))
        with pytest.raises(ValueError):
            function(invalid_shape)

    @pytest.mark.parametrize('function,invalid',[(Rotation.from_quaternion, np.array([-1,0,0,0])),
                                                 (Rotation.from_quaternion, np.array([1,1,1,0])),
                                                 (Rotation.from_Eulers,     np.array([1,4,0])),
                                                 (Rotation.from_axis_angle, np.array([1,0,0,4])),
                                                 (Rotation.from_axis_angle, np.array([1,1,0,1])),
                                                 (Rotation.from_matrix,     np.random.rand(3,3)),
                                                 (Rotation.from_Rodrigues,  np.array([1,0,0,-1])),
                                                 (Rotation.from_Rodrigues,  np.array([1,1,0,1])),
                                                 (Rotation.from_homochoric, np.array([2,2,2]))  ])
    def test_invalid_value(self,function,invalid):
        with pytest.raises(ValueError):
            function(invalid)

    @pytest.mark.parametrize('vectorized, single',[(Rotation.qu2om,rotation_conversion.qu2om),
                                                   (Rotation.qu2eu,rotation_conversion.qu2eu),
                                                   (Rotation.qu2ax,rotation_conversion.qu2ax),
                                                   (Rotation.qu2ro,rotation_conversion.qu2ro),
                                                   (Rotation.qu2ho,rotation_conversion.qu2ho)])
    def test_quaternion_vectorization(self,default,vectorized,single):
        qu = np.array([rot.as_quaternion() for rot in default])
        vectorized(qu.reshape(qu.shape[0]//2,-1,4))
        co = vectorized(qu)
        for q,c in zip(qu,co):
            print(q,c)
            assert np.allclose(single(q),c) and np.allclose(single(q),vectorized(q))


    @pytest.mark.parametrize('vectorized, single',[(Rotation.om2qu,rotation_conversion.om2qu),
                                                   (Rotation.om2eu,rotation_conversion.om2eu),
                                                   (Rotation.om2ax,rotation_conversion.om2ax)])
    def test_matrix_vectorization(self,default,vectorized,single):
        om = np.array([rot.as_matrix() for rot in default])
        vectorized(om.reshape(om.shape[0]//2,-1,3,3))
        co = vectorized(om)
        for o,c in zip(om,co):
            print(o,c)
            assert np.allclose(single(o),c) and np.allclose(single(o),vectorized(o))

    @pytest.mark.parametrize('vectorized, single',[(Rotation.eu2qu,rotation_conversion.eu2qu),
                                                   (Rotation.eu2om,rotation_conversion.eu2om),
                                                   (Rotation.eu2ax,rotation_conversion.eu2ax),
                                                   (Rotation.eu2ro,rotation_conversion.eu2ro)])
    def test_Euler_vectorization(self,default,vectorized,single):
        eu = np.array([rot.as_Eulers() for rot in default])
        vectorized(eu.reshape(eu.shape[0]//2,-1,3))
        co = vectorized(eu)
        for e,c in zip(eu,co):
            print(e,c)
            assert np.allclose(single(e),c) and np.allclose(single(e),vectorized(e))

    @pytest.mark.parametrize('vectorized, single',[(Rotation.ax2qu,rotation_conversion.ax2qu),
                                                   (Rotation.ax2om,rotation_conversion.ax2om),
                                                   (Rotation.ax2ro,rotation_conversion.ax2ro),
                                                   (Rotation.ax2ho,rotation_conversion.ax2ho)])
    def test_axisAngle_vectorization(self,default,vectorized,single):
        ax = np.array([rot.as_axis_angle() for rot in default])
        vectorized(ax.reshape(ax.shape[0]//2,-1,4))
        co = vectorized(ax)
        for a,c in zip(ax,co):
            print(a,c)
            assert np.allclose(single(a),c) and np.allclose(single(a),vectorized(a))


    @pytest.mark.parametrize('vectorized, single',[(Rotation.ro2ax,rotation_conversion.ro2ax),
                                                   (Rotation.ro2ho,rotation_conversion.ro2ho)])
    def test_Rodrigues_vectorization(self,default,vectorized,single):
        ro = np.array([rot.as_Rodrigues() for rot in default])
        vectorized(ro.reshape(ro.shape[0]//2,-1,4))
        co = vectorized(ro)
        for r,c in zip(ro,co):
            print(r,c)
            assert np.allclose(single(r),c) and np.allclose(single(r),vectorized(r))

    @pytest.mark.parametrize('vectorized, single',[(Rotation.ho2ax,rotation_conversion.ho2ax),
                                                   (Rotation.ho2cu,rotation_conversion.ho2cu)])
    def test_homochoric_vectorization(self,default,vectorized,single):
        ho = np.array([rot.as_homochoric() for rot in default])
        vectorized(ho.reshape(ho.shape[0]//2,-1,3))
        co = vectorized(ho)
        for h,c in zip(ho,co):
            print(h,c)
            assert np.allclose(single(h),c) and np.allclose(single(h),vectorized(h))

    @pytest.mark.parametrize('vectorized, single',[(Rotation.cu2ho,rotation_conversion.cu2ho)])
    def test_cubochoric_vectorization(self,default,vectorized,single):
        cu = np.array([rot.as_cubochoric() for rot in default])
        vectorized(cu.reshape(cu.shape[0]//2,-1,3))
        co = vectorized(cu)
        for u,c in zip(cu,co):
            print(u,c)
            assert np.allclose(single(u),c) and np.allclose(single(u),vectorized(u))

    @pytest.mark.parametrize('direction',['forward',
                                          'backward'])
    def test_pyramid_vectorization(self,direction):
        p = np.random.rand(n,3)
        o = Rotation._get_pyramid_order(p,direction)
        for i,o_i in enumerate(o):
            assert np.all(o_i==Rotation._get_pyramid_order(p[i],direction))

    def test_pyramid_invariant(self):
        a = np.random.rand(n,3)
        f = Rotation._get_pyramid_order(a,'forward')
        b = Rotation._get_pyramid_order(a,'backward')
        assert np.all(np.take_along_axis(np.take_along_axis(a,f,-1),b,-1) == a)


    @pytest.mark.parametrize('data',[np.random.rand(3),
                                     np.random.rand(3,3),
                                     np.random.rand(3,3,3,3)])
    def test_rotate_identity(self,data):
        R = Rotation()
        assert np.allclose(data,R*data)

    @pytest.mark.parametrize('data',[np.random.rand(3),
                                     np.random.rand(3,3),
                                     np.random.rand(3,3,3,3)])
    def test_rotate_360deg(self,data):
        phi_1 = np.random.random() * np.pi
        phi_2 = 2*np.pi - phi_1
        R_1 = Rotation.from_Eulers(np.array([phi_1,0.,0.]))
        R_2 = Rotation.from_Eulers(np.array([0.,0.,phi_2]))
        assert np.allclose(data,R_2*(R_1*data))

    @pytest.mark.parametrize('data',[np.random.rand(3),
                                     np.random.rand(3,3),
                                     np.random.rand(3,3,3,3)])
    def test_rotate_inverse(self,data):
        R = Rotation.from_random()
        assert np.allclose(data,R.inversed()*(R*data))

    @pytest.mark.parametrize('data',[np.random.rand(4),
                                     np.random.rand(3,2),
                                     np.random.rand(3,2,3,3)])
    def test_rotate_invalid_shape(self,data):
        R = Rotation.from_random()
        with pytest.raises(ValueError):
            R*data

    @pytest.mark.parametrize('data',['does_not_work',
                                     (1,2),
                                     5])
    def test_rotate_invalid_type(self,data):
        R = Rotation.from_random()
        with pytest.raises(TypeError):
            R*data
