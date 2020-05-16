import os

import pytest
import numpy as np

from damask import Rotation

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

    def test_Eulers(self,default):
        for rot in default:
            m = rot.as_quaternion()
            o = Rotation.from_Eulers(rot.as_Eulers()).as_quaternion()
            ok = np.allclose(m,o,atol=atol)
            if np.isclose(rot.as_quaternion()[0],0.0,atol=atol):
                ok = ok or np.allclose(m*-1.,o,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.isclose(np.linalg.norm(o),1.0)

    def test_AxisAngle(self,default):
        for rot in default:
            m = rot.as_Eulers()
            o = Rotation.from_axis_angle(rot.as_axis_angle()).as_Eulers()
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

    def test_Rodrigues(self,default):
        for rot in default:
            m = rot.as_matrix()
            o = Rotation.from_Rodrigues(rot.as_Rodrigues()).as_matrix()
            ok = np.allclose(m,o,atol=atol)
            print(m,o)
            assert ok and np.isclose(np.linalg.det(o),1.0)

    def test_Homochoric(self,default):
        cutoff = np.tan(np.pi*.5*(1.-1e-4))
        for rot in default:
            m = rot.as_Rodrigues()
            o = Rotation.from_homochoric(rot.as_homochoric()).as_Rodrigues()
            ok = np.allclose(np.clip(m,None,cutoff),np.clip(o,None,cutoff),atol=atol)
            ok = ok or np.isclose(m[3],0.0,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.isclose(np.linalg.norm(o[:3]),1.0)

    def test_Cubochoric(self,default):
        for rot in default:
            m = rot.as_homochoric()
            o = Rotation.from_cubochoric(rot.as_cubochoric()).as_homochoric()
            ok = np.allclose(m,o,atol=atol)
            print(m,o,rot.as_quaternion())
            assert ok and np.linalg.norm(o) < (3.*np.pi/4.)**(1./3.) + 1.e-9

    def test_Quaternion(self,default):
        for rot in default:
            m = rot.as_cubochoric()
            o = Rotation.from_quaternion(rot.as_quaternion()).as_cubochoric()
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

    @pytest.mark.parametrize('conversion',[Rotation.qu2om,
                                           Rotation.qu2eu,
                                           Rotation.qu2ax,
                                           Rotation.qu2ro,
                                           Rotation.qu2ho,
                                           Rotation.qu2cu
                                          ])
    def test_quaternion_vectorization(self,default,conversion):
        qu = np.array([rot.as_quaternion() for rot in default])
        conversion(qu.reshape(qu.shape[0]//2,-1,4))
        co = conversion(qu)
        for q,c in zip(qu,co):
            print(q,c)
            assert np.allclose(conversion(q),c)

    @pytest.mark.parametrize('conversion',[Rotation.om2qu,
                                           Rotation.om2eu,
                                           Rotation.om2ax,
                                           Rotation.om2ro,
                                           Rotation.om2ho,
                                           Rotation.om2cu
                                          ])
    def test_matrix_vectorization(self,default,conversion):
        om = np.array([rot.as_matrix() for rot in default])
        conversion(om.reshape(om.shape[0]//2,-1,3,3))
        co = conversion(om)
        for o,c in zip(om,co):
            print(o,c)
            assert np.allclose(conversion(o),c)

    @pytest.mark.parametrize('conversion',[Rotation.eu2qu,
                                           Rotation.eu2om,
                                           Rotation.eu2ax,
                                           Rotation.eu2ro,
                                           Rotation.eu2ho,
                                           Rotation.eu2cu
                                          ])
    def test_Euler_vectorization(self,default,conversion):
        eu = np.array([rot.as_Eulers() for rot in default])
        conversion(eu.reshape(eu.shape[0]//2,-1,3))
        co = conversion(eu)
        for e,c in zip(eu,co):
            print(e,c)
            assert np.allclose(conversion(e),c)

    @pytest.mark.parametrize('conversion',[Rotation.ax2qu,
                                           Rotation.ax2om,
                                           Rotation.ax2eu,
                                           Rotation.ax2ro,
                                           Rotation.ax2ho,
                                           Rotation.ax2cu
                                          ])
    def test_axisAngle_vectorization(self,default,conversion):
        ax = np.array([rot.as_axis_angle() for rot in default])
        conversion(ax.reshape(ax.shape[0]//2,-1,4))
        co = conversion(ax)
        for a,c in zip(ax,co):
            print(a,c)
            assert np.allclose(conversion(a),c)


    @pytest.mark.parametrize('conversion',[Rotation.ro2qu,
                                           Rotation.ro2om,
                                           Rotation.ro2eu,
                                           Rotation.ro2ax,
                                           Rotation.ro2ho,
                                           Rotation.ro2cu
                                          ])
    def test_Rodrigues_vectorization(self,default,conversion):
        ro = np.array([rot.as_Rodrigues() for rot in default])
        conversion(ro.reshape(ro.shape[0]//2,-1,4))
        co = conversion(ro)
        for r,c in zip(ro,co):
            print(r,c)
            assert np.allclose(conversion(r),c)

    @pytest.mark.parametrize('conversion',[Rotation.ho2qu,
                                           Rotation.ho2om,
                                           Rotation.ho2eu,
                                           Rotation.ho2ax,
                                           Rotation.ho2ro,
                                           Rotation.ho2cu
                                          ])
    def test_homochoric_vectorization(self,default,conversion):
        ho = np.array([rot.as_homochoric() for rot in default])
        conversion(ho.reshape(ho.shape[0]//2,-1,3))
        co = conversion(ho)
        for h,c in zip(ho,co):
            print(h,c)
            assert np.allclose(conversion(h),c)

    @pytest.mark.parametrize('conversion',[Rotation.cu2qu,
                                           Rotation.cu2om,
                                           Rotation.cu2eu,
                                           Rotation.cu2ax,
                                           Rotation.cu2ro,
                                           Rotation.cu2ho
                                          ])
    def test_cubochoric_vectorization(self,default,conversion):
        cu = np.array([rot.as_cubochoric() for rot in default])
        conversion(cu.reshape(cu.shape[0]//2,-1,3))
        co = conversion(cu)
        for u,c in zip(cu,co):
            print(u,c)
            assert np.allclose(conversion(u),c)

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
