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

    return [Rotation.fromQuaternion(s) for s in specials] + \
           [Rotation.fromQuaternion(s) for s in specials_scatter] + \
           [Rotation.fromRandom() for _ in range(n-len(specials)-len(specials_scatter))]

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Rotation')


class TestRotation:

    def test_Eulers(self,default):
        for rot in default:
            m = rot.asQuaternion()
            o = Rotation.fromEulers(rot.asEulers()).asQuaternion()
            ok = np.allclose(m,o,atol=atol)
            if np.isclose(rot.asQuaternion()[0],0.0,atol=atol):
                ok = ok or np.allclose(m*-1.,o,atol=atol)
            print(m,o,rot.asQuaternion())
            assert ok and np.isclose(np.linalg.norm(o),1.0)

    def test_AxisAngle(self,default):
        for rot in default:
            m = rot.asEulers()
            o = Rotation.fromAxisAngle(rot.asAxisAngle()).asEulers()
            u = np.array([np.pi*2,np.pi,np.pi*2])
            ok = np.allclose(m,o,atol=atol)
            ok = ok or np.allclose(np.where(np.isclose(m,u),m-u,m),np.where(np.isclose(o,u),o-u,o),atol=atol)
            if np.isclose(m[1],0.0,atol=atol) or np.isclose(m[1],np.pi,atol=atol):
                sum_phi = np.unwrap([m[0]+m[2],o[0]+o[2]])
                ok = ok or np.isclose(sum_phi[0],sum_phi[1],atol=atol)
            print(m,o,rot.asQuaternion())
            assert ok and (np.zeros(3)-1.e-9 <= o).all() and (o <= np.array([np.pi*2.,np.pi,np.pi*2.])+1.e-9).all()

    def test_Matrix(self,default):
        for rot in default:
            m = rot.asAxisAngle()
            o = Rotation.fromAxisAngle(rot.asAxisAngle()).asAxisAngle()
            ok = np.allclose(m,o,atol=atol)
            if np.isclose(m[3],np.pi,atol=atol):
                ok = ok or np.allclose(m*np.array([-1.,-1.,-1.,1.]),o,atol=atol)
            print(m,o,rot.asQuaternion())
            assert ok and np.isclose(np.linalg.norm(o[:3]),1.0) and o[3]<=np.pi++1.e-9

    def test_Rodriques(self,default):
        for rot in default:
            m = rot.asMatrix()
            o = Rotation.fromRodrigues(rot.asRodrigues()).asMatrix()
            ok = np.allclose(m,o,atol=atol)
            print(m,o)
            assert ok and np.isclose(np.linalg.det(o),1.0)

    def test_Homochoric(self,default):
        cutoff = np.tan(np.pi*.5*(1.-1e-4))
        for rot in default:
            m = rot.asRodrigues()
            o = Rotation.fromHomochoric(rot.asHomochoric()).asRodrigues()
            ok = np.allclose(np.clip(m,None,cutoff),np.clip(o,None,cutoff),atol=atol)
            ok = ok or np.isclose(m[3],0.0,atol=atol)
            print(m,o,rot.asQuaternion())
            assert ok and np.isclose(np.linalg.norm(o[:3]),1.0)

    def test_Cubochoric(self,default):
        for rot in default:
            m = rot.asHomochoric()
            o = Rotation.fromCubochoric(rot.asCubochoric()).asHomochoric()
            ok = np.allclose(m,o,atol=atol)
            print(m,o,rot.asQuaternion())
            assert ok and np.linalg.norm(o) < (3.*np.pi/4.)**(1./3.) + 1.e-9

    def test_Quaternion(self,default):
        for rot in default:
            m = rot.asCubochoric()
            o = Rotation.fromQuaternion(rot.asQuaternion()).asCubochoric()
            ok = np.allclose(m,o,atol=atol)
            print(m,o,rot.asQuaternion())
            assert ok and o.max() < np.pi**(2./3.)*0.5+1.e-9

    @pytest.mark.parametrize('conversion',[Rotation.qu2om,
                                           Rotation.qu2eu,
                                           Rotation.qu2ax,
                                           Rotation.qu2ro,
                                           Rotation.qu2ho])
    def test_quaternion_vectorization(self,default,conversion):
        qu = np.array([rot.asQuaternion() for rot in default])
        dev_null = conversion(qu.reshape(qu.shape[0]//2,-1,4))
        co = conversion(qu)
        for q,c in zip(qu,co):
            print(q,c)
            assert np.allclose(conversion(q),c)

    @pytest.mark.parametrize('conversion',[Rotation.om2eu,
                                           Rotation.om2ax,
                                          ])
    def test_matrix_vectorization(self,default,conversion):
        om = np.array([rot.asMatrix() for rot in default])
        dev_null = conversion(om.reshape(om.shape[0]//2,-1,3,3))
        co = conversion(om)
        for o,c in zip(om,co):
            print(o,c)
            assert np.allclose(conversion(o),c)

    @pytest.mark.parametrize('conversion',[Rotation.eu2qu,
                                           Rotation.eu2om,
                                           Rotation.eu2ax,
                                           Rotation.eu2ro,
                                          ])
    def test_Euler_vectorization(self,default,conversion):
        eu = np.array([rot.asEulers() for rot in default])
        dev_null = conversion(eu.reshape(eu.shape[0]//2,-1,3))
        co = conversion(eu)
        for e,c in zip(eu,co):
            print(e,c)
            assert np.allclose(conversion(e),c)

    @pytest.mark.parametrize('conversion',[Rotation.ax2qu,
                                           Rotation.ax2om,
                                           Rotation.ax2ro,
                                           Rotation.ax2ho,
                                          ])
    def test_axisAngle_vectorization(self,default,conversion):
        ax = np.array([rot.asAxisAngle() for rot in default])
        dev_null = conversion(ax.reshape(ax.shape[0]//2,-1,4))
        co = conversion(ax)
        for a,c in zip(ax,co):
            print(a,c)
            assert np.allclose(conversion(a),c)
