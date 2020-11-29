import pytest
import numpy as np

from damask import tensor
from damask import mechanics
from damask import Rotation

def stress_Cauchy(P,F):
    sigma = 1.0/np.linalg.det(F) * np.dot(P,F.T)
    return symmetric(sigma)


def eigenvalues(T_sym):
    return np.linalg.eigvalsh(symmetric(T_sym))


def maximum_shear(T_sym):
    w = eigenvalues(T_sym)
    return (w[0] - w[2])*0.5


def equivalent_strain_Mises(epsilon):
    return equivalent_Mises(epsilon,2.0/3.0)


def equivalent_stress_Mises(sigma):
    return equivalent_Mises(sigma,3.0/2.0)


def stress_second_Piola_Kirchhoff(P,F):
    S = np.dot(np.linalg.inv(F),P)
    return symmetric(S)


def rotation(T):
    return polar_decomposition(T,'R')[0]


def strain(F,t,m):

    if   t == 'V':
        B   = np.matmul(F,F.T)
        w,n = np.linalg.eigh(B)
    elif t == 'U':
        C   = np.matmul(F.T,F)
        w,n = np.linalg.eigh(C)

    if   m > 0.0:
        eps = 1.0/(2.0*abs(m)) * (+ np.matmul(n,np.einsum('j,kj->jk',w**m,n))
                                  - np.eye(3))
    elif m < 0.0:
        eps = 1.0/(2.0*abs(m)) * (- np.matmul(n,np.einsum('j,kj->jk',w**m,n))
                                  + np.eye(3))
    else:
        eps = np.matmul(n,np.einsum('j,kj->jk',0.5*np.log(w),n))

    return eps


def stretch_left(T):
    return polar_decomposition(T,'V')[0]


def stretch_right(T):
    return polar_decomposition(T,'U')[0]


def symmetric(T):
    return (T+T.T)*0.5


def polar_decomposition(T,requested):
    u, s, vh = np.linalg.svd(T)
    R = np.dot(u,vh)

    output = []
    if 'R' in requested:
        output.append(R)
    if 'V' in requested:
        output.append(np.dot(T,R.T))
    if 'U' in requested:
        output.append(np.dot(R.T,T))

    return tuple(output)

def equivalent_Mises(T_sym,s):
    return np.sqrt(s*(np.sum(deviatoric(T_sym)**2.0)))

def deviatoric(T):
    return T - np.eye(3)*np.trace(T)/3.0

class TestMechanics:

    n = 1000
    c = np.random.randint(n)

    @pytest.mark.parametrize('vectorized,single',[(mechanics.maximum_shear,           maximum_shear),
                                                  (mechanics.equivalent_stress_Mises, equivalent_stress_Mises),
                                                  (mechanics.equivalent_strain_Mises, equivalent_strain_Mises),
                                                  (mechanics.stretch_left,            stretch_left),
                                                  (mechanics.stretch_right,           stretch_right),
                                                 ])
    def test_vectorize_1_arg(self,vectorized,single):
        epsilon     = np.random.rand(self.n,3,3)
        epsilon_vec = np.reshape(epsilon,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(epsilon_vec),vectorized(epsilon).shape)):
            assert np.allclose(single(epsilon[i]),v)

    def test_vectorize_rotation(self):
        epsilon     = Rotation.from_random(self.n).as_matrix()
        epsilon_vec = np.reshape(epsilon,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(mechanics.rotation(epsilon_vec).as_matrix(),
                                        mechanics.rotation(epsilon).as_matrix().shape)):
            assert np.allclose(rotation(epsilon[i]),v)


    @pytest.mark.parametrize('vectorized,single',[(mechanics.stress_Cauchy,                 stress_Cauchy),
                                                  (mechanics.stress_second_Piola_Kirchhoff, stress_second_Piola_Kirchhoff)
                                                 ])
    def test_vectorize_2_arg(self,vectorized,single):
        P     = np.random.rand(self.n,3,3)
        F     = np.random.rand(self.n,3,3)
        P_vec = np.reshape(P,(self.n//10,10,3,3))
        F_vec = np.reshape(F,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(P_vec,F_vec),vectorized(P,F).shape)):
            assert np.allclose(single(P[i],F[i]),v)


    @pytest.mark.parametrize('vectorized,single',[(mechanics.strain,strain)])
    def test_vectorize_strain(self,vectorized,single):
        F     = np.random.rand(self.n,3,3)
        F_vec = np.reshape(F,(self.n//10,10,3,3))
        t     = ['V','U'][np.random.randint(0,2)]
        m     = np.random.random()*10.0 -5.0
        for i,v in enumerate(np.reshape(vectorized(F_vec,t,m),vectorized(F,t,m).shape)):
            assert np.allclose(single(F[i],t,m),v)

    @pytest.mark.parametrize('function',[mechanics.stress_Cauchy,
                                         mechanics.stress_second_Piola_Kirchhoff,
                                        ])
    def test_stress_measures(self,function):
        """Ensure that all stress measures are equivalent for no deformation."""
        P = np.random.rand(self.n,3,3)
        assert np.allclose(function(P,np.broadcast_to(np.eye(3),(self.n,3,3))),tensor.symmetric(P))

    def test_polar_decomposition(self):
        """F = RU = VR."""
        F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.rand(self.n,3,3)
        R = mechanics.rotation(F).as_matrix()
        V = mechanics.stretch_left(F)
        U = mechanics.stretch_right(F)
        assert np.allclose(np.matmul(R,U),
                           np.matmul(V,R))

    @pytest.mark.parametrize('m',[0.0,np.random.random()*10.,np.random.random()*-10.])
    def test_strain_no_rotation(self,m):
        """Ensure that left and right stretch give same results for no rotation."""
        F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.rand(self.n,3,3)
        assert np.allclose(mechanics.strain(F,'U',m),
                           mechanics.strain(F,'V',m))

    @pytest.mark.parametrize('m',[0.0,np.random.random()*2.5,np.random.random()*-2.5])
    def test_strain_rotation_equivalence(self,m):
        """Ensure that left and right strain differ only by a rotation."""
        F = np.broadcast_to(np.eye(3),[self.n,3,3]) + (np.random.rand(self.n,3,3)*0.5 - 0.25)
        assert np.allclose(np.linalg.det(mechanics.strain(F,'U',m)),
                           np.linalg.det(mechanics.strain(F,'V',m)))

    @pytest.mark.parametrize('m',[0.0,np.random.random(),np.random.random()*-1.])
    @pytest.mark.parametrize('t',['V','U'])
    def test_strain_rotation(self,m,t):
        """Ensure that pure rotation results in no strain."""
        F = Rotation.from_random(self.n).as_matrix()
        assert np.allclose(mechanics.strain(F,t,m),
                           0.0)

    def test_rotation_determinant(self):
        """
        Ensure that the determinant of the rotational part is +- 1.

        Should be +1, but random F might contain a reflection.
        """
        x = np.random.rand(self.n,3,3)
        assert np.allclose(np.abs(np.linalg.det(mechanics._polar_decomposition(x,'R')[0])),
                           1.0)

    def test_deviatoric_Mises(self):
        """Ensure that Mises equivalent stress depends only on deviatoric part."""
        x = np.random.rand(self.n,3,3)
        full = mechanics.equivalent_stress_Mises(x)
        dev  = mechanics.equivalent_stress_Mises(tensor.deviatoric(x))
        assert np.allclose(full,
                           dev)

    @pytest.mark.parametrize('Mises_equivalent',[mechanics.equivalent_strain_Mises,
                                                 mechanics.equivalent_stress_Mises])
    def test_spherical_Mises(self,Mises_equivalent):
        """Ensure that Mises equivalent strain/stress of spherical strain is 0."""
        x = np.random.rand(self.n,3,3)
        assert np.allclose(Mises_equivalent(tensor.spherical(x,True)),
                           0.0)


    def test_Mises(self):
        """Ensure that equivalent stress is 3/2 of equivalent strain."""
        x = np.random.rand(self.n,3,3)
        assert np.allclose(mechanics.equivalent_stress_Mises(x)/mechanics.equivalent_strain_Mises(x),
                           1.5)

    def test_spherical_no_shear(self):
        """Ensure that sherical stress has max shear of 0.0."""
        A = tensor.spherical(tensor.symmetric(np.random.rand(self.n,3,3)),True)
        assert np.allclose(mechanics.maximum_shear(A),0.0)

    def test_invalid_decomposition(self):
        with pytest.raises(ValueError):
            mechanics._polar_decomposition(np.random.rand(10,3,3),'A')
