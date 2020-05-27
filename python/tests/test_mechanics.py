import pytest
import numpy as np

from damask import mechanics

def deviatoric_part(T):
    return T - np.eye(3)*spherical_part(T)

def spherical_part(T,tensor=False):
    sph = np.trace(T)/3.0
    return sph if not tensor else np.eye(3)*sph

class TestMechanics:

    n = 1000
    c = np.random.randint(n)


    @pytest.mark.parametrize('vectorized,single',[(mechanics.deviatoric_part, deviatoric_part),
                                                  (mechanics.spherical_part,  spherical_part)
                                                 ])
    def test_vectorize_1_arg_(self,vectorized,single):
        test_data_flat = np.random.rand(self.n,3,3)
        test_data = np.reshape(test_data_flat,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(test_data),vectorized(test_data_flat).shape)):
            assert np.allclose(single(test_data_flat[i]),v)

    @pytest.mark.parametrize('function',[mechanics.deviatoric_part,
                                         mechanics.eigenvalues,
                                         mechanics.eigenvectors,
                                         mechanics.left_stretch,
                                         mechanics.maximum_shear,
                                         mechanics.Mises_strain,
                                         mechanics.Mises_stress,
                                         mechanics.right_stretch,
                                         mechanics.rotational_part,
                                         mechanics.spherical_part,
                                         mechanics.symmetric,
                                         mechanics.transpose,
                                        ])
    def test_vectorize_1_arg(self,function):
        epsilon = np.random.rand(self.n,3,3)
        assert np.allclose(function(epsilon)[self.c],function(epsilon[self.c]))

    @pytest.mark.parametrize('function',[mechanics.Cauchy,
                                         mechanics.PK2,
                                        ])
    def test_vectorize_2_arg(self,function):
        P = np.random.rand(self.n,3,3)
        F = np.random.rand(self.n,3,3)
        assert np.allclose(function(P,F)[self.c],function(P[self.c],F[self.c]))

    def test_vectorize_strain_tensor(self):
        F = np.random.rand(self.n,3,3)
        t = ['V','U'][np.random.randint(0,2)]
        m = np.random.random()*10. -5.0
        assert np.allclose(mechanics.strain_tensor(F,t,m)[self.c],
                           mechanics.strain_tensor(F[self.c],t,m))

    @pytest.mark.parametrize('function',[mechanics.Cauchy,
                                         mechanics.PK2,
                                        ])
    def test_stress_measures(self,function):
        """Ensure that all stress measures are equivalent for no deformation."""
        P = np.random.rand(self.n,3,3)
        assert np.allclose(function(P,np.broadcast_to(np.eye(3),(self.n,3,3))),mechanics.symmetric(P))

    def test_deviatoric_part(self):
        I_n = np.broadcast_to(np.eye(3),(self.n,3,3))
        r   = np.logical_not(I_n)*np.random.rand(self.n,3,3)
        assert np.allclose(mechanics.deviatoric_part(I_n+r),r)

    def test_polar_decomposition(self):
        """F = RU = VR."""
        F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.rand(self.n,3,3)
        R = mechanics.rotational_part(F)
        V = mechanics.left_stretch(F)
        U = mechanics.right_stretch(F)
        assert np.allclose(np.matmul(R,U),
                           np.matmul(V,R))

    def test_strain_tensor_no_rotation(self):
        """Ensure that left and right stretch give same results for no rotation."""
        F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.rand(self.n,3,3)
        m = np.random.random()*20.0-10.0
        assert np.allclose(mechanics.strain_tensor(F,'U',m),
                           mechanics.strain_tensor(F,'V',m))

    def test_strain_tensor_rotation_equivalence(self):
        """Ensure that left and right strain differ only by a rotation."""
        F = np.broadcast_to(np.eye(3),[self.n,3,3]) + (np.random.rand(self.n,3,3)*0.5 - 0.25)
        m = np.random.random()*5.0-2.5
        assert np.allclose(np.linalg.det(mechanics.strain_tensor(F,'U',m)),
                           np.linalg.det(mechanics.strain_tensor(F,'V',m)))

    def test_strain_tensor_rotation(self):
        """Ensure that pure rotation results in no strain."""
        F = mechanics.rotational_part(np.random.rand(self.n,3,3))
        t = ['V','U'][np.random.randint(0,2)]
        m = np.random.random()*2.0 - 1.0
        assert np.allclose(mechanics.strain_tensor(F,t,m),
                           0.0)

    def test_rotation_determinant(self):
        """
        Ensure that the determinant of the rotational part is +- 1.

        Should be +1, but random F might contain a reflection.
        """
        x = np.random.rand(self.n,3,3)
        assert np.allclose(np.abs(np.linalg.det(mechanics.rotational_part(x))),
                           1.0)

    def test_spherical_deviatoric_part(self):
        """Ensure that full tensor is sum of spherical and deviatoric part."""
        x = np.random.rand(self.n,3,3)
        sph = mechanics.spherical_part(x,True)
        assert np.allclose(sph + mechanics.deviatoric_part(x),
                           x)

    def test_deviatoric_Mises(self):
        """Ensure that Mises equivalent stress depends only on deviatoric part."""
        x = np.random.rand(self.n,3,3)
        full = mechanics.Mises_stress(x)
        dev  = mechanics.Mises_stress(mechanics.deviatoric_part(x))
        assert np.allclose(full,
                           dev)

    def test_spherical_mapping(self):
        """Ensure that mapping to tensor is correct."""
        x = np.random.rand(self.n,3,3)
        tensor = mechanics.spherical_part(x,True)
        scalar = mechanics.spherical_part(x)
        assert np.allclose(np.linalg.det(tensor),
                           scalar**3.0)

    def test_spherical_Mises(self):
        """Ensure that Mises equivalent strrain of spherical strain is 0."""
        x = np.random.rand(self.n,3,3)
        sph = mechanics.spherical_part(x,True)
        assert np.allclose(mechanics.Mises_strain(sph),
                           0.0)

    def test_symmetric(self):
        """Ensure that a symmetric tensor is half of the sum of a tensor and its transpose."""
        x = np.random.rand(self.n,3,3)
        assert np.allclose(mechanics.symmetric(x)*2.0,
                           mechanics.transpose(x)+x)

    def test_transpose(self):
        """Ensure that a symmetric tensor equals its transpose."""
        x = mechanics.symmetric(np.random.rand(self.n,3,3))
        assert np.allclose(mechanics.transpose(x),
                           x)

    def test_Mises(self):
        """Ensure that equivalent stress is 3/2 of equivalent strain."""
        x = np.random.rand(self.n,3,3)
        assert np.allclose(mechanics.Mises_stress(x)/mechanics.Mises_strain(x),
                           1.5)

    def test_eigenvalues(self):
        """Ensure that the characteristic polynomial can be solved."""
        A = mechanics.symmetric(np.random.rand(self.n,3,3))
        lambd = mechanics.eigenvalues(A)
        s = np.random.randint(self.n)
        for i in range(3):
           assert np.allclose(np.linalg.det(A[s]-lambd[s,i]*np.eye(3)),.0)

    def test_eigenvalues_and_vectors(self):
        """Ensure that eigenvalues and -vectors are the solution to the characteristic polynomial."""
        A = mechanics.symmetric(np.random.rand(self.n,3,3))
        lambd = mechanics.eigenvalues(A)
        x     = mechanics.eigenvectors(A)
        s = np.random.randint(self.n)
        for i in range(3):
           assert np.allclose(np.dot(A[s]-lambd[s,i]*np.eye(3),x[s,:,i]),.0)

    def test_eigenvectors_RHS(self):
        """Ensure that RHS coordinate system does only change sign of determinant."""
        A = mechanics.symmetric(np.random.rand(self.n,3,3))
        LRHS = np.linalg.det(mechanics.eigenvectors(A,RHS=False))
        RHS  = np.linalg.det(mechanics.eigenvectors(A,RHS=True))
        assert np.allclose(np.abs(LRHS),RHS)

    def test_spherical_no_shear(self):
        """Ensure that sherical stress has max shear of 0.0."""
        A = mechanics.spherical_part(mechanics.symmetric(np.random.rand(self.n,3,3)),True)
        assert np.allclose(mechanics.maximum_shear(A),0.0)
