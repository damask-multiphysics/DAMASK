import numpy as np
from damask import mechanics

class TestMechanics:

    n = 1000
    c = np.random.randint(n)


    def test_vectorize_Cauchy(self):
         P = np.random.random((self.n,3,3))
         F = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Cauchy(F,P)[self.c],
                            mechanics.Cauchy(F[self.c],P[self.c]))

    def test_vectorize_deviatoric_part(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.deviatoric_part(x)[self.c],
                            mechanics.deviatoric_part(x[self.c]))

    def test_vectorize_eigenvalues(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.eigenvalues(x)[self.c],
                            mechanics.eigenvalues(x[self.c]))

    def test_vectorize_eigenvectors(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.eigenvectors(x)[self.c],
                            mechanics.eigenvectors(x[self.c]))

    def test_vectorize_left_stretch(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.left_stretch(x)[self.c],
                            mechanics.left_stretch(x[self.c]))

    def test_vectorize_maximum_shear(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.maximum_shear(x)[self.c],
                            mechanics.maximum_shear(x[self.c]))

    def test_vectorize_Mises_strain(self):
         epsilon = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Mises_strain(epsilon)[self.c],
                            mechanics.Mises_strain(epsilon[self.c]))

    def test_vectorize_Mises_stress(self):
         sigma = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Mises_stress(sigma)[self.c],
                            mechanics.Mises_stress(sigma[self.c]))

    def test_vectorize_PK2(self):
         F = np.random.random((self.n,3,3))
         P = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.PK2(F,P)[self.c],
                            mechanics.PK2(F[self.c],P[self.c]))

    def test_vectorize_right_stretch(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.right_stretch(x)[self.c],
                            mechanics.right_stretch(x[self.c]))

    def test_vectorize_rotational_part(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.rotational_part(x)[self.c],
                            mechanics.rotational_part(x[self.c]))

    def test_vectorize_spherical_part(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.spherical_part(x,True)[self.c],
                            mechanics.spherical_part(x[self.c],True))

    def test_vectorize_strain_tensor(self):
         F = np.random.random((self.n,3,3))
         t = ['V','U'][np.random.randint(0,2)]
         m = np.random.random()*10. -5.0
         assert np.allclose(mechanics.strain_tensor(F,t,m)[self.c],
                            mechanics.strain_tensor(F[self.c],t,m))

    def test_vectorize_symmetric(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.symmetric(x)[self.c],
                            mechanics.symmetric(x[self.c]))

    def test_vectorize_transpose(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.transpose(x)[self.c],
                            mechanics.transpose(x[self.c]))


    def test_Cauchy(self):
         """Ensure Cauchy stress is symmetrized 1. Piola-Kirchhoff stress for no deformation."""
         P = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Cauchy(np.broadcast_to(np.eye(3),(self.n,3,3)),P),
                            mechanics.symmetric(P))


    def test_polar_decomposition(self):
         """F = RU = VR."""
         F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.random((self.n,3,3))
         R = mechanics.rotational_part(F)
         V = mechanics.left_stretch(F)
         U = mechanics.right_stretch(F)
         assert np.allclose(np.matmul(R,U),
                            np.matmul(V,R))


    def test_PK2(self):
         """Ensure 2. Piola-Kirchhoff stress is symmetrized 1. Piola-Kirchhoff stress for no deformation."""
         P = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.PK2(np.broadcast_to(np.eye(3),(self.n,3,3)),P),
                            mechanics.symmetric(P))


    def test_strain_tensor_no_rotation(self):
         """Ensure that left and right stretch give same results for no rotation."""
         F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.random((self.n,3,3))
         m = np.random.random()*20.0-10.0
         assert np.allclose(mechanics.strain_tensor(F,'U',m),
                            mechanics.strain_tensor(F,'V',m))

    def test_strain_tensor_rotation_equivalence(self):
         """Ensure that left and right strain differ only by a rotation."""
         F = np.broadcast_to(np.eye(3),[self.n,3,3]) + (np.random.random((self.n,3,3))*0.5 - 0.25)
         m = np.random.random()*5.0-2.5
         assert np.allclose(np.linalg.det(mechanics.strain_tensor(F,'U',m)),
                            np.linalg.det(mechanics.strain_tensor(F,'V',m)))

    def test_strain_tensor_rotation(self):
         """Ensure that pure rotation results in no strain."""
         F = mechanics.rotational_part(np.random.random((self.n,3,3)))
         t = ['V','U'][np.random.randint(0,2)]
         m = np.random.random()*2.0 - 1.0
         assert np.allclose(mechanics.strain_tensor(F,t,m),
                            0.0)

    def test_rotation_determinant(self):
         """
         Ensure that the determinant of the rotational part is +- 1.

         Should be +1, but random F might contain a reflection.
         """
         x = np.random.random((self.n,3,3))
         assert np.allclose(np.abs(np.linalg.det(mechanics.rotational_part(x))),
                            1.0)


    def test_spherical_deviatoric_part(self):
         """Ensure that full tensor is sum of spherical and deviatoric part."""
         x = np.random.random((self.n,3,3))
         sph = mechanics.spherical_part(x,True)
         assert np.allclose(sph + mechanics.deviatoric_part(x),
                            x)

    def test_deviatoric_Mises(self):
         """Ensure that Mises equivalent stress depends only on deviatoric part."""
         x = np.random.random((self.n,3,3))
         full = mechanics.Mises_stress(x)
         dev  = mechanics.Mises_stress(mechanics.deviatoric_part(x))
         assert np.allclose(full,
                            dev)

    def test_spherical_mapping(self):
         """Ensure that mapping to tensor is correct."""
         x = np.random.random((self.n,3,3))
         tensor = mechanics.spherical_part(x,True)
         scalar = mechanics.spherical_part(x)
         assert np.allclose(np.linalg.det(tensor),
                            scalar**3.0)

    def test_spherical_Mises(self):
         """Ensure that Mises equivalent strrain of spherical strain is 0."""
         x = np.random.random((self.n,3,3))
         sph = mechanics.spherical_part(x,True)
         assert np.allclose(mechanics.Mises_strain(sph),
                            0.0)

    def test_symmetric(self):
         """Ensure that a symmetric tensor is half of the sum of a tensor and its transpose."""
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.symmetric(x)*2.0,
                            mechanics.transpose(x)+x)


    def test_transpose(self):
         """Ensure that a symmetric tensor equals its transpose."""
         x = mechanics.symmetric(np.random.random((self.n,3,3)))
         assert np.allclose(mechanics.transpose(x),
                            x)


    def test_Mises(self):
         """Ensure that equivalent stress is 3/2 of equivalent strain."""
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Mises_stress(x)/mechanics.Mises_strain(x),
                            1.5)


    def test_eigenvalues(self):
         """Ensure that the characteristic polynomial can be solved."""
         A = mechanics.symmetric(np.random.random((self.n,3,3)))
         lambd = mechanics.eigenvalues(A)
         s = np.random.randint(self.n)
         for i in range(3):
            assert np.allclose(np.linalg.det(A[s]-lambd[s,i]*np.eye(3)),.0)

    def test_eigenvalues_and_vectors(self):
         """Ensure that eigenvalues and -vectors are the solution to the characteristic polynomial."""
         A = mechanics.symmetric(np.random.random((self.n,3,3)))
         lambd = mechanics.eigenvalues(A)
         x     = mechanics.eigenvectors(A)
         s = np.random.randint(self.n)
         for i in range(3):
            assert np.allclose(np.dot(A[s]-lambd[s,i]*np.eye(3),x[s,:,i]),.0)

    def test_eigenvectors_RHS(self):
         """Ensure that RHS coordinate system does only change sign of determinant."""
         A = mechanics.symmetric(np.random.random((self.n,3,3)))
         LRHS = np.linalg.det(mechanics.eigenvectors(A,RHS=False))
         RHS  = np.linalg.det(mechanics.eigenvectors(A,RHS=True))
         s = np.random.randint(self.n)
         assert np.allclose(np.abs(LRHS),RHS)

    def test_spherical_no_shear(self):
         """Ensure that sherical stress has max shear of 0.0."""
         A = mechanics.spherical_part(mechanics.symmetric(np.random.random((self.n,3,3))),True)
         assert np.allclose(mechanics.maximum_shear(A),0.0)
