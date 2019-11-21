import numpy as np
from damask import mechanics

class TestMechanics:
   
    n = 9
    c = np.random.randint(n)
   

    def test_vectorize_Cauchy(self):
         P = np.random.random((self.n,3,3))
         F = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Cauchy(F,P)[self.c],
                            mechanics.Cauchy(F[self.c],P[self.c]))


    def test_vectorize_strain_tensor(self):
         F = np.random.random((self.n,3,3))
         t = ['V','U'][np.random.randint(0,2)]
         m = np.random.random()*10. -5.0
         assert np.allclose(mechanics.strain_tensor(F,t,m)[self.c],
                            mechanics.strain_tensor(F[self.c],t,m))


    def test_vectorize_deviatoric_part(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.deviatoric_part(x)[self.c],
                            mechanics.deviatoric_part(x[self.c]))


    def test_vectorize_spherical_part(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.spherical_part(x)[self.c],
                            mechanics.spherical_part(x[self.c]))


    def test_vectorize_Mises_stress(self):
         sigma = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Mises_stress(sigma)[self.c],
                            mechanics.Mises_stress(sigma[self.c]))


    def test_vectorize_Mises_strain(self):
         epsilon = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Mises_strain(epsilon)[self.c],
                            mechanics.Mises_strain(epsilon[self.c]))


    def test_vectorize_symmetric(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.symmetric(x)[self.c],
                            mechanics.symmetric(x[self.c]))


    def test_vectorize_maximum_shear(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.maximum_shear(x)[self.c],
                            mechanics.maximum_shear(x[self.c]))


    def test_vectorize_principal_components(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.principal_components(x)[self.c],
                            mechanics.principal_components(x[self.c]))


    def test_vectorize_transpose(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.transpose(x)[self.c],
                            mechanics.transpose(x[self.c]))


    def test_vectorize_rotational_part(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.rotational_part(x)[self.c],
                            mechanics.rotational_part(x[self.c]))


    def test_vectorize_left_stretch(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.left_stretch(x)[self.c],
                            mechanics.left_stretch(x[self.c]))


    def test_vectorize_right_stretch(self):
         x = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.right_stretch(x)[self.c],
                            mechanics.right_stretch(x[self.c]))


    def test_Cauchy(self):
         """Ensure Cauchy stress is symmetrized 1. Piola-Kirchhoff stress for no deformation."""
         P = np.random.random((self.n,3,3))
         assert np.allclose(mechanics.Cauchy(np.broadcast_to(np.eye(3),(self.n,3,3)),P),
                            mechanics.symmetric(P))


    def test_strain_tensor_no_rotation(self):
         """Ensure that left and right stretch give same results for no rotation."""
         F = np.broadcast_to(np.eye(3),[self.n,3,3])*np.random.random((self.n,3,3))
         m = np.random.random()*20.0-10.0
         assert np.allclose(mechanics.strain_tensor(F,'U',m),
                            mechanics.strain_tensor(F,'V',m))
   

    def test_strain_tensor_rotation(self):
         """Ensure that pure rotation results in no strain."""
         F = mechanics.rotational_part(np.random.random((self.n,3,3)))
         t = ['V','U'][np.random.randint(0,2)]
         m = np.random.random()*2.0 - 1.0
         assert np.allclose(mechanics.strain_tensor(F,t,m),
                            0.0)
   

    def test_spherical_deviatoric_part(self):
         """Ensure that full tensor is sum of spherical and deviatoric part."""
         x = np.random.random((self.n,3,3))
         sph = np.broadcast_to(np.eye(3),(self.n,3,3))\
             * np.repeat(mechanics.spherical_part(x),9).reshape(self.n,3,3)
         assert np.allclose(sph + mechanics.deviatoric_part(x),
                            x)


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
