import pytest
import numpy as np

from damask import tensor

def deviatoric(T):
    return T - spherical(T)

def eigenvalues(T_sym):
    return np.linalg.eigvalsh(symmetric(T_sym))

def eigenvectors(T_sym,RHS=False):
    (u,v) = np.linalg.eigh(symmetric(T_sym))

    if RHS:
        if np.linalg.det(v) < 0.0: v[:,2] *= -1.0
    return v

def symmetric(T):
    return (T+transpose(T))*0.5

def transpose(T):
    return T.T

def spherical(T,tensor=True):
    sph = np.trace(T)/3.0
    return sph if not tensor else np.eye(3)*sph


class TestTensor:

    n = 1000
    c = np.random.randint(n)


    @pytest.mark.parametrize('vectorized,single',[(tensor.deviatoric,   deviatoric),
                                                  (tensor.eigenvalues,  eigenvalues),
                                                  (tensor.eigenvectors, eigenvectors),
                                                  (tensor.symmetric,    symmetric),
                                                  (tensor.transpose,    transpose),
                                                  (tensor.spherical,    spherical),
                                        ])
    def test_vectorize_1_arg(self,vectorized,single):
        epsilon     = np.random.rand(self.n,3,3)
        epsilon_vec = np.reshape(epsilon,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(epsilon_vec),vectorized(epsilon).shape)):
            assert np.allclose(single(epsilon[i]),v)

    def test_symmetric(self):
        """Ensure that a symmetric tensor is half of the sum of a tensor and its transpose."""
        x = np.random.rand(self.n,3,3)
        assert np.allclose(tensor.symmetric(x)*2.0,tensor.transpose(x)+x)

    def test_transpose(self):
        """Ensure that a symmetric tensor equals its transpose."""
        x = tensor.symmetric(np.random.rand(self.n,3,3))
        assert np.allclose(tensor.transpose(x),x)

    def test_eigenvalues(self):
        """Ensure that the characteristic polynomial can be solved."""
        A = tensor.symmetric(np.random.rand(self.n,3,3))
        lambd = tensor.eigenvalues(A)
        s = np.random.randint(self.n)
        for i in range(3):
           assert np.allclose(np.linalg.det(A[s]-lambd[s,i]*np.eye(3)),.0)

    def test_eigenvalues_and_vectors(self):
        """Ensure that eigenvalues and -vectors are the solution to the characteristic polynomial."""
        A = tensor.symmetric(np.random.rand(self.n,3,3))
        lambd = tensor.eigenvalues(A)
        x     = tensor.eigenvectors(A)
        s = np.random.randint(self.n)
        for i in range(3):
           assert np.allclose(np.dot(A[s]-lambd[s,i]*np.eye(3),x[s,:,i]),.0)

    def test_eigenvectors_RHS(self):
        """Ensure that RHS coordinate system does only change sign of determinant."""
        A = tensor.symmetric(np.random.rand(self.n,3,3))
        LRHS = np.linalg.det(tensor.eigenvectors(A,RHS=False))
        RHS  = np.linalg.det(tensor.eigenvectors(A,RHS=True))
        assert np.allclose(np.abs(LRHS),RHS)

    def test_spherical_deviatoric_part(self):
        """Ensure that full tensor is sum of spherical and deviatoric part."""
        x = np.random.rand(self.n,3,3)
        assert np.allclose(tensor.spherical(x,True) + tensor.deviatoric(x),
                           x)
    def test_spherical_mapping(self):
        """Ensure that mapping to tensor is correct."""
        x = np.random.rand(self.n,3,3)
        tnsr   = tensor.spherical(x,True)
        scalar = tensor.spherical(x,False)
        assert np.allclose(np.linalg.det(tnsr),
                           scalar**3.0)

    def test_deviatoric(self):
        I_n = np.broadcast_to(np.eye(3),(self.n,3,3))
        r   = np.logical_not(I_n)*np.random.rand(self.n,3,3)
        assert np.allclose(tensor.deviatoric(I_n+r),r)
