import pytest
import numpy as np

from damask import mechanics

def Cauchy(P,F):
    sigma = 1.0/np.linalg.det(F) * np.dot(P,F.T)
    return mechanics.symmetric(sigma)


def deviatoric_part(T):
    return T - np.eye(3)*spherical_part(T)


def eigenvalues(T_sym):
    return np.linalg.eigvalsh(symmetric(T_sym))


def eigenvectors(T_sym,RHS=False):
    (u,v) = np.linalg.eigh(symmetric(T_sym))

    if RHS:
        if np.linalg.det(v) < 0.0: v[:,2] *= -1.0
    return v


def left_stretch(T):
    return polar_decomposition(T,'V')[0]


def maximum_shear(T_sym):
    w = eigenvalues(T_sym)
    return (w[0] - w[2])*0.5


def Mises_strain(epsilon):
    return Mises(epsilon,2.0/3.0)


def Mises_stress(sigma):
    return Mises(sigma,3.0/2.0)


def PK2(P,F):
    S = np.dot(np.linalg.inv(F),P)
    return symmetric(S)


def right_stretch(T):
    return polar_decomposition(T,'U')[0]


def rotational_part(T):
    return polar_decomposition(T,'R')[0]

def spherical_part(T,tensor=False):
    sph = np.trace(T)/3.0
    return sph if not tensor else np.eye(3)*sph


def strain_tensor(F,t,m):
    F_  = F.reshape(1,3,3)
    
    if   t == 'V':
        B   = np.matmul(F_,transpose(F_))
        w,n = np.linalg.eigh(B)
    elif t == 'U':
        C   = np.matmul(transpose(F_),F_)
        w,n = np.linalg.eigh(C)

    if   m > 0.0:
        eps = 1.0/(2.0*abs(m)) * (+ np.matmul(n,np.einsum('ij,ikj->ijk',w**m,n))
                                  - np.einsum('ijk->ijk',np.eye(3)))
    elif m < 0.0:
        eps = 1.0/(2.0*abs(m)) * (- np.matmul(n,np.einsum('ij,ikj->ijk',w**m,n))
                                  + np.einsum('ijk->ijk',np.eye(3)))
    else:
        eps = np.matmul(n,np.einsum('ij,ikj->ijk',0.5*np.log(w),n))

    return eps.reshape(3,3)


def symmetric(T):
    return (T+transpose(T))*0.5


def transpose(T):
    return T.T


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

def Mises(T_sym,s):
    d = deviatoric_part(T_sym)
    return np.sqrt(s*(np.sum(d**2.0)))


class TestMechanics:

    n = 1000
    c = np.random.randint(n)


    @pytest.mark.parametrize('vectorized,single',[(mechanics.deviatoric_part, deviatoric_part),
                                                  (mechanics.spherical_part,  spherical_part)
                                                 ])
    def test_vectorize_1_arg_(self,vectorized,single):
        print("done")
        test_data_flat = np.random.rand(self.n,3,3)
        test_data = np.reshape(test_data_flat,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(test_data),vectorized(test_data_flat).shape)):
            assert np.allclose(single(test_data_flat[i]),v)

    @pytest.mark.parametrize('vectorized,single',[
        (mechanics.deviatoric_part, deviatoric_part),
        (mechanics.eigenvalues    , eigenvalues    ),
        (mechanics.eigenvectors   , eigenvectors   ),
        (mechanics.left_stretch   , left_stretch   ),
        (mechanics.maximum_shear  , maximum_shear  ),
        (mechanics.Mises_strain   , Mises_strain   ),
        (mechanics.Mises_stress   , Mises_stress   ),
        (mechanics.right_stretch  , right_stretch  ),
        (mechanics.rotational_part, rotational_part),
        (mechanics.spherical_part , spherical_part ),
        (mechanics.symmetric      , symmetric      ),
        (mechanics.transpose      , transpose      ),
                                        ])
    def test_vectorize_1_arg(self,vectorized,single):
        epsilon     = np.random.rand(self.n,3,3)
        epsilon_vec = np.reshape(epsilon,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(epsilon_vec),vectorized(epsilon).shape)):
            assert np.allclose(single(epsilon[i]),v)

    @pytest.mark.parametrize('vectorized,single',[
        (mechanics.Cauchy,Cauchy),
        (mechanics.PK2   ,PK2   )
                                                 ])
    def test_vectorize_2_arg(self,vectorized,single):
        P     = np.random.rand(self.n,3,3)
        F     = np.random.rand(self.n,3,3)
        P_vec = np.reshape(P,(self.n//10,10,3,3))
        F_vec = np.reshape(F,(self.n//10,10,3,3))
        for i,v in enumerate(np.reshape(vectorized(P_vec,F_vec),vectorized(P,F).shape)):
            assert np.allclose(single(P[i],F[i]),v)


    @pytest.mark.parametrize('vectorized,single',[(mechanics.strain_tensor,strain_tensor)])
    def test_vectorize_strain_tensor(self,vectorized,single):
        F     = np.random.rand(self.n,3,3)
        F_vec = np.reshape(F,(self.n//10,10,3,3))
        t     = ['V','U'][np.random.randint(0,2)]
        m     = np.random.random()*10.0 -5.0
        for i,v in enumerate(np.reshape(vectorized(F_vec,t,m),vectorized(F,t,m).shape)):
            assert np.allclose(single(F[i],t,m),v)

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
