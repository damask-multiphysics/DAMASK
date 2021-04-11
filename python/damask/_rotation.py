import numpy as np

from . import tensor
from . import util
from . import grid_filters

_P = -1

# parameters for conversion from/to cubochoric
_sc   = np.pi**(1./6.)/6.**(1./6.)
_beta = np.pi**(5./6.)/6.**(1./6.)/2.
_R1   = (3.*np.pi/4.)**(1./3.)

class Rotation:
    u"""
    Rotation with functionality for conversion between different representations.

    The following conventions apply:

    - Coordinate frames are right-handed.
    - A rotation angle ω is taken to be positive for a counterclockwise rotation
      when viewing from the end point of the rotation axis towards the origin.
    - Rotations will be interpreted in the passive sense.
    - Euler angle triplets are implemented using the Bunge convention,
      with angular ranges of [0,2π], [0,π], [0,2π].
    - The rotation angle ω is limited to the interval [0,π].
    - The real part of a quaternion is positive, Re(q) > 0
    - P = -1 (as default).

    Examples
    --------
    Rotate vector 'a' (defined in coordinate system 'A') to
    coordinates 'b' expressed in system 'B':

    >>> import damask
    >>> import numpy as np
    >>> Q = damask.Rotation.from_random()
    >>> a = np.random.rand(3)
    >>> b = R @ a
    >>> np.allclose(np.dot(Q.as_matrix(),a),b)
    True

    Compound rotations R1 (first) and R2 (second):

    >>> import damask
    >>> import numpy as np
    >>> R1 = damask.Rotation.from_random()
    >>> R2 = damask.Rotation.from_random()
    >>> R = R2 * R1
    >>> np.allclose(R.as_matrix(), np.dot(R2.as_matrix(),R1.as_matrix()))
    True

    References
    ----------
    D. Rowenhorst et al., Modelling and Simulation in Materials Science and Engineering 23:083501, 2015
    https://doi.org/10.1088/0965-0393/23/8/083501

    """

    __slots__ = ['quaternion']

    def __init__(self,rotation = np.array([1.0,0.0,0.0,0.0])):
        """
        New rotation.

        Parameters
        ----------
        rotation : list, numpy.ndarray, Rotation, optional
            Unit quaternion in positive real hemisphere.
            Use .from_quaternion to perform a sanity check.
            Defaults to no rotation.

        """
        if isinstance(rotation,Rotation):
            self.quaternion = rotation.quaternion.copy()
        elif np.array(rotation).shape[-1] == 4:
            self.quaternion = np.array(rotation)
        else:
            raise TypeError('"rotation" is neither a Rotation nor a quaternion')


    def __repr__(self):
        """Represent rotation as unit quaternion(s)."""
        return f'Quaternion{" " if self.quaternion.shape == (4,) else "s of shape "+str(self.quaternion.shape[:-1])+chr(10)}'\
               + str(self.quaternion)


    def __copy__(self,**kwargs):
        """Create deep copy."""
        return self.__class__(rotation=kwargs['rotation'] if 'rotation' in kwargs else self.quaternion)

    copy = __copy__


    def __getitem__(self,item):
        """Return slice according to item."""
        return self.copy() \
               if self.shape == () else \
               self.copy(rotation=self.quaternion[item+(slice(None),)] if isinstance(item,tuple) else self.quaternion[item])


    def __eq__(self,other):
        """
        Equal to other.

        Parameters
        ----------
        other : Rotation
            Rotation to check for equality.

        """
        return np.logical_or(np.all(self.quaternion ==      other.quaternion,axis=-1),
                             np.all(self.quaternion == -1.0*other.quaternion,axis=-1))


    def __ne__(self,other):
        """
        Not equal to other.

        Parameters
        ----------
        other : Rotation
            Rotation to check for equality.

        """
        return np.logical_not(self==other)


    def isclose(self,other,rtol=1e-5,atol=1e-8,equal_nan=True):
        """
        Report where values are approximately equal to corresponding ones of other Rotation.

        Parameters
        ----------
        other : Rotation
            Rotation to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        mask : numpy.ndarray bool
            Mask indicating where corresponding rotations are close.

        """
        s = self.quaternion
        o = other.quaternion
        return np.logical_or(np.all(np.isclose(s,     o,rtol,atol,equal_nan),axis=-1),
                             np.all(np.isclose(s,-1.0*o,rtol,atol,equal_nan),axis=-1))


    def allclose(self,other,rtol=1e-5,atol=1e-8,equal_nan=True):
        """
        Test whether all values are approximately equal to corresponding ones of other Rotation.

        Parameters
        ----------
        other : Rotation
            Rotation to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        answer : bool
            Whether all values are close between both rotations.

        """
        return np.all(self.isclose(other,rtol,atol,equal_nan))


    def __array__(self):
        """Initializer for numpy."""
        return self.quaternion


    @property
    def size(self):
        return self.quaternion[...,0].size

    @property
    def shape(self):
        return self.quaternion[...,0].shape


    def __len__(self):
        """Length of leading/leftmost dimension of array."""
        return 0 if self.shape == () else self.shape[0]


    def __invert__(self):
        """Inverse rotation (backward rotation)."""
        dup = self.copy()
        dup.quaternion[...,1:] *= -1
        return dup


    def __pow__(self,exp):
        """
        Perform the rotation 'exp' times.

        Parameters
        ----------
        exp : float
            Exponent.

        """
        phi = np.arccos(self.quaternion[...,0:1])
        p = self.quaternion[...,1:]/np.linalg.norm(self.quaternion[...,1:],axis=-1,keepdims=True)
        return self.copy(rotation=Rotation(np.block([np.cos(exp*phi),np.sin(exp*phi)*p]))._standardize())

    def __ipow__(self,exp):
        """
        Perform the rotation 'exp' times (in-place).

        Parameters
        ----------
        exp : float
            Exponent.

        """
        return self**exp


    def __mul__(self,other):
        """
        Compose with other.

        Parameters
        ----------
        other : Rotation of shape(self.shape)
            Rotation for composition.

        Returns
        -------
        composition : Rotation
            Compound rotation self*other, i.e. first other then self rotation.

        """
        if isinstance(other,Rotation):
            q_m = self.quaternion[...,0:1]
            p_m = self.quaternion[...,1:]
            q_o = other.quaternion[...,0:1]
            p_o = other.quaternion[...,1:]
            q = (q_m*q_o - np.einsum('...i,...i',p_m,p_o).reshape(self.shape+(1,)))
            p = q_m*p_o + q_o*p_m + _P * np.cross(p_m,p_o)
            return Rotation(np.block([q,p]))._standardize()
        else:
            raise TypeError('Use "R@b", i.e. matmul, to apply rotation "R" to object "b"')

    def __imul__(self,other):
        """
        Compose with other (in-place).

        Parameters
        ----------
        other : Rotation of shape(self.shape)
            Rotation for composition.

        """
        return self*other


    def __truediv__(self,other):
        """
        Compose with inverse of other.

        Parameters
        ----------
        other : damask.Rotation of shape (self.shape)
            Rotation to inverse composition.

        Returns
        -------
        composition : Rotation
            Compound rotation self*(~other), i.e. first inverse of other then self rotation.

        """
        if isinstance(other,Rotation):
            return self*~other
        else:
            raise TypeError('Use "R@b", i.e. matmul, to apply rotation "R" to object "b"')

    def __itruediv__(self,other):
        """
        Compose with inverse of other (in-place).

        Parameters
        ----------
        other : Rotation of shape (self.shape)
            Rotation to inverse composition.

        """
        return self/other


    def __matmul__(self,other):
        """
        Rotate vector, second order tensor, or fourth order tensor.

        Parameters
        ----------
        other : numpy.ndarray of shape (...,3), (...,3,3), or (...,3,3,3,3)
            Vector or tensor on which to apply the rotation.

        Returns
        -------
        rotated : numpy.ndarray of shape (...,3), (...,3,3), or (...,3,3,3,3)
            Rotated vector or tensor, i.e. transformed to frame defined by rotation.

        """
        if isinstance(other,np.ndarray):
            if self.shape + (3,) == other.shape:
                q_m = self.quaternion[...,0]
                p_m = self.quaternion[...,1:]
                A = q_m**2.0 - np.einsum('...i,...i',p_m,p_m)
                B = 2.0 * np.einsum('...i,...i',p_m,other)
                C = 2.0 * _P * q_m
                return np.block([(A * other[...,i]).reshape(self.shape+(1,)) +
                                 (B * p_m[...,i]).reshape(self.shape+(1,)) +
                                 (C * (  p_m[...,(i+1)%3]*other[...,(i+2)%3]\
                                       - p_m[...,(i+2)%3]*other[...,(i+1)%3])).reshape(self.shape+(1,))
                                 for i in [0,1,2]])
            if self.shape + (3,3) == other.shape:
                R = self.as_matrix()
                return np.einsum('...im,...jn,...mn',R,R,other)
            if self.shape + (3,3,3,3) == other.shape:
                R = self.as_matrix()
                return np.einsum('...im,...jn,...ko,...lp,...mnop',R,R,R,R,other)
            else:
                raise ValueError('Can only rotate vectors, 2nd order tensors, and 4th order tensors')
        elif isinstance(other,Rotation):
            raise TypeError('Use "R1*R2", i.e. multiplication, to compose rotations "R1" and "R2"')
        else:
            raise TypeError(f'Cannot rotate {type(other)}')

    apply = __matmul__


    def _standardize(self):
        """Standardize quaternion (ensure positive real hemisphere)."""
        self.quaternion[self.quaternion[...,0] < 0.0] *= -1
        return self


    def append(self,other):
        """
        Extend array along first dimension with other array(s).

        Parameters
        ----------
            other : Rotation or list of Rotations.

        """
        return self.copy(rotation=np.vstack(tuple(map(lambda x:x.quaternion,
                                                      [self]+other if isinstance(other,list) else [self,other]))))


    def flatten(self,order = 'C'):
        """Flatten array."""
        return self.copy(rotation=self.quaternion.reshape((-1,4),order=order))


    def reshape(self,shape,order = 'C'):
        """Reshape array."""
        if isinstance(shape,(int,np.integer)): shape = (shape,)
        return self.copy(rotation=self.quaternion.reshape(tuple(shape)+(4,),order=order))


    def broadcast_to(self,shape,mode = 'right'):
        """
        Broadcast array.

        Parameters
        ----------
        shape : tuple
            Shape of broadcasted array.
        mode : str, optional
            Where to preferentially locate missing dimensions.
            Either 'left' or 'right' (default).

        """
        if isinstance(shape,(int,np.integer)): shape = (shape,)
        return self.copy(rotation=np.broadcast_to(self.quaternion.reshape(util.shapeshifter(self.shape,shape,mode)+(4,)),
                                                  shape+(4,)))


    def average(self,weights = None):
        """
        Average along last array dimension.

        Parameters
        ----------
        weights : list of floats, optional
            Relative weight of each rotation.

        Returns
        -------
        average : Rotation
            Weighted average of original Rotation field.

        References
        ----------
        F. Landis Markley et al., Journal of Guidance, Control, and Dynamics 30(4):1193-1197, 2007
        https://doi.org/10.2514/1.28949

        """
        def _M(quat):
            """Intermediate representation supporting quaternion averaging."""
            return np.einsum('...i,...j',quat,quat)

        if weights is None:
            weights = np.ones(self.shape,dtype=float)

        eig, vec = np.linalg.eig(np.sum(_M(self.quaternion) * weights[...,np.newaxis,np.newaxis],axis=-3) \
                                /np.sum(                      weights[...,np.newaxis,np.newaxis],axis=-3))

        return Rotation.from_quaternion(np.real(
                                        np.squeeze(
                                        np.take_along_axis(vec,
                                                           eig.argmax(axis=-1)[...,np.newaxis,np.newaxis],
                                                           axis=-1),
                                        axis=-1)),
                                        accept_homomorph = True)


    def misorientation(self,other):
        """
        Calculate misorientation to other Rotation.

        Parameters
        ----------
        other : Rotation
            Rotation to which the misorientation is computed.

        """
        return other*~self


    ################################################################################################
    # convert to different orientation representations (numpy arrays)

    def as_quaternion(self):
        """
        Represent as unit quaternion.

        Returns
        -------
        q : numpy.ndarray of shape (...,4)
            Unit quaternion in positive real hemisphere: (q_0, q_1, q_2, q_3), ǀqǀ=1, q_0 ≥ 0.

        """
        return self.quaternion.copy()

    def as_Euler_angles(self,
                        degrees = False):
        """
        Represent as Bunge-Euler angles.

        Parameters
        ----------
        degrees : bool, optional
            Return angles in degrees.

        Returns
        -------
        phi : numpy.ndarray of shape (...,3)
            Bunge-Euler angles: (φ_1, ϕ, φ_2), φ_1 ∈ [0,2π], ϕ ∈ [0,π], φ_2 ∈ [0,2π]
            unless degrees == True: φ_1 ∈ [0,360], ϕ ∈ [0,180], φ_2 ∈ [0,360]

        """
        eu = Rotation._qu2eu(self.quaternion)
        if degrees: eu = np.degrees(eu)
        return eu

    def as_axis_angle(self,
                      degrees = False,
                      pair = False):
        """
        Represent as axis angle pair.

        Parameters
        ----------
        degrees : bool, optional
            Return rotation angle in degrees. Defaults to False.
        pair : bool, optional
            Return tuple of axis and angle. Defaults to False.

        Returns
        -------
        axis_angle : numpy.ndarray of shape (...,4) unless pair == True:
            tuple containing numpy.ndarray of shapes (...,3) and (...)
            Axis angle pair: (n_1, n_2, n_3, ω), ǀnǀ = 1 and ω ∈ [0,π]
            unless degrees = True: ω ∈ [0,180].

        """
        ax = Rotation._qu2ax(self.quaternion)
        if degrees: ax[...,3] = np.degrees(ax[...,3])
        return (ax[...,:3],ax[...,3]) if pair else ax

    def as_matrix(self):
        """
        Represent as rotation matrix.

        Returns
        -------
        R : numpy.ndarray of shape (...,3,3)
            Rotation matrix R, det(R) = 1, R.T∙R=I.

        """
        return Rotation._qu2om(self.quaternion)

    def as_Rodrigues_vector(self,
                            compact = False):
        """
        Represent as Rodrigues-Frank vector with separated axis and angle argument.

        Parameters
        ----------
        compact : bool, optional
            Return as actual Rodrigues-Frank vector,
            i.e. axis and angle argument are not separated.

        Returns
        -------
        rho : numpy.ndarray of shape (...,4) containing 
            [n_1, n_2, n_3, tan(ω/2)], ǀnǀ = 1 and ω ∈ [0,π]
            unless compact == True:
            numpy.ndarray of shape (...,3) containing 
            tan(ω/2) [n_1, n_2, n_3], ω ∈ [0,π].

        """
        ro = Rotation._qu2ro(self.quaternion)
        if compact:
            with np.errstate(invalid='ignore'):
                return ro[...,:3]*ro[...,3:4]
        else:
            return ro

    def as_homochoric(self):
        """
        Represent as homochoric vector.

        Returns
        -------
        h : numpy.ndarray of shape (...,3)
            Homochoric vector: (h_1, h_2, h_3), ǀhǀ < (3/4*π)^(1/3).

        """
        return Rotation._qu2ho(self.quaternion)

    def as_cubochoric(self):
        """
        Represent as cubochoric vector.

        Returns
        -------
        x : numpy.ndarray of shape (...,3)
              Cubochoric vector: (x_1, x_2, x_3), max(x_i) < 1/2*π^(2/3).

        """
        return Rotation._qu2cu(self.quaternion)

    ################################################################################################
    # Static constructors. The input data needs to follow the conventions, options allow to
    # relax the conventions.
    @staticmethod
    def from_quaternion(q,
                        accept_homomorph = False,
                        P = -1):
        """
        Initialize from quaternion.

        Parameters
        ----------
        q : numpy.ndarray of shape (...,4)
            Unit quaternion in positive real hemisphere: (q_0, q_1, q_2, q_3),
            ǀqǀ=1, q_0 ≥ 0.
        accept_homomorph : boolean, optional
            Allow homomorphic variants, i.e. q_0 < 0 (negative real hemisphere).
            Defaults to False.
        P : int ∈ {-1,1}, optional
            Convention used. Defaults to -1.

        """
        qu = np.array(q,dtype=float)
        if qu.shape[:-2:-1] != (4,):
            raise ValueError('Invalid shape.')
        if abs(P) != 1:
            raise ValueError('P ∉ {-1,1}')

        qu[...,1:4] *= -P
        if accept_homomorph:
            qu[qu[...,0] < 0.0] *= -1
        else:
            if np.any(qu[...,0] < 0.0):
                raise ValueError('Quaternion with negative first (real) component.')
        if not np.all(np.isclose(np.linalg.norm(qu,axis=-1), 1.0)):
            raise ValueError('Quaternion is not of unit length.')

        return Rotation(qu)

    @staticmethod
    def from_Euler_angles(phi,
                          degrees = False):
        """
        Initialize from Bunge-Euler angles.

        Parameters
        ----------
        phi : numpy.ndarray of shape (...,3)
            Bunge-Euler angles: (φ_1, ϕ, φ_2), φ_1 ∈ [0,2π], ϕ ∈ [0,π], φ_2 ∈ [0,2π]
            unless degrees == True: φ_1 ∈ [0,360], ϕ ∈ [0,180], φ_2 ∈ [0,360].
        degrees : boolean, optional
            Bunge-Euler angles are given in degrees. Defaults to False.

        """
        eu = np.array(phi,dtype=float)
        if eu.shape[:-2:-1] != (3,):
            raise ValueError('Invalid shape.')

        eu = np.radians(eu) if degrees else eu
        if np.any(eu < 0.0) or np.any(eu > 2.0*np.pi) or np.any(eu[...,1] > np.pi):                 # ToDo: No separate check for PHI
            raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π].')

        return Rotation(Rotation._eu2qu(eu))

    @staticmethod
    def from_axis_angle(axis_angle,
                        degrees = False,
                        normalize = False,
                        P = -1):
        """
        Initialize from Axis angle pair.

        Parameters
        ----------
        axis_angle : numpy.ndarray of shape (...,4)
            Axis angle pair: [n_1, n_2, n_3, ω], ǀnǀ = 1 and ω ∈ [0,π]
            unless degrees = True: ω ∈ [0,180].
        degrees : boolean, optional
            Angle ω is given in degrees. Defaults to False.
        normalize: boolean, optional
            Allow ǀnǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Convention used. Defaults to -1.

        """
        ax = np.array(axis_angle,dtype=float)
        if ax.shape[:-2:-1] != (4,):
            raise ValueError('Invalid shape.')
        if abs(P) != 1:
            raise ValueError('P ∉ {-1,1}')

        ax[...,0:3] *= -P
        if degrees:   ax[...,  3]  = np.radians(ax[...,3])
        if normalize: ax[...,0:3] /= np.linalg.norm(ax[...,0:3],axis=-1,keepdims=True)
        if np.any(ax[...,3] < 0.0) or np.any(ax[...,3] > np.pi):
            raise ValueError('Axis angle rotation angle outside of [0..π].')
        if not np.all(np.isclose(np.linalg.norm(ax[...,0:3],axis=-1), 1.0)):
            print(np.linalg.norm(ax[...,0:3],axis=-1))
            raise ValueError('Axis angle rotation axis is not of unit length.')

        return Rotation(Rotation._ax2qu(ax))

    @staticmethod
    def from_basis(basis,
                   orthonormal = True,
                   reciprocal = False):
        """
        Initialize from lattice basis vectors.

        Parameters
        ----------
        basis : numpy.ndarray of shape (...,3,3)
            Three three-dimensional lattice basis vectors.
        orthonormal : boolean, optional
            Basis is strictly orthonormal, i.e. is free of stretch components. Defaults to True.
        reciprocal : boolean, optional
            Basis vectors are given in reciprocal (instead of real) space. Defaults to False.

        """
        om = np.array(basis,dtype=float)
        if om.shape[-2:] != (3,3):
            raise ValueError('Invalid shape.')

        if reciprocal:
            om = np.linalg.inv(tensor.transpose(om)/np.pi)                                          # transform reciprocal basis set
            orthonormal = False                                                                     # contains stretch
        if not orthonormal:
            (U,S,Vh) = np.linalg.svd(om)                                                            # singular value decomposition
            om = np.einsum('...ij,...jl',U,Vh)
        if not np.all(np.isclose(np.linalg.det(om),1.0)):
            raise ValueError('Orientation matrix has determinant ≠ 1.')
        if    not np.all(np.isclose(np.einsum('...i,...i',om[...,0],om[...,1]), 0.0)) \
           or not np.all(np.isclose(np.einsum('...i,...i',om[...,1],om[...,2]), 0.0)) \
           or not np.all(np.isclose(np.einsum('...i,...i',om[...,2],om[...,0]), 0.0)):
            raise ValueError('Orientation matrix is not orthogonal.')

        return Rotation(Rotation._om2qu(om))

    @staticmethod
    def from_matrix(R):
        """
        Initialize from rotation matrix.

        Parameters
        ----------
        R : numpy.ndarray of shape (...,3,3)
            Rotation matrix: det(R) = 1, R.T∙R=I.

        """
        return Rotation.from_basis(R)

    @staticmethod
    def from_parallel(a,b,
                      **kwargs):
        """
        Initialize from pairs of two orthogonal lattice basis vectors.

        Parameters
        ----------
        a : numpy.ndarray of shape (...,2,3)
            Two three-dimensional lattice vectors of first orthogonal basis.
        b : numpy.ndarray of shape (...,2,3)
            Corresponding three-dimensional lattice vectors of second basis.

        """
        a_ = np.array(a)
        b_ = np.array(b)
        if a_.shape[-2:] != (2,3) or b_.shape[-2:] != (2,3) or a_.shape != b_.shape:
            raise ValueError('Invalid shape.')
        am = np.stack([          a_[...,0,:],
                                             a_[...,1,:],
                        np.cross(a_[...,0,:],a_[...,1,:]) ],axis=-2)
        bm = np.stack([          b_[...,0,:],
                                             b_[...,1,:],
                        np.cross(b_[...,0,:],b_[...,1,:]) ],axis=-2)

        return              Rotation.from_basis(np.swapaxes(am/np.linalg.norm(am,axis=-1,keepdims=True),-1,-2))\
            .misorientation(Rotation.from_basis(np.swapaxes(bm/np.linalg.norm(bm,axis=-1,keepdims=True),-1,-2)))


    @staticmethod
    def from_Rodrigues_vector(rho,
                              normalize = False,
                              P = -1):
        """
        Initialize from Rodrigues-Frank vector (angle separated from axis).

        Parameters
        ----------
        rho : numpy.ndarray of shape (...,4)
            Rodrigues-Frank vector. (n_1, n_2, n_3, tan(ω/2)), ǀnǀ = 1  and ω ∈ [0,π].
        normalize : boolean, optional
            Allow ǀnǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Convention used. Defaults to -1.

        """
        ro = np.array(rho,dtype=float)
        if ro.shape[:-2:-1] != (4,):
            raise ValueError('Invalid shape.')
        if abs(P) != 1:
            raise ValueError('P ∉ {-1,1}')

        ro[...,0:3] *= -P
        if normalize: ro[...,0:3] /= np.linalg.norm(ro[...,0:3],axis=-1,keepdims=True)
        if np.any(ro[...,3] < 0.0):
            raise ValueError('Rodrigues vector rotation angle not positive.')
        if not np.all(np.isclose(np.linalg.norm(ro[...,0:3],axis=-1), 1.0)):
            raise ValueError('Rodrigues vector rotation axis is not of unit length.')

        return Rotation(Rotation._ro2qu(ro))

    @staticmethod
    def from_homochoric(h,
                        P = -1):
        """
        Initialize from homochoric vector.

        Parameters
        ----------
        h : numpy.ndarray of shape (...,3)
            Homochoric vector: (h_1, h_2, h_3), ǀhǀ < (3/4*π)^(1/3).
        P : int ∈ {-1,1}, optional
            Convention used. Defaults to -1.

        """
        ho = np.array(h,dtype=float)
        if ho.shape[:-2:-1] != (3,):
            raise ValueError('Invalid shape.')
        if abs(P) != 1:
            raise ValueError('P ∉ {-1,1}')

        ho *= -P

        if np.any(np.linalg.norm(ho,axis=-1) >_R1+1e-9):
            raise ValueError('Homochoric coordinate outside of the sphere.')

        return Rotation(Rotation._ho2qu(ho))

    @staticmethod
    def from_cubochoric(x,
                        P = -1):
        """
        Initialize from cubochoric vector.

        Parameters
        ----------
        x : numpy.ndarray of shape (...,3)
            Cubochoric vector: (x_1, x_2, x_3), max(x_i) < 1/2*π^(2/3).
        P : int ∈ {-1,1}, optional
            Convention used. Defaults to -1.

        """
        cu = np.array(x,dtype=float)
        if cu.shape[:-2:-1] != (3,):
            raise ValueError('Invalid shape.')
        if abs(P) != 1:
            raise ValueError('P ∉ {-1,1}')

        if np.abs(np.max(cu)) > np.pi**(2./3.) * 0.5+1e-9:
            raise ValueError('Cubochoric coordinate outside of the cube.')

        ho = -P * Rotation._cu2ho(cu)

        return Rotation(Rotation._ho2qu(ho))


    @staticmethod
    def from_random(shape = None,
                    rng_seed = None):
        """
        Draw random rotation.

        Rotations are uniformly distributed.

        Parameters
        ----------
        shape : tuple of ints, optional
            Shape of the sample. Defaults to None which gives a
            single rotation
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator. Defaults to None.
            If None, then fresh, unpredictable entropy will be pulled from the OS.

        """
        rng = np.random.default_rng(rng_seed)
        r = rng.random(3 if shape is None else tuple(shape)+(3,) if hasattr(shape, '__iter__') else (shape,3))

        A = np.sqrt(r[...,2])
        B = np.sqrt(1.0-r[...,2])
        q = np.stack([np.cos(2.0*np.pi*r[...,0])*A,
                      np.sin(2.0*np.pi*r[...,1])*B,
                      np.cos(2.0*np.pi*r[...,1])*B,
                      np.sin(2.0*np.pi*r[...,0])*A],axis=-1)

        return Rotation(q if shape is None else q.reshape(r.shape[:-1]+(4,)))._standardize()


    @staticmethod
    def from_ODF(weights,
                 phi,
                 N = 500,
                 degrees = True,
                 fractions = True,
                 rng_seed = None,
                 **kwargs):
        """
        Sample discrete values from a binned ODF.

        Parameters
        ----------
        weights : numpy.ndarray of shape (n)
            Texture intensity values (probability density or volume fraction) at Euler grid points.
        phi : numpy.ndarray of shape (n,3)
            Grid coordinates in Euler space at which weights are defined.
        N : integer, optional
            Number of discrete orientations to be sampled from the given ODF.
            Defaults to 500.
        degrees : boolean, optional
            Euler grid values are in degrees. Defaults to True.
        fractions : boolean, optional
            ODF values correspond to volume fractions, not probability density.
            Defaults to True.
        rng_seed: {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator. Defaults to None, i.e. unpredictable entropy
            will be pulled from the OS.

        Returns
        -------
        samples : damask.Rotation of shape (N)
            Array of sampled rotations closely representing the input ODF.

        Notes
        -----
        Due to the distortion of Euler space in the vicinity of ϕ = 0, probability densities, p, defined on
        grid points with ϕ = 0 will never result in reconstructed orientations as their dV/V = p dγ = p × 0.
        Hence, it is recommended to transform any such dataset to cell centers that avoid grid points at ϕ = 0.

        References
        ----------
        P. Eisenlohr and F. Roters, Computational Materials Science 42(4):670-678, 2008
        https://doi.org/10.1016/j.commatsci.2007.09.015

        """
        def _dg(eu,deg):
            """Return infinitesimal Euler space volume of bin(s)."""
            phi_sorted = eu[np.lexsort((eu[:,0],eu[:,1],eu[:,2]))]
            steps,size,_ = grid_filters.cellsSizeOrigin_coordinates0_point(phi_sorted)
            delta = np.radians(size/steps) if deg else size/steps
            return delta[0]*2.0*np.sin(delta[1]/2.0)*delta[2] / 8.0 / np.pi**2 * np.sin(np.radians(eu[:,1]) if deg else eu[:,1])

        dg = 1.0 if fractions else _dg(phi,degrees)
        dV_V = dg * np.maximum(0.0,weights.squeeze())

        return Rotation.from_Euler_angles(phi[util.hybrid_IA(dV_V,N,rng_seed)],degrees)


    @staticmethod
    def from_spherical_component(center,
                                 sigma,
                                 N = 500,
                                 degrees = True,
                                 rng_seed = None):
        """
        Calculate set of rotations with Gaussian distribution around center.

        Parameters
        ----------
        center : Rotation
            Central Rotation.
        sigma : float
            Standard deviation of (Gaussian) misorientation distribution.
        N : int, optional
            Number of samples, defaults to 500.
        degrees : boolean, optional
            sigma is given in degrees.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator. Defaults to None, i.e. unpredictable entropy
            will be pulled from the OS.

        """
        rng = np.random.default_rng(rng_seed)
        sigma = np.radians(sigma) if degrees else sigma
        u,Theta  = (rng.random((N,2)) * 2.0 * np.array([1,np.pi]) - np.array([1.0, 0])).T
        omega = abs(rng.normal(scale=sigma,size=N))
        p = np.column_stack([np.sqrt(1-u**2)*np.cos(Theta),
                             np.sqrt(1-u**2)*np.sin(Theta),
                             u, omega])

        return  Rotation.from_axis_angle(p) * center


    @staticmethod
    def from_fiber_component(alpha,
                             beta,
                             sigma = 0.0,
                             N = 500,
                             degrees = True,
                             rng_seed = None):
        """
        Calculate set of rotations with Gaussian distribution around direction.

        Parameters
        ----------
        alpha : numpy.ndarray of size 2
            Polar coordinates (phi from x,theta from z) of fiber direction in crystal frame.
        beta : numpy.ndarray of size 2
            Polar coordinates (phi from x,theta from z) of fiber direction in sample frame.
        sigma : float, optional
            Standard deviation of (Gaussian) misorientation distribution.
            Defaults to 0.
        N : int, optional
            Number of samples, defaults to 500.
        degrees : boolean, optional
            sigma, alpha, and beta are given in degrees.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator. Defaults to None, i.e. unpredictable entropy
            will be pulled from the OS.

        """
        rng = np.random.default_rng(rng_seed)
        sigma_,alpha_,beta_ = map(np.radians,(sigma,alpha,beta)) if degrees else (sigma,alpha,beta)

        d_cr  = np.array([np.sin(alpha_[0])*np.cos(alpha_[1]), np.sin(alpha_[0])*np.sin(alpha_[1]), np.cos(alpha_[0])])
        d_lab = np.array([np.sin( beta_[0])*np.cos( beta_[1]), np.sin( beta_[0])*np.sin( beta_[1]), np.cos( beta_[0])])
        ax_align = np.append(np.cross(d_lab,d_cr), np.arccos(np.dot(d_lab,d_cr)))
        if np.isclose(ax_align[3],0.0): ax_align[:3] = np.array([1,0,0])
        R_align  = Rotation.from_axis_angle(ax_align if ax_align[3] > 0.0 else -ax_align,normalize=True) # rotate fiber axis from sample to crystal frame

        u,Theta  = (rng.random((N,2)) * 2.0 * np.array([1,np.pi]) - np.array([1.0, 0])).T
        omega  = abs(rng.normal(scale=sigma_,size=N))
        p = np.column_stack([np.sqrt(1-u**2)*np.cos(Theta),
                             np.sqrt(1-u**2)*np.sin(Theta),
                             u, omega])
        p[:,:3] = np.einsum('ij,...j',np.eye(3)-np.outer(d_lab,d_lab),p[:,:3])                      # remove component along fiber axis
        f = np.column_stack((np.broadcast_to(d_lab,(N,3)),rng.random(N)*np.pi))
        f[::2,:3] *= -1                                                                             # flip half the rotation axes to negative sense

        return R_align.broadcast_to(N) \
             * Rotation.from_axis_angle(p,normalize=True) \
             * Rotation.from_axis_angle(f)


####################################################################################################
# Code below available according to the following conditions on https://github.com/MarDiehl/3Drotations
####################################################################################################
# Copyright (c) 2017-2020, Martin Diehl/Max-Planck-Institut für Eisenforschung GmbH
# Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
#     - Redistributions of source code must retain the above copyright notice, this list
#        of conditions and the following disclaimer.
#     - Redistributions in binary form must reproduce the above copyright notice, this
#        list of conditions and the following disclaimer in the documentation and/or
#        other materials provided with the distribution.
#     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names
#        of its contributors may be used to endorse or promote products derived from
#        this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
####################################################################################################
    #---------- Quaternion ----------
    @staticmethod
    def _qu2om(qu):
        qq = qu[...,0:1]**2-(qu[...,1:2]**2 + qu[...,2:3]**2 + qu[...,3:4]**2)
        om = np.block([qq + 2.0*qu[...,1:2]**2,
                       2.0*(qu[...,2:3]*qu[...,1:2]-_P*qu[...,0:1]*qu[...,3:4]),
                       2.0*(qu[...,3:4]*qu[...,1:2]+_P*qu[...,0:1]*qu[...,2:3]),
                       2.0*(qu[...,1:2]*qu[...,2:3]+_P*qu[...,0:1]*qu[...,3:4]),
                       qq + 2.0*qu[...,2:3]**2,
                       2.0*(qu[...,3:4]*qu[...,2:3]-_P*qu[...,0:1]*qu[...,1:2]),
                       2.0*(qu[...,1:2]*qu[...,3:4]-_P*qu[...,0:1]*qu[...,2:3]),
                       2.0*(qu[...,2:3]*qu[...,3:4]+_P*qu[...,0:1]*qu[...,1:2]),
                       qq + 2.0*qu[...,3:4]**2,
                      ]).reshape(qu.shape[:-1]+(3,3))
        return om

    @staticmethod
    def _qu2eu(qu):
        """Quaternion to Bunge-Euler angles."""
        q02   = qu[...,0:1]*qu[...,2:3]
        q13   = qu[...,1:2]*qu[...,3:4]
        q01   = qu[...,0:1]*qu[...,1:2]
        q23   = qu[...,2:3]*qu[...,3:4]
        q03_s = qu[...,0:1]**2+qu[...,3:4]**2
        q12_s = qu[...,1:2]**2+qu[...,2:3]**2
        chi = np.sqrt(q03_s*q12_s)

        eu = np.where(np.abs(q12_s) < 1.0e-8,
                np.block([np.arctan2(-_P*2.0*qu[...,0:1]*qu[...,3:4],qu[...,0:1]**2-qu[...,3:4]**2),
                          np.zeros(qu.shape[:-1]+(2,))]),
                      np.where(np.abs(q03_s) < 1.0e-8,
                          np.block([np.arctan2(   2.0*qu[...,1:2]*qu[...,2:3],qu[...,1:2]**2-qu[...,2:3]**2),
                                    np.broadcast_to(np.pi,qu[...,0:1].shape),
                                    np.zeros(qu.shape[:-1]+(1,))]),
                          np.block([np.arctan2((-_P*q02+q13)*chi, (-_P*q01-q23)*chi),
                                    np.arctan2( 2.0*chi,          q03_s-q12_s    ),
                                    np.arctan2(( _P*q02+q13)*chi, (-_P*q01+q23)*chi)])
                              )
                     )
        # reduce Euler angles to definition range
        eu[np.abs(eu)<1.e-6] = 0.0
        eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)                # needed?
        return eu

    @staticmethod
    def _qu2ax(qu):
        """
        Quaternion to axis angle pair.

        Modified version of the original formulation, should be numerically more stable
        """
        with np.errstate(invalid='ignore',divide='ignore'):
            s = np.sign(qu[...,0:1])/np.sqrt(qu[...,1:2]**2+qu[...,2:3]**2+qu[...,3:4]**2)
            omega = 2.0 * np.arccos(np.clip(qu[...,0:1],-1.0,1.0))
            ax = np.where(np.broadcast_to(qu[...,0:1] < 1.0e-8,qu.shape),
                          np.block([qu[...,1:4],np.broadcast_to(np.pi,qu[...,0:1].shape)]),
                          np.block([qu[...,1:4]*s,omega]))
        ax[np.isclose(qu[...,0],1.,rtol=0.0)] = [0.0, 0.0, 1.0, 0.0]
        return ax

    @staticmethod
    def _qu2ro(qu):
        """Quaternion to Rodrigues-Frank vector."""
        with np.errstate(invalid='ignore',divide='ignore'):
            s  = np.linalg.norm(qu[...,1:4],axis=-1,keepdims=True)
            ro = np.where(np.broadcast_to(np.abs(qu[...,0:1]) < 1.0e-12,qu.shape),
                          np.block([qu[...,1:2], qu[...,2:3], qu[...,3:4], np.broadcast_to(np.inf,qu[...,0:1].shape)]),
                          np.block([qu[...,1:2]/s,qu[...,2:3]/s,qu[...,3:4]/s,
                                    np.tan(np.arccos(np.clip(qu[...,0:1],-1.0,1.0)))
                                   ])
                       )
        ro[np.abs(s).squeeze(-1) < 1.0e-12] = [0.0,0.0,_P,0.0]
        return ro

    @staticmethod
    def _qu2ho(qu):
        """Quaternion to homochoric vector."""
        with np.errstate(invalid='ignore'):
            omega = 2.0 * np.arccos(np.clip(qu[...,0:1],-1.0,1.0))
            ho = np.where(np.abs(omega) < 1.0e-12,
                          np.zeros(3),
                          qu[...,1:4]/np.linalg.norm(qu[...,1:4],axis=-1,keepdims=True) \
                          * (0.75*(omega - np.sin(omega)))**(1./3.))
        return ho

    @staticmethod
    def _qu2cu(qu):
        """Quaternion to cubochoric vector."""
        return Rotation._ho2cu(Rotation._qu2ho(qu))


    #---------- Rotation matrix ----------
    @staticmethod
    def _om2qu(om):
        """
        Rotation matrix to quaternion.

        This formulation is from  www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion.
        The original formulation had issues.
        """
        trace = om[...,0,0:1]+om[...,1,1:2]+om[...,2,2:3]

        with np.errstate(invalid='ignore',divide='ignore'):
            s = [
                 0.5 / np.sqrt( 1.0 + trace),
                 2.0 * np.sqrt( 1.0 + om[...,0,0:1] - om[...,1,1:2] - om[...,2,2:3]),
                 2.0 * np.sqrt( 1.0 + om[...,1,1:2] - om[...,2,2:3] - om[...,0,0:1]),
                 2.0 * np.sqrt( 1.0 + om[...,2,2:3] - om[...,0,0:1] - om[...,1,1:2] )
                 ]
            qu= np.where(trace>0,
                         np.block([0.25 / s[0],
                                  (om[...,2,1:2] - om[...,1,2:3] ) * s[0],
                                  (om[...,0,2:3] - om[...,2,0:1] ) * s[0],
                                  (om[...,1,0:1] - om[...,0,1:2] ) * s[0]]),
                         np.where(om[...,0,0:1] > np.maximum(om[...,1,1:2],om[...,2,2:3]),
                                  np.block([(om[...,2,1:2] - om[...,1,2:3]) / s[1],
                                            0.25 * s[1],
                                            (om[...,0,1:2] + om[...,1,0:1]) / s[1],
                                            (om[...,0,2:3] + om[...,2,0:1]) / s[1]]),
                                  np.where(om[...,1,1:2] > om[...,2,2:3],
                                           np.block([(om[...,0,2:3] - om[...,2,0:1]) / s[2],
                                                     (om[...,0,1:2] + om[...,1,0:1]) / s[2],
                                                     0.25 * s[2],
                                                     (om[...,1,2:3] + om[...,2,1:2]) / s[2]]),
                                           np.block([(om[...,1,0:1] - om[...,0,1:2]) / s[3],
                                                     (om[...,0,2:3] + om[...,2,0:1]) / s[3],
                                                     (om[...,1,2:3] + om[...,2,1:2]) / s[3],
                                                     0.25 * s[3]]),
                                          )
                                 )
                        )*np.array([1,_P,_P,_P])
            qu[qu[...,0]<0] *=-1
        return qu

    @staticmethod
    def _om2eu(om):
        """Rotation matrix to Bunge-Euler angles."""
        with np.errstate(invalid='ignore',divide='ignore'):
            zeta = 1.0/np.sqrt(1.0-om[...,2,2:3]**2)
            eu = np.where(np.isclose(np.abs(om[...,2,2:3]),1.0,0.0),
                          np.block([np.arctan2(om[...,0,1:2],om[...,0,0:1]),
                                    np.pi*0.5*(1-om[...,2,2:3]),
                                    np.zeros(om.shape[:-2]+(1,)),
                                   ]),
                          np.block([np.arctan2(om[...,2,0:1]*zeta,-om[...,2,1:2]*zeta),
                                    np.arccos( om[...,2,2:3]),
                                    np.arctan2(om[...,0,2:3]*zeta,+om[...,1,2:3]*zeta)
                                   ])
                          )
        eu[np.abs(eu)<1.e-8] = 0.0
        eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
        return eu

    @staticmethod
    def _om2ax(om):
        """Rotation matrix to axis angle pair."""
        diag_delta = -_P*np.block([om[...,1,2:3]-om[...,2,1:2],
                                   om[...,2,0:1]-om[...,0,2:3],
                                   om[...,0,1:2]-om[...,1,0:1]
                                 ])
        t = 0.5*(om.trace(axis2=-2,axis1=-1) -1.0).reshape(om.shape[:-2]+(1,))
        w,vr = np.linalg.eig(om)
        # mask duplicated real eigenvalues
        w[np.isclose(w[...,0],1.0+0.0j),1:] = 0.
        w[np.isclose(w[...,1],1.0+0.0j),2:] = 0.
        vr = np.swapaxes(vr,-1,-2)
        ax = np.where(np.abs(diag_delta)<1e-12,
                             np.real(vr[np.isclose(w,1.0+0.0j)]).reshape(om.shape[:-2]+(3,)),
                      np.abs(np.real(vr[np.isclose(w,1.0+0.0j)]).reshape(om.shape[:-2]+(3,))) \
                      *np.sign(diag_delta))
        ax = np.block([ax,np.arccos(np.clip(t,-1.0,1.0))])
        ax[np.abs(ax[...,3])<1.e-8] = [ 0.0, 0.0, 1.0, 0.0]
        return ax

    @staticmethod
    def _om2ro(om):
        """Rotation matrix to Rodrigues-Frank vector."""
        return Rotation._eu2ro(Rotation._om2eu(om))

    @staticmethod
    def _om2ho(om):
        """Rotation matrix to homochoric vector."""
        return Rotation._ax2ho(Rotation._om2ax(om))

    @staticmethod
    def _om2cu(om):
        """Rotation matrix to cubochoric vector."""
        return Rotation._ho2cu(Rotation._om2ho(om))


    #---------- Bunge-Euler angles ----------
    @staticmethod
    def _eu2qu(eu):
        """Bunge-Euler angles to quaternion."""
        ee = 0.5*eu
        cPhi = np.cos(ee[...,1:2])
        sPhi = np.sin(ee[...,1:2])
        qu = np.block([    cPhi*np.cos(ee[...,0:1]+ee[...,2:3]),
                       -_P*sPhi*np.cos(ee[...,0:1]-ee[...,2:3]),
                       -_P*sPhi*np.sin(ee[...,0:1]-ee[...,2:3]),
                       -_P*cPhi*np.sin(ee[...,0:1]+ee[...,2:3])])
        qu[qu[...,0]<0.0]*=-1
        return qu

    @staticmethod
    def _eu2om(eu):
        """Bunge-Euler angles to rotation matrix."""
        c = np.cos(eu)
        s = np.sin(eu)
        om = np.block([+c[...,0:1]*c[...,2:3]-s[...,0:1]*s[...,2:3]*c[...,1:2],
                       +s[...,0:1]*c[...,2:3]+c[...,0:1]*s[...,2:3]*c[...,1:2],
                       +s[...,2:3]*s[...,1:2],
                       -c[...,0:1]*s[...,2:3]-s[...,0:1]*c[...,2:3]*c[...,1:2],
                       -s[...,0:1]*s[...,2:3]+c[...,0:1]*c[...,2:3]*c[...,1:2],
                       +c[...,2:3]*s[...,1:2],
                       +s[...,0:1]*s[...,1:2],
                       -c[...,0:1]*s[...,1:2],
                       +c[...,1:2]
                       ]).reshape(eu.shape[:-1]+(3,3))
        om[np.abs(om)<1.e-12] = 0.0
        return om

    @staticmethod
    def _eu2ax(eu):
        """Bunge-Euler angles to axis angle pair."""
        t = np.tan(eu[...,1:2]*0.5)
        sigma = 0.5*(eu[...,0:1]+eu[...,2:3])
        delta = 0.5*(eu[...,0:1]-eu[...,2:3])
        tau   = np.linalg.norm(np.block([t,np.sin(sigma)]),axis=-1,keepdims=True)
        alpha = np.where(np.abs(np.cos(sigma))<1.e-12,np.pi,2.0*np.arctan(tau/np.cos(sigma)))
        with np.errstate(invalid='ignore',divide='ignore'):
            ax = np.where(np.broadcast_to(np.abs(alpha)<1.0e-12,eu.shape[:-1]+(4,)),
                          [0.0,0.0,1.0,0.0],
                          np.block([-_P/tau*t*np.cos(delta),
                                    -_P/tau*t*np.sin(delta),
                                    -_P/tau*  np.sin(sigma),
                                     alpha
                                    ]))
        ax[(alpha<0.0).squeeze()] *=-1
        return ax

    @staticmethod
    def _eu2ro(eu):
        """Bunge-Euler angles to Rodrigues-Frank vector."""
        ax = Rotation._eu2ax(eu)
        ro = np.block([ax[...,:3],np.tan(ax[...,3:4]*.5)])
        ro[ax[...,3]>=np.pi,3] = np.inf
        ro[np.abs(ax[...,3])<1.e-16] = [ 0.0, 0.0, _P, 0.0 ]
        return ro

    @staticmethod
    def _eu2ho(eu):
        """Bunge-Euler angles to homochoric vector."""
        return Rotation._ax2ho(Rotation._eu2ax(eu))

    @staticmethod
    def _eu2cu(eu):
        """Bunge-Euler angles to cubochoric vector."""
        return Rotation._ho2cu(Rotation._eu2ho(eu))


    #---------- Axis angle pair ----------
    @staticmethod
    def _ax2qu(ax):
        """Axis angle pair to quaternion."""
        c = np.cos(ax[...,3:4]*.5)
        s = np.sin(ax[...,3:4]*.5)
        qu = np.where(np.abs(ax[...,3:4])<1.e-6,[1.0, 0.0, 0.0, 0.0],np.block([c, ax[...,:3]*s]))
        return qu

    @staticmethod
    def _ax2om(ax):
        """Axis angle pair to rotation matrix."""
        c = np.cos(ax[...,3:4])
        s = np.sin(ax[...,3:4])
        omc = 1. -c
        om = np.block([c+omc*ax[...,0:1]**2,
                       omc*ax[...,0:1]*ax[...,1:2] + s*ax[...,2:3],
                       omc*ax[...,0:1]*ax[...,2:3] - s*ax[...,1:2],
                       omc*ax[...,0:1]*ax[...,1:2] - s*ax[...,2:3],
                       c+omc*ax[...,1:2]**2,
                       omc*ax[...,1:2]*ax[...,2:3] + s*ax[...,0:1],
                       omc*ax[...,0:1]*ax[...,2:3] + s*ax[...,1:2],
                       omc*ax[...,1:2]*ax[...,2:3] - s*ax[...,0:1],
                       c+omc*ax[...,2:3]**2]).reshape(ax.shape[:-1]+(3,3))
        return om if _P < 0.0 else np.swapaxes(om,-1,-2)

    @staticmethod
    def _ax2eu(ax):
        """Rotation matrix to Bunge Euler angles."""
        return Rotation._om2eu(Rotation._ax2om(ax))

    @staticmethod
    def _ax2ro(ax):
        """Axis angle pair to Rodrigues-Frank vector."""
        ro = np.block([ax[...,:3],
                       np.where(np.isclose(ax[...,3:4],np.pi,atol=1.e-15,rtol=.0),
                                np.inf,
                                np.tan(ax[...,3:4]*0.5))
                      ])
        ro[np.abs(ax[...,3])<1.e-6] = [.0,.0,_P,.0]
        return ro

    @staticmethod
    def _ax2ho(ax):
        """Axis angle pair to homochoric vector."""
        f = (0.75 * ( ax[...,3:4] - np.sin(ax[...,3:4]) ))**(1.0/3.0)
        ho = ax[...,:3] * f
        return ho

    @staticmethod
    def _ax2cu(ax):
        """Axis angle pair to cubochoric vector."""
        return Rotation._ho2cu(Rotation._ax2ho(ax))


    #---------- Rodrigues-Frank vector ----------
    @staticmethod
    def _ro2qu(ro):
        """Rodrigues-Frank vector to quaternion."""
        return Rotation._ax2qu(Rotation._ro2ax(ro))

    @staticmethod
    def _ro2om(ro):
        """Rodgrigues-Frank vector to rotation matrix."""
        return Rotation._ax2om(Rotation._ro2ax(ro))

    @staticmethod
    def _ro2eu(ro):
        """Rodrigues-Frank vector to Bunge-Euler angles."""
        return Rotation._om2eu(Rotation._ro2om(ro))

    @staticmethod
    def _ro2ax(ro):
        """Rodrigues-Frank vector to axis angle pair."""
        with np.errstate(invalid='ignore',divide='ignore'):
            ax = np.where(np.isfinite(ro[...,3:4]),
                 np.block([ro[...,0:3]*np.linalg.norm(ro[...,0:3],axis=-1,keepdims=True),2.*np.arctan(ro[...,3:4])]),
                 np.block([ro[...,0:3],np.broadcast_to(np.pi,ro[...,3:4].shape)]))
        ax[np.abs(ro[...,3]) < 1.e-8]  = np.array([ 0.0, 0.0, 1.0, 0.0 ])
        return ax

    @staticmethod
    def _ro2ho(ro):
        """Rodrigues-Frank vector to homochoric vector."""
        f = np.where(np.isfinite(ro[...,3:4]),2.0*np.arctan(ro[...,3:4]) -np.sin(2.0*np.arctan(ro[...,3:4])),np.pi)
        ho = np.where(np.broadcast_to(np.sum(ro[...,0:3]**2.0,axis=-1,keepdims=True) < 1.e-8,ro[...,0:3].shape),
                      np.zeros(3), ro[...,0:3]* (0.75*f)**(1.0/3.0))
        return ho

    @staticmethod
    def _ro2cu(ro):
        """Rodrigues-Frank vector to cubochoric vector."""
        return Rotation._ho2cu(Rotation._ro2ho(ro))


    #---------- Homochoric vector----------
    @staticmethod
    def _ho2qu(ho):
        """Homochoric vector to quaternion."""
        return Rotation._ax2qu(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2om(ho):
        """Homochoric vector to rotation matrix."""
        return Rotation._ax2om(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2eu(ho):
        """Homochoric vector to Bunge-Euler angles."""
        return Rotation._ax2eu(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2ax(ho):
        """Homochoric vector to axis angle pair."""
        tfit = np.array([+1.0000000000018852,      -0.5000000002194847,
                         -0.024999992127593126,    -0.003928701544781374,
                         -0.0008152701535450438,   -0.0002009500426119712,
                         -0.00002397986776071756,  -0.00008202868926605841,
                         +0.00012448715042090092,  -0.0001749114214822577,
                         +0.0001703481934140054,   -0.00012062065004116828,
                         +0.000059719705868660826, -0.00001980756723965647,
                         +0.000003953714684212874, -0.00000036555001439719544])
        hmag_squared = np.sum(ho**2.,axis=-1,keepdims=True)
        hm = hmag_squared.copy()
        s = tfit[0] + tfit[1] * hmag_squared
        for i in range(2,16):
            hm *= hmag_squared
            s  += tfit[i] * hm
        with np.errstate(invalid='ignore'):
            ax = np.where(np.broadcast_to(np.abs(hmag_squared)<1.e-8,ho.shape[:-1]+(4,)),
                          [ 0.0, 0.0, 1.0, 0.0 ],
                          np.block([ho/np.sqrt(hmag_squared),2.0*np.arccos(np.clip(s,-1.0,1.0))]))
        return ax

    @staticmethod
    def _ho2ro(ho):
        """Axis angle pair to Rodrigues-Frank vector."""
        return Rotation._ax2ro(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2cu(ho):
        """
        Homochoric vector to cubochoric vector.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        rs = np.linalg.norm(ho,axis=-1,keepdims=True)

        xyz3 = np.take_along_axis(ho,Rotation._get_pyramid_order(ho,'forward'),-1)

        with np.errstate(invalid='ignore',divide='ignore'):
            # inverse M_3
            xyz2 = xyz3[...,0:2] * np.sqrt( 2.0*rs/(rs+np.abs(xyz3[...,2:3])) )
            qxy = np.sum(xyz2**2,axis=-1,keepdims=True)

            q2 = qxy + np.max(np.abs(xyz2),axis=-1,keepdims=True)**2
            sq2 = np.sqrt(q2)
            q = (_beta/np.sqrt(2.0)/_R1) * np.sqrt(q2*qxy/(q2-np.max(np.abs(xyz2),axis=-1,keepdims=True)*sq2))
            tt = np.clip((np.min(np.abs(xyz2),axis=-1,keepdims=True)**2\
                +np.max(np.abs(xyz2),axis=-1,keepdims=True)*sq2)/np.sqrt(2.0)/qxy,-1.0,1.0)
            T_inv = np.where(np.abs(xyz2[...,1:2]) <= np.abs(xyz2[...,0:1]),
                                np.block([np.ones_like(tt),np.arccos(tt)/np.pi*12.0]),
                                np.block([np.arccos(tt)/np.pi*12.0,np.ones_like(tt)]))*q
            T_inv[xyz2<0.0] *= -1.0
            T_inv[np.broadcast_to(np.isclose(qxy,0.0,rtol=0.0,atol=1.0e-12),T_inv.shape)] = 0.0
            cu = np.block([T_inv, np.where(xyz3[...,2:3]<0.0,-np.ones_like(xyz3[...,2:3]),np.ones_like(xyz3[...,2:3])) \
                                  * rs/np.sqrt(6.0/np.pi),
                          ])/ _sc

        cu[np.isclose(np.sum(np.abs(ho),axis=-1),0.0,rtol=0.0,atol=1.0e-16)] = 0.0
        cu = np.take_along_axis(cu,Rotation._get_pyramid_order(ho,'backward'),-1)

        return cu

    #---------- Cubochoric ----------
    @staticmethod
    def _cu2qu(cu):
        """Cubochoric vector to quaternion."""
        return Rotation._ho2qu(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2om(cu):
        """Cubochoric vector to rotation matrix."""
        return Rotation._ho2om(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2eu(cu):
        """Cubochoric vector to Bunge-Euler angles."""
        return Rotation._ho2eu(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2ax(cu):
        """Cubochoric vector to axis angle pair."""
        return Rotation._ho2ax(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2ro(cu):
        """Cubochoric vector to Rodrigues-Frank vector."""
        return Rotation._ho2ro(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2ho(cu):
        """
        Cubochoric vector to homochoric vector.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        with np.errstate(invalid='ignore',divide='ignore'):
            # get pyramide and scale by grid parameter ratio
            XYZ = np.take_along_axis(cu,Rotation._get_pyramid_order(cu,'forward'),-1) * _sc
            order = np.abs(XYZ[...,1:2]) <= np.abs(XYZ[...,0:1])
            q = np.pi/12.0 * np.where(order,XYZ[...,1:2],XYZ[...,0:1]) \
                           / np.where(order,XYZ[...,0:1],XYZ[...,1:2])
            c = np.cos(q)
            s = np.sin(q)
            q = _R1*2.0**0.25/_beta/ np.sqrt(np.sqrt(2.0)-c) \
              * np.where(order,XYZ[...,0:1],XYZ[...,1:2])

            T = np.block([ (np.sqrt(2.0)*c - 1.0), np.sqrt(2.0) * s]) * q

            # transform to sphere grid (inverse Lambert)
            c = np.sum(T**2,axis=-1,keepdims=True)
            s = c *         np.pi/24.0 /XYZ[...,2:3]**2
            c = c * np.sqrt(np.pi/24.0)/XYZ[...,2:3]
            q = np.sqrt( 1.0 - s)

            ho = np.where(np.isclose(np.sum(np.abs(XYZ[...,0:2]),axis=-1,keepdims=True),0.0,rtol=0.0,atol=1.0e-16),
                          np.block([np.zeros_like(XYZ[...,0:2]),np.sqrt(6.0/np.pi) *XYZ[...,2:3]]),
                          np.block([np.where(order,T[...,0:1],T[...,1:2])*q,
                                    np.where(order,T[...,1:2],T[...,0:1])*q,
                                    np.sqrt(6.0/np.pi) * XYZ[...,2:3] - c])
                          )

        ho[np.isclose(np.sum(np.abs(cu),axis=-1),0.0,rtol=0.0,atol=1.0e-16)] = 0.0
        ho = np.take_along_axis(ho,Rotation._get_pyramid_order(cu,'backward'),-1)

        return ho


    @staticmethod
    def _get_pyramid_order(xyz,direction=None):
        """
        Get order of the coordinates.

        Depending on the pyramid in which the point is located, the order need to be adjusted.

        Parameters
        ----------
        xyz : numpy.ndarray
           coordinates of a point on a uniform refinable grid on a ball or
           in a uniform refinable cubical grid.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        order = {'forward': np.array([[0,1,2],[1,2,0],[2,0,1]]),
                 'backward':np.array([[0,1,2],[2,0,1],[1,2,0]])}

        p = np.where(np.maximum(np.abs(xyz[...,0]),np.abs(xyz[...,1])) <= np.abs(xyz[...,2]),0,
                     np.where(np.maximum(np.abs(xyz[...,1]),np.abs(xyz[...,2])) <= np.abs(xyz[...,0]),1,2))

        return order[direction][p]
