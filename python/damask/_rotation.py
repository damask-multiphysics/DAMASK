import copy
from typing import Optional, Union, Sequence, Tuple, Literal, List, TypeVar

import numpy as np

from ._typehints import FloatSequence, IntSequence, NumpyRngSeed
from . import tensor
from . import util
from . import grid_filters

_P = -1

# parameters for conversion from/to cubochoric
_sc   = np.pi**(1./6.)/6.**(1./6.)
_beta = np.pi**(5./6.)/6.**(1./6.)/2.
_R1   = (3.*np.pi/4.)**(1./3.)

MyType = TypeVar('MyType', bound='Rotation')

class Rotation:
    u"""
    Rotation with functionality for conversion between different representations.

    The following conventions apply:

    - Coordinate frames are right-handed.
    - A rotation angle ω is taken to be positive for a counterclockwise rotation
      when viewing from the end point of the rotation axis towards the origin.
    - Rotations will be interpreted in the passive sense, i.e. as rotation of
      the coordinate frame.
    - P = -1 (as default).

    Examples
    --------
    Rotate vector 'a' (defined in coordinate system 'A') to
    coordinates 'b' expressed in system 'B':

    >>> import numpy as np
    >>> import damask
    >>> Q = damask.Rotation.from_random()
    >>> a = np.random.rand(3)
    >>> b = Q @ a
    >>> np.allclose(np.dot(Q.as_matrix(),a),b)
    True

    Compound rotations R1 (first) and R2 (second):

    >>> import numpy as np
    >>> import damask
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

    def __init__(self,
                 rotation: Union[FloatSequence, 'Rotation'] = np.array([1.,0.,0.,0.])):
        """
        New rotation.

        Parameters
        ----------
        rotation : list, numpy.ndarray, or Rotation, optional
            Unit quaternion in positive real hemisphere.
            Use .from_quaternion to perform a sanity check.
            Defaults to no rotation.

        """
        self.quaternion: np.ndarray
        if isinstance(rotation,Rotation):
            self.quaternion = rotation.quaternion.copy()
        elif np.array(rotation).shape[-1] == 4:
            self.quaternion = np.array(rotation,dtype=float)
        else:
            raise TypeError('"rotation" is neither a Rotation nor a quaternion')


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.

        """
        return f'Quaternion{" " if self.quaternion.shape == (4,) else "s of shape "+str(self.quaternion.shape[:-1])+chr(10)}'\
               + str(self.quaternion)


    def __copy__(self: MyType,
                 rotation: Union[None, FloatSequence, 'Rotation'] = None) -> MyType:
        """
        Return deepcopy(self).

        Create deep copy.

        """
        dup = copy.deepcopy(self)
        if rotation is not None:
            dup.quaternion = Rotation(rotation).quaternion
        return dup

    copy = __copy__


    def __getitem__(self,
                    item: Union[Tuple[Union[None, int, slice]], int, bool, np.bool_, np.ndarray]):
        """
        Return self[item].

        Return slice according to item.

        """
        return self.copy() if self.shape == () else \
               self.copy(self.quaternion[item+(slice(None),)] if isinstance(item,tuple) else self.quaternion[item])


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : Rotation
            Rotation to check for equality.

        """
        return NotImplemented if not isinstance(other, Rotation) else \
               np.logical_or(np.all(self.quaternion ==     other.quaternion,axis=-1),
                             np.all(self.quaternion == -1.*other.quaternion,axis=-1))


    def __ne__(self,
               other: object) -> bool:
        """
        Return self!=other.

        Test inequality of other.

        Parameters
        ----------
        other : Rotation
            Rotation to check for inequality.

        """
        return np.logical_not(self==other) if isinstance(other, Rotation) else NotImplemented

    def isclose(self: MyType,
                other: MyType,
                rtol: float = 1.e-5,
                atol: float = 1.e-8,
                equal_nan: bool = True) -> bool:
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
        mask : numpy.ndarray of bool, shape (self.shape)
            Mask indicating where corresponding rotations are close.

        """
        s = self.quaternion
        o = other.quaternion
        return np.logical_or(np.all(np.isclose(s,    o,rtol,atol,equal_nan),axis=-1),
                             np.all(np.isclose(s,-1.*o,rtol,atol,equal_nan),axis=-1))


    def allclose(self: MyType,
                 other: MyType,
                 rtol: float = 1.e-5,
                 atol: float = 1.e-8,
                 equal_nan: bool = True) -> Union[np.bool_, bool]:
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
    def size(self) -> int:
        return self.quaternion[...,0].size

    @property
    def shape(self) -> Tuple[int, ...]:
        return self.quaternion[...,0].shape


    def __len__(self) -> int:
        """
        Return len(self).

        Length of leading/leftmost dimension of array.

        """
        return 0 if self.shape == () else self.shape[0]


    def __invert__(self: MyType) -> MyType:
        """
        Return ~self.

        Inverse rotation (backward rotation).

        """
        dup = self.copy()
        dup.quaternion[...,1:] *= -1.
        return dup


    def __pow__(self: MyType,
                exp: Union[float, int]) -> MyType:
        """
        Return self**exp.

        Perform the rotation 'exp' times.

        Parameters
        ----------
        exp : float
            Exponent.

        """
        phi = np.arccos(self.quaternion[...,0:1])
        p = self.quaternion[...,1:]/np.linalg.norm(self.quaternion[...,1:],axis=-1,keepdims=True)
        return self.copy(Rotation(np.block([np.cos(exp*phi),np.sin(exp*phi)*p]))._standardize())

    def __ipow__(self: MyType,
                 exp: Union[float, int]) -> MyType:
        """
        Return self**=exp.

        Perform the rotation 'exp' times (in-place).

        Parameters
        ----------
        exp : float
            Exponent.

        """
        return self**exp


    def __mul__(self: MyType,
                other: MyType) -> MyType:
        """
        Return self*other.

        Compose with other.

        Parameters
        ----------
        other : Rotation, shape (self.shape)
            Rotation for composition.
            Compatible innermost dimensions will blend.

        Returns
        -------
        composition : Rotation
            Compound rotation self*other, i.e. first other then self rotation.

        """
        if isinstance(other,Rotation):
            blend = util.shapeblender( self.shape,other.shape)
            s_m   = util.shapeshifter( self.shape,blend,mode='right')
            s_o   = util.shapeshifter(other.shape,blend,mode='left')

            q_m = self.broadcast_to(s_m).quaternion[...,0:1]
            p_m = self.broadcast_to(s_m).quaternion[...,1:]
            q_o = other.broadcast_to(s_o).quaternion[...,0:1]
            p_o = other.broadcast_to(s_o).quaternion[...,1:]

            qmo = q_m*q_o
            q = (qmo - np.einsum('...i,...i',p_m,p_o).reshape(qmo.shape))
            p = q_m*p_o + q_o*p_m + _P * np.cross(p_m,p_o)
            return self.copy(Rotation(np.block([q,p]))._standardize())
        else:
            raise TypeError('use "R@b", i.e. matmul, to apply rotation "R" to object "b"')

    def __imul__(self: MyType,
                 other: MyType) -> MyType:
        """
        Return self*=other.

        Compose with other (in-place).

        Parameters
        ----------
        other : Rotation, shape (self.shape)
            Rotation for composition.
            Compatible innermost dimensions will blend.

        """
        return self*other


    def __truediv__(self: MyType,
                    other: MyType) -> MyType:
        """
        Return self/other.

        Compose with inverse of other.

        Parameters
        ----------
        other : damask.Rotation, shape (self.shape)
            Rotation to invert for composition.
            Compatible innermost dimensions will blend.

        Returns
        -------
        composition : Rotation
            Compound rotation self*(~other), i.e. first inverse of other then self rotation.

        """
        if isinstance(other,Rotation):
            return self*~other
        else:
            raise TypeError('use "R@b", i.e. matmul, to apply rotation "R" to object "b"')

    def __itruediv__(self: MyType,
                     other: MyType) -> MyType:
        """
        Return self/=other.

        Compose with inverse of other (in-place).

        Parameters
        ----------
        other : Rotation, shape (self.shape)
            Rotation to invert for composition.

        """
        return self/other


    def __matmul__(self,
                   other: np.ndarray) -> np.ndarray:
        """
        Return self@other.

        Rotate vector, second-order tensor, or fourth-order tensor.
        `other` is interpreted as an array of tensor quantities with the highest-possible order
        considering the shape of `self`. Compatible innermost dimensions will blend.
        For instance, shapes of (2,) and (3,3) for `self` and `other` prompt interpretation of
        `other` as a second-rank tensor and result in (2,) rotated tensors, whereas
        shapes of (2,1) and (3,3) for `self` and `other` result in (2,3) rotated vectors.

        Parameters
        ----------
        other : numpy.ndarray, shape (...,3), (...,3,3), or (...,3,3,3,3)
            Vector or tensor on which to apply the rotation.

        Returns
        -------
        rotated : numpy.ndarray, shape (...,3), (...,3,3), or (...,3,3,3,3)
            Rotated vector or tensor, i.e. transformed to frame defined by rotation.

        Examples
        --------
        All below examples rely on imported modules:
        >>> import numpy as np
        >>> import damask

        Application of twelve (random) rotations to a set of five vectors.

        >>> r = damask.Rotation.from_random(shape=(12))
        >>> o = np.ones((5,3))
        >>> (r@o).shape                                    # (12) @ (5, 3)
        (12,5, 3)

        Application of a (random) rotation to all twelve second-rank tensors.

        >>> r = damask.Rotation.from_random()
        >>> o = np.ones((12,3,3))
        >>> (r@o).shape                                    # (1) @ (12, 3,3)
        (12,3,3)

        Application of twelve (random) rotations to the corresponding twelve second-rank tensors.

        >>> r = damask.Rotation.from_random(shape=(12))
        >>> o = np.ones((12,3,3))
        >>> (r@o).shape                                    # (12) @ (3,3)
        (12,3,3)

        Application of each of three (random) rotations to all three vectors.

        >>> r = damask.Rotation.from_random(shape=(3))
        >>> o = np.ones((3,3))
        >>> (r[...,np.newaxis]@o[np.newaxis,...]).shape    # (3,1) @ (1,3, 3)
        (3,3,3)

        Application of twelve (random) rotations to all twelve second-rank tensors.

        >>> r = damask.Rotation.from_random(shape=(12))
        >>> o = np.ones((12,3,3))
        >>> (r@o[np.newaxis,...]).shape                    # (12) @ (1,12, 3,3)
        (12,3,3,3)

        """
        if isinstance(other, np.ndarray):
            obs = util.shapeblender(self.shape,other.shape)[len(self.shape):]
            for l in [4,2,1]:
                if obs[-l:] == l*(3,):
                    bs = util.shapeblender(self.shape,other.shape[:-l],False)
                    self_ = self.broadcast_to(bs) if self.shape != bs else self
                    if l==1:
                        q_m = self_.quaternion[...,0]
                        p_m = self_.quaternion[...,1:]
                        A = q_m**2 - np.einsum('...i,...i',p_m,p_m)
                        B = 2. * np.einsum('...i,...i',p_m,other)
                        C = 2. * _P * q_m
                        return np.block([(A * other[...,i]) +
                                         (B *   p_m[...,i]) +
                                         (C * ( p_m[...,(i+1)%3]*other[...,(i+2)%3]
                                              - p_m[...,(i+2)%3]*other[...,(i+1)%3]))
                                        for i in [0,1,2]]).reshape(bs+(3,),order='F')
                    else:
                        return np.einsum({2: '...im,...jn,...mn',
                                          4: '...im,...jn,...ko,...lp,...mnop'}[l],
                                         *l*[self_.as_matrix()],
                                         other)
            raise ValueError('can only rotate vectors, second-order tensors, and fourth-order tensors')
        elif isinstance(other, Rotation):
            raise TypeError('use "R2*R1", i.e. multiplication, to compose rotations "R1" and "R2"')
        else:
            raise TypeError(f'cannot rotate "{type(other)}"')

    apply = __matmul__


    def _standardize(self: MyType) -> MyType:
        """Standardize quaternion (ensure positive real hemisphere)."""
        self.quaternion[self.quaternion[...,0] < 0.] *= -1.
        return self


    def append(self: MyType,
               other: Union[MyType, List[MyType]]) -> MyType:
        """
        Extend array along first dimension with other array(s).

        Parameters
        ----------
        other : (list of) damask.Rotation

        """
        return self.copy(np.vstack(tuple(map(lambda x:x.quaternion,
                                             [self]+other if isinstance(other,list) else [self,other]))))


    def flatten(self: MyType,
                order: Literal['C','F','A'] = 'C') -> MyType:
        """
        Flatten array.

        Parameters
        ----------
        order : {'C', 'F', 'A'}, optional
            'C' flattens in row-major (C-style) order.
            'F' flattens in column-major (Fortran-style) order.
            'A' flattens in column-major order if object is Fortran contiguous in memory,
            row-major order otherwise.
            Defaults to 'C'.

        Returns
        -------
        flattened : damask.Rotation
            Rotation flattened to single dimension.

        """
        return self.copy(self.quaternion.reshape((-1,4),order=order))


    def reshape(self: MyType,
                shape: Union[int, IntSequence],
                order: Literal['C','F','A'] = 'C') -> MyType:
        """
        Reshape array.

        Parameters
        ----------
        shape : (sequence of) int
            New shape, number of elements needs to match the original shape.
            If an integer is supplied, then the result will be a 1-D array of that length.
        order : {'C', 'F', 'A'}, optional
            'C' flattens in row-major (C-style) order.
            'F' flattens in column-major (Fortran-style) order.
            'A' flattens in column-major order if object is Fortran contiguous in memory,
            row-major order otherwise.
            Defaults to 'C'.

        Returns
        -------
        reshaped : damask.Rotation
            Rotation of given shape.

        """
        if isinstance(shape,(int,np.integer)): shape = (shape,)
        return self.copy(self.quaternion.reshape(tuple(shape)+(4,),order=order))


    def broadcast_to(self: MyType,
                     shape: Union[int, IntSequence],
                     mode: Literal['left', 'right'] = 'right') -> MyType:
        """
        Broadcast array.

        Parameters
        ----------
        shape : (sequence of) int
            Shape of broadcasted array, needs to be compatible with the original shape.
        mode : str, optional
            Where to preferentially locate missing dimensions.
            Either 'left' or 'right' (default).

        Returns
        -------
        broadcasted : damask.Rotation
            Rotation broadcasted to given shape.

        """
        shape_ = (shape,) if isinstance(shape,(int,np.integer)) else tuple(shape)
        return self.copy(np.broadcast_to(self.quaternion.reshape(util.shapeshifter(self.shape,shape_,mode)+(4,)),
                                         shape_+(4,)))


    def average(self: MyType,
                weights: Optional[FloatSequence] = None) -> MyType:
        """
        Average along last array dimension.

        Parameters
        ----------
        weights : numpy.ndarray, shape (self.shape), optional
            Relative weight of each rotation.

        Returns
        -------
        average : damask.Rotation
            Weighted average of original Rotation field.

        References
        ----------
        F. Landis Markley et al., Journal of Guidance, Control, and Dynamics 30(4):1193-1197, 2007
        https://doi.org/10.2514/1.28949

        """
        def _M(quat):
            """Intermediate representation supporting quaternion averaging."""
            return np.einsum('...i,...j',quat,quat)

        weights_ = np.ones(self.shape,dtype=float) if weights is None else np.array(weights,float)

        eig, vec = np.linalg.eig(np.sum(_M(self.quaternion) * weights_[...,np.newaxis,np.newaxis],axis=-3)
                                /np.sum(                      weights_[...,np.newaxis,np.newaxis],axis=-3))

        return self.copy(Rotation.from_quaternion(np.real(
                                                  np.squeeze(
                                                  np.take_along_axis(vec,
                                                                     eig.argmax(axis=-1)[...,np.newaxis,np.newaxis],
                                                                     axis=-1),
                                                  axis=-1)),
                                                  accept_homomorph = True))


    def misorientation(self: MyType,
                       other: MyType) -> MyType:
        """
        Calculate misorientation to other Rotation.

        Parameters
        ----------
        other : damask.Rotation
            Rotation to which the misorientation is computed.
            Compatible innermost dimensions will blend.

        Returns
        -------
        g : damask.Rotation
            Misorientation.

        """
        return ~(self*~other)


    ################################################################################################
    # convert to different orientation representations (numpy arrays)

    def as_quaternion(self) -> np.ndarray:
        """
        Represent as unit quaternion.

        Returns
        -------
        q : numpy.ndarray, shape (...,4)
            Unit quaternion (q_0, q_1, q_2, q_3) in positive real hemisphere, i.e. ǀqǀ = 1, q_0 ≥ 0.

        """
        return self.quaternion.copy()

    def as_Euler_angles(self,
                        degrees: bool = False) -> np.ndarray:
        """
        Represent as Bunge Euler angles.

        Parameters
        ----------
        degrees : bool, optional
            Return angles in degrees. Defaults to False.

        Returns
        -------
        phi : numpy.ndarray, shape (...,3)
            Bunge Euler angles (φ_1 ∈ [0,2π], ϕ ∈ [0,π], φ_2 ∈ [0,2π])
            or (φ_1 ∈ [0,360], ϕ ∈ [0,180], φ_2 ∈ [0,360]) if degrees == True.

        Notes
        -----
        Bunge Euler angles correspond to a rotation axis sequence of z–x'–z''.

        Examples
        --------
        Cube orientation as Bunge Euler angles.

        >>> import damask
        >>> damask.Rotation([1,0,0,0]).as_Euler_angles()
        array([0., 0., 0.])

        """
        eu = Rotation._qu2eu(self.quaternion)
        return np.degrees(eu) if degrees else eu

    def as_axis_angle(self,
                      degrees: bool = False,
                      pair: bool = False) -> Union[Tuple[np.ndarray, np.ndarray], np.ndarray]:
        """
        Represent as axis–angle pair.

        Parameters
        ----------
        degrees : bool, optional
            Return rotation angle in degrees. Defaults to False.
        pair : bool, optional
            Return tuple of axis and angle. Defaults to False.

        Returns
        -------
        n_omega : numpy.ndarray, shape (...,4) or tuple ((...,3), (...)) if pair == True
            Axis and angle [n_1, n_2, n_3, ω] with ǀnǀ = 1 and ω ∈ [0,π]
            or ω ∈ [0,180] if degrees == True.

        Examples
        --------
        Cube orientation as axis–angle pair.

        >>> import damask
        >>> damask.Rotation([1,0,0,0]).as_axis_angle(pair=True)
        (array([0., 0., 1.]), array(0.))

        """
        ax: np.ndarray = Rotation._qu2ax(self.quaternion)
        if degrees: ax[...,3] = np.degrees(ax[...,3])
        return (ax[...,:3],ax[...,3]) if pair else ax

    def as_matrix(self) -> np.ndarray:
        """
        Represent as rotation matrix.

        Returns
        -------
        R : numpy.ndarray, shape (...,3,3)
            Rotation matrix R with det(R) = 1, R.T ∙ R = I.

        Examples
        --------
        Cube orientation as rotation matrix.

        >>> import damask
        >>> damask.Rotation([1,0,0,0]).as_matrix()
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])

        """
        return Rotation._qu2om(self.quaternion)

    def as_Rodrigues_vector(self,
                            compact: bool = False) -> np.ndarray:
        """
        Represent as Rodrigues–Frank vector with separate axis and angle argument.

        Parameters
        ----------
        compact : bool, optional
            Return three-component Rodrigues–Frank vector,
            i.e. axis and angle argument are not separated.

        Returns
        -------
        rho : numpy.ndarray, shape (...,4) or (...,3) if compact == True
            Rodrigues–Frank vector [n_1, n_2, n_3, tan(ω/2)] with ǀnǀ = 1 and ω ∈ [0,π]
            or [n_1, n_2, n_3] with ǀnǀ = tan(ω/2) and ω ∈ [0,π] if compact == True.

        Examples
        --------
        Cube orientation as three-component Rodrigues–Frank vector.

        >>> import damask
        >>> damask.Rotation([1,0,0,0]).as_Rodrigues_vector(compact=True)
        array([ 0.,  0., 0.])

        """
        ro = Rotation._qu2ro(self.quaternion)
        if compact:
            with np.errstate(invalid='ignore'):
                return ro[...,:3]*ro[...,3:4]
        else:
            return ro

    def as_homochoric(self) -> np.ndarray:
        """
        Represent as homochoric vector.

        Returns
        -------
        h : numpy.ndarray, shape (...,3)
            Homochoric vector (h_1, h_2, h_3) with ǀhǀ < (3/4*π)^(1/3).

        Examples
        --------
        Cube orientation as homochoric vector.

        >>> import damask
        >>> damask.Rotation([1,0,0,0]).as_homochoric()
        array([0., 0., 0.])

        """
        return Rotation._qu2ho(self.quaternion)

    def as_cubochoric(self) -> np.ndarray:
        """
        Represent as cubochoric vector.

        Returns
        -------
        x : numpy.ndarray, shape (...,3)
            Cubochoric vector (x_1, x_2, x_3) with max(x_i) < 1/2*π^(2/3).

        Examples
        --------
        Cube orientation as cubochoric vector.

        >>> import damask
        >>> damask.Rotation([1,0,0,0]).as_cubochoric()
        array([0., 0., 0.])

        """
        return Rotation._qu2cu(self.quaternion)

    ################################################################################################
    # Static constructors. The input data needs to follow the conventions, options allow to
    # relax the conventions.
    @staticmethod
    def from_quaternion(q: Union[Sequence[FloatSequence], np.ndarray],
                        accept_homomorph: bool = False,
                        normalize: bool = False,
                        P: Literal[1, -1] = -1) -> 'Rotation':
        """
        Initialize from quaternion.

        Parameters
        ----------
        q : numpy.ndarray, shape (...,4)
            Unit quaternion (q_0, q_1, q_2, q_3) in positive real hemisphere, i.e. ǀqǀ = 1 and q_0 ≥ 0.
        accept_homomorph : bool, optional
            Allow homomorphic variants, i.e. q_0 < 0 (negative real hemisphere).
            Defaults to False.
        normalize: bool, optional
            Allow ǀqǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.

        Returns
        -------
        new : damask.Rotation

        """
        qu = np.array(q,dtype=float)
        if qu.shape[:-2:-1] != (4,): raise ValueError('invalid shape')
        if abs(P) != 1: raise ValueError('P ∉ {-1,1}')

        qu[...,1:4] *= -P

        if accept_homomorph:
            qu[qu[...,0]<0.] *= -1.
        elif np.any(qu[...,0] < 0.):
            raise ValueError('quaternion with negative first (real) component')
        if normalize:
            qu /= np.linalg.norm(qu,axis=-1,keepdims=True)
        elif not np.allclose(np.linalg.norm(qu,axis=-1),1.,rtol=1.e-8):
            raise ValueError('quaternion is not of unit length')

        return Rotation(qu)

    @staticmethod
    def from_Euler_angles(phi: np.ndarray,
                          degrees: bool = False) -> 'Rotation':
        """
        Initialize from Bunge Euler angles.

        Parameters
        ----------
        phi : numpy.ndarray, shape (...,3)
            Euler angles (φ_1 ∈ [0,2π], ϕ ∈ [0,π], φ_2 ∈ [0,2π])
            or (φ_1 ∈ [0,360], ϕ ∈ [0,180], φ_2 ∈ [0,360]) if degrees == True.
        degrees : bool, optional
            Euler angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Rotation

        Notes
        -----
        Bunge Euler angles correspond to a rotation axis sequence of z–x'–z''.

        """
        eu = np.array(phi,dtype=float)
        if eu.shape[:-2:-1] != (3,): raise ValueError('invalid shape')

        eu = np.radians(eu) if degrees else eu
        if np.any(eu < 0.) or np.any(eu > np.pi*np.array([2.,1.,2.])):
            raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π]')

        return Rotation(Rotation._eu2qu(eu))

    @staticmethod
    def from_axis_angle(n_omega: np.ndarray,
                        degrees: bool = False,
                        normalize: bool = False,
                        P: Literal[1, -1] = -1) -> 'Rotation':
        """
        Initialize from axis–angle pair.

        Parameters
        ----------
        n_omega : numpy.ndarray, shape (...,4)
            Axis and angle (n_1, n_2, n_3, ω) with ǀnǀ = 1 and ω ∈ [0,π]
            or ω ∈ [0,180] if degrees == True.
        degrees : bool, optional
            Angle ω is given in degrees. Defaults to False.
        normalize: bool, optional
            Allow ǀnǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.

        Returns
        -------
        new : damask.Rotation

        """
        ax = np.array(n_omega,dtype=float)
        if ax.shape[:-2:-1] != (4,): raise ValueError('invalid shape')
        if abs(P) != 1: raise ValueError('P ∉ {-1,1}')

        ax[...,0:3] *= -P
        if degrees: ax[...,  3] = np.radians(ax[...,3])
        if np.any(ax[...,3] < 0.) or np.any(ax[...,3] > np.pi):
            raise ValueError('axis–angle rotation angle outside of [0..π]')

        if normalize:
            ax[...,0:3] /= np.linalg.norm(ax[...,0:3],axis=-1,keepdims=True)
        elif not np.allclose(np.linalg.norm(ax[...,0:3],axis=-1),1.):
            raise ValueError('axis–angle rotation axis is not of unit length')

        return Rotation(Rotation._ax2qu(ax))

    @staticmethod
    def from_basis(basis: np.ndarray,
                   orthonormal: bool = True,
                   reciprocal: bool = False) -> 'Rotation':
        """
        Initialize from basis vector triplet.

        Parameters
        ----------
        basis : numpy.ndarray, shape (...,3,3)
            Three three-dimensional basis vectors.
        orthonormal : bool, optional
            Basis is strictly orthonormal, i.e. is free of stretch components. Defaults to True.
        reciprocal : bool, optional
            Basis vectors are given in reciprocal (instead of real) space. Defaults to False.

        Returns
        -------
        new : damask.Rotation

        """
        om = np.array(basis,dtype=float)
        if om.shape[-2:] != (3,3): raise ValueError('invalid shape')

        if reciprocal:
            om = np.linalg.inv(tensor.transpose(om)/np.pi)                                          # transform reciprocal basis set
            orthonormal = False                                                                     # contains stretch

        if not orthonormal:
            (U,S,Vh) = np.linalg.svd(om)                                                            # singular value decomposition
            om = np.einsum('...ij,...jl',U,Vh)
        elif  not np.allclose(np.einsum('...i,...i',om[...,0],om[...,1]),0.) \
           or not np.allclose(np.einsum('...i,...i',om[...,1],om[...,2]),0.) \
           or not np.allclose(np.einsum('...i,...i',om[...,2],om[...,0]),0.):
            raise ValueError('orientation matrix is not orthogonal')

        if not np.allclose(np.linalg.det(om),1.):
            raise ValueError('orientation matrix has determinant ≠ 1')

        return Rotation(Rotation._om2qu(om))

    @staticmethod
    def from_matrix(R: np.ndarray,
                    normalize: bool = False) -> 'Rotation':
        """
        Initialize from rotation matrix.

        Parameters
        ----------
        R : numpy.ndarray, shape (...,3,3)
            Rotation matrix with det(R) = 1 and R.T ∙ R = I.
        normalize : bool, optional
            Rescales rotation matrix to unit determinant. Defaults to False.

        Returns
        -------
        new : damask.Rotation

        """
        return Rotation.from_basis(np.array(R,dtype=float) * (np.linalg.det(R)**(-1./3.))[...,np.newaxis,np.newaxis]
                                   if normalize else
                                   R)

    @staticmethod
    def from_parallel(a: np.ndarray,
                      b: np.ndarray ) -> 'Rotation':
        """
        Initialize from pairs of two orthogonal basis vectors.

        Parameters
        ----------
        a : numpy.ndarray, shape (...,2,3)
            Two three-dimensional vectors of first orthogonal basis.
        b : numpy.ndarray, shape (...,2,3)
            Corresponding three-dimensional vectors of second basis.

        Returns
        -------
        new : damask.Rotation

        """
        a_ = np.array(a,dtype=float)
        b_ = np.array(b,dtype=float)
        if a_.shape[-2:] != (2,3) or b_.shape[-2:] != (2,3) or a_.shape != b_.shape:
            raise ValueError('invalid shape')
        am = np.stack([          a_[...,0,:],
                                             a_[...,1,:],
                        np.cross(a_[...,0,:],a_[...,1,:]) ],axis=-2)
        bm = np.stack([          b_[...,0,:],
                                             b_[...,1,:],
                        np.cross(b_[...,0,:],b_[...,1,:]) ],axis=-2)

        return              Rotation.from_basis(np.swapaxes(am/np.linalg.norm(am,axis=-1,keepdims=True),-1,-2))\
            .misorientation(Rotation.from_basis(np.swapaxes(bm/np.linalg.norm(bm,axis=-1,keepdims=True),-1,-2)))


    @staticmethod
    def from_Rodrigues_vector(rho: np.ndarray,
                              normalize: bool = False,
                              P: Literal[1, -1] = -1) -> 'Rotation':
        """
        Initialize from Rodrigues–Frank vector (with angle separated from axis).

        Parameters
        ----------
        rho : numpy.ndarray, shape (...,4)
            Rodrigues–Frank vector (n_1, n_2, n_3, tan(ω/2)) with ǀnǀ = 1  and ω ∈ [0,π].
        normalize : bool, optional
            Allow ǀnǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.

        Returns
        -------
        new : damask.Rotation

        """
        ro = np.array(rho,dtype=float)
        if ro.shape[:-2:-1] != (4,): raise ValueError('invalid shape')
        if abs(P) != 1: raise ValueError('P ∉ {-1,1}')

        ro[...,0:3] *= -P
        if np.any(ro[...,3] < 0.): raise ValueError('Rodrigues vector rotation angle is negative')

        if normalize:
            ro[...,0:3] /= np.linalg.norm(ro[...,0:3],axis=-1,keepdims=True)
        elif not np.allclose(np.linalg.norm(ro[...,0:3],axis=-1),1.):
            raise ValueError('Rodrigues vector rotation axis is not of unit length')

        return Rotation(Rotation._ro2qu(ro))

    @staticmethod
    def from_homochoric(h: np.ndarray,
                        P: Literal[1, -1] = -1) -> 'Rotation':
        """
        Initialize from homochoric vector.

        Parameters
        ----------
        h : numpy.ndarray, shape (...,3)
            Homochoric vector (h_1, h_2, h_3) with ǀhǀ < (3/4*π)^(1/3).
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.

        Returns
        -------
        new : damask.Rotation

        """
        ho = np.array(h,dtype=float)
        if ho.shape[:-2:-1] != (3,): raise ValueError('invalid shape')
        if abs(P) != 1: raise ValueError('P ∉ {-1,1}')

        ho *= -P

        if np.any(np.linalg.norm(ho,axis=-1) > _R1+1.e-9):
            raise ValueError('homochoric coordinate outside of the sphere')

        return Rotation(Rotation._ho2qu(ho))

    @staticmethod
    def from_cubochoric(x: np.ndarray,
                        P: Literal[1, -1] = -1) -> 'Rotation':
        """
        Initialize from cubochoric vector.

        Parameters
        ----------
        x : numpy.ndarray, shape (...,3)
            Cubochoric vector (x_1, x_2, x_3) with max(x_i) < 1/2*π^(2/3).
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.

        Returns
        -------
        new : damask.Rotation

        """
        cu = np.array(x,dtype=float)
        if cu.shape[:-2:-1] != (3,): raise ValueError('invalid shape')
        if abs(P) != 1: raise ValueError('P ∉ {-1,1}')
        if np.abs(np.max(cu)) > np.pi**(2./3.) * 0.5+1.e-9:
            raise ValueError('cubochoric coordinate outside of the cube')

        ho = -P * Rotation._cu2ho(cu)

        return Rotation(Rotation._ho2qu(ho))


    @staticmethod
    def from_random(shape: Union[None, int, IntSequence] = None,
                    rng_seed: Optional[NumpyRngSeed] = None) -> 'Rotation':
        """
        Initialize with samples from a uniform distribution.

        Parameters
        ----------
        shape : (sequence of) int, optional
            Shape of the returned array. Defaults to None, which gives a scalar.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.

        Returns
        -------
        new : damask.Rotation

        """
        rng = np.random.default_rng(rng_seed)
        r = rng.random(3 if shape is None else tuple(shape)+(3,) if hasattr(shape, '__iter__') else (shape,3)) # type: ignore

        A = np.sqrt(r[...,2])
        B = np.sqrt(1.-r[...,2])
        q = np.stack([np.cos(2.*np.pi*r[...,0])*A,
                      np.sin(2.*np.pi*r[...,1])*B,
                      np.cos(2.*np.pi*r[...,1])*B,
                      np.sin(2.*np.pi*r[...,0])*A],axis=-1)

        return Rotation(q if shape is None else q.reshape(r.shape[:-1]+(4,)))._standardize()


    @staticmethod
    def from_ODF(weights: np.ndarray,
                 phi: np.ndarray,
                 shape: Union[None, int, IntSequence] = None,
                 degrees: bool = False,
                 fractions: bool = True,
                 rng_seed: Optional[NumpyRngSeed] = None) -> 'Rotation':
        """
        Initialize with samples from a binned orientation distribution function (ODF).

        Parameters
        ----------
        weights : numpy.ndarray, shape (n)
            Texture intensity values (probability density or volume fraction) at Euler space grid points.
        phi : numpy.ndarray, shape (n,3)
            Grid coordinates in Euler space at which weights are defined.
        shape : (sequence of) int, optional
            Shape of the returned array. Defaults to None, which gives a scalar.
        degrees : bool, optional
            Euler space grid coordinates are in degrees. Defaults to True.
        fractions : bool, optional
            ODF values correspond to volume fractions, not probability densities.
            Defaults to True.
        rng_seed: {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.

        Returns
        -------
        new : damask.Rotation

        Notes
        -----
        Due to the distortion of Euler space in the vicinity of ϕ = 0, probability densities, p, defined on
        grid points with ϕ = 0 will never result in reconstructed orientations as their dV/V = p dγ = p × 0.
        Hence, it is recommended to transform any such dataset to a cell-centered version, which avoids grid points at ϕ = 0.

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
            return delta[0]*2.*np.sin(delta[1]/2.)*delta[2] / 8. / np.pi**2 * np.sin(np.radians(eu[:,1]) if deg else eu[:,1])

        dg = 1. if fractions else _dg(phi,degrees)
        dV_V = dg * np.maximum(0.,weights.squeeze())

        N = 1 if shape is None else np.prod(shape).astype(int)
        return Rotation.from_Euler_angles(phi[util.hybrid_IA(dV_V,N,rng_seed)],degrees).reshape(() if shape is None else shape)


    @staticmethod
    def from_spherical_component(center: 'Rotation',
                                 sigma: float,
                                 shape: Union[None, int, IntSequence] = None,
                                 degrees: bool = False,
                                 rng_seed: Optional[NumpyRngSeed] = None) -> 'Rotation':
        """
        Initialize with samples from a Gaussian distribution around a given center.

        Parameters
        ----------
        center : Rotation or Orientation
            Central rotation.
        sigma : float
            Standard deviation of (Gaussian) misorientation distribution.
        shape : (sequence of) int, optional
            Shape of the returned array. Defaults to None, which gives a scalar.
        degrees : bool, optional
            sigma is given in degrees. Defaults to False.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.

        Examples
        --------
        Create a brass texture consisting of
        200 orientations:

        >>> import damask
        >>> center = damask.Rotation.from_Euler_angles([35.,45.,0.],degrees=True)
        >>> brass = damask.Rotation.from_spherical_component(center=center,sigma=3.,shape=200,degrees=True)

        Create a Goss texture consisting of
        100 orientations:

        >>> import damask
        >>> center = damask.Rotation.from_Euler_angles([0.,45.,0.],degrees=True)
        >>> goss = damask.Rotation.from_spherical_component(center=center,sigma=3.,shape=100,degrees=True)

        """
        rng = np.random.default_rng(rng_seed)
        sigma = np.radians(sigma) if degrees else sigma
        N = 1 if shape is None else np.prod(shape)
        u,Theta  = (rng.random((N,2)) * 2. * np.array([1.,np.pi]) - np.array([1.,0.])).T
        omega = abs(rng.normal(scale=sigma,size=N))
        p = np.column_stack([np.sqrt(1.-u**2)*np.cos(Theta),
                             np.sqrt(1.-u**2)*np.sin(Theta),
                             u, omega])

        return Rotation.from_axis_angle(p).reshape(() if shape is None else shape) * center


    @staticmethod
    def from_fiber_component(crystal: IntSequence,
                             sample: IntSequence,
                             sigma: float = 0.,
                             shape: Union[None, int, IntSequence] = None,
                             degrees: bool = False,
                             rng_seed: Optional[NumpyRngSeed] = None) -> 'Rotation':
        """
        Initialize with samples from a Gaussian distribution around a given direction.

        Parameters
        ----------
        crystal : numpy.ndarray, shape (2)
            Polar coordinates (polar angle θ from [0 0 1], azimuthal angle φ from [1 0 0])
            of fiber direction in crystal frame.
        sample : numpy.ndarray, shape (2)
            Polar coordinates (polar angle θ from z, azimuthal angle φ from x)
            of fiber direction in sample frame.
        sigma : float, optional
            Standard deviation of (Gaussian) misorientation distribution.
            Defaults to 0.
        shape : (sequence of) int, optional
            Shape of the returned array. Defaults to None, which gives a scalar.
        degrees : bool, optional
            sigma and polar coordinates are given in degrees. Defaults to False.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.

        Returns
        -------
        new : damask.Rotation

        Notes
        -----
        The crystal direction for (θ=0,φ=0) is [0 0 1],
        the sample direction for (θ=0,φ=0) is z.

        Polar coordinates follow the ISO 80000-2:2019 convention
        typically used in physics.
        See https://en.wikipedia.org/wiki/Spherical_coordinate_system.

        Ranges 0≤θ≤π and 0≤φ≤2π give a unique set of coordinates.

        Examples
        --------
        Create an ideal α-fiber texture (<1 1 0> ǀǀ RD=x) consisting of
        200 orientations:

        >>> import damask
        >>> import numpy as np
        >>> alpha = damask.Rotation.from_fiber_component([np.pi/4.,0.],[np.pi/2.,0.],shape=200)

        Create an ideal γ-fiber texture (<1 1 1> ǀǀ ND=z) consisting of
        100 orientations:

        >>> import damask
        >>> gamma = damask.Rotation.from_fiber_component([54.7,45.],[0.,0.],shape=100,degrees=True)

        """
        rng = np.random.default_rng(rng_seed)
        sigma_,alpha,beta = (np.radians(coordinate) for coordinate in (sigma,crystal,sample)) if degrees else \
                             map(np.array, (sigma,crystal,sample))

        d_cr  = np.array([np.sin(alpha[0])*np.cos(alpha[1]), np.sin(alpha[0])*np.sin(alpha[1]), np.cos(alpha[0])])
        d_lab = np.array([np.sin( beta[0])*np.cos( beta[1]), np.sin( beta[0])*np.sin( beta[1]), np.cos( beta[0])])
        ax_align = np.append(np.cross(d_lab,d_cr), np.arccos(np.dot(d_lab,d_cr)))
        if np.isclose(ax_align[3],0.): ax_align[:3] = np.array([1.,0.,0.])
        R_align  = Rotation.from_axis_angle(ax_align if ax_align[3] > 0. else -ax_align,normalize=True) # rotate fiber axis from sample to crystal frame

        N = 1 if shape is None else np.prod(shape).astype(int)
        u,Theta  = (rng.random((N,2)) * 2. * np.array([1.,np.pi]) - np.array([1.,0.])).T
        omega  = abs(rng.normal(scale=sigma_,size=N))
        p = np.column_stack([np.sqrt(1.-u**2)*np.cos(Theta),
                             np.sqrt(1.-u**2)*np.sin(Theta),
                             u, omega])
        p[:,:3] = np.einsum('ij,...j',np.eye(3)-np.outer(d_lab,d_lab),p[:,:3])                      # remove component along fiber axis
        f = np.column_stack((np.broadcast_to(d_lab,(N,3)),rng.random(N)*np.pi))
        f[::2,:3] *= -1.                                                                            # flip half the rotation axes to negative sense

        return (R_align.broadcast_to(N)
              * Rotation.from_axis_angle(p,normalize=True)
              * Rotation.from_axis_angle(f)).reshape(() if shape is None else shape)


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
    def _qu2om(qu: np.ndarray) -> np.ndarray:
        qq = qu[...,0:1]**2-(qu[...,1:2]**2 + qu[...,2:3]**2 + qu[...,3:4]**2)
        om = np.block([qq + 2.*qu[...,1:2]**2,
                       2.*(qu[...,2:3]*qu[...,1:2]-_P*qu[...,0:1]*qu[...,3:4]),
                       2.*(qu[...,3:4]*qu[...,1:2]+_P*qu[...,0:1]*qu[...,2:3]),
                       2.*(qu[...,1:2]*qu[...,2:3]+_P*qu[...,0:1]*qu[...,3:4]),
                       qq + 2.*qu[...,2:3]**2,
                       2.*(qu[...,3:4]*qu[...,2:3]-_P*qu[...,0:1]*qu[...,1:2]),
                       2.*(qu[...,1:2]*qu[...,3:4]-_P*qu[...,0:1]*qu[...,2:3]),
                       2.*(qu[...,2:3]*qu[...,3:4]+_P*qu[...,0:1]*qu[...,1:2]),
                       qq + 2.*qu[...,3:4]**2,
                      ]).reshape(qu.shape[:-1]+(3,3))
        return om

    @staticmethod
    def _qu2eu(qu: np.ndarray) -> np.ndarray:
        """
        Quaternion to Bunge Euler angles.

        References
        ----------
        E. Bernardes and S. Viollet, PLoS ONE 17(11):e0276302, 2022
        https://doi.org/10.1371/journal.pone.0276302

        """
        a =     qu[...,0:1]
        b = -_P*qu[...,3:4]
        c = -_P*qu[...,1:2]
        d = -_P*qu[...,2:3]

        eu = np.block([
            np.arctan2(b,a),
            np.arccos(2*(a**2+b**2)/(a**2+b**2+c**2+d**2)-1),
            np.arctan2(-d,c),
        ])

        eu_sum  = eu[...,0] + eu[...,2]
        eu_diff = eu[...,0] - eu[...,2]

        is_zero  = np.isclose(eu[...,1],0.0)
        is_pi    = np.isclose(eu[...,1],np.pi)
        is_ok    = ~np.logical_or(is_zero,is_pi)

        eu[...,0][is_zero] =  2*eu[...,0][is_zero]
        eu[...,0][is_pi]   = -2*eu[...,2][is_pi]
        eu[...,2][~is_ok]  = 0.0
        eu[...,0][is_ok]   = eu_diff[is_ok]
        eu[...,2][is_ok]   = eu_sum [is_ok]

        eu[np.logical_or(np.abs(eu)         < 1.e-6,
                         np.abs(eu-2*np.pi) < 1.e-6)] = 0.
        return np.where(eu < 0., eu%(np.pi*np.array([2.,1.,2.])),eu)

    @staticmethod
    def _qu2ax(qu: np.ndarray) -> np.ndarray:
        """
        Quaternion to axis–angle pair.

        Modified version of the original formulation, should be numerically more stable.
        """
        with np.errstate(invalid='ignore',divide='ignore'):
            s = np.sign(qu[...,0:1])/np.sqrt(qu[...,1:2]**2+qu[...,2:3]**2+qu[...,3:4]**2)
            omega = 2. * np.arccos(np.clip(qu[...,0:1],-1.,1.))
            ax = np.where(np.broadcast_to(qu[...,0:1] < 1.e-8,qu.shape),
                          np.block([qu[...,1:4],np.broadcast_to(np.pi,qu[...,0:1].shape)]),
                          np.block([qu[...,1:4]*s,omega]))
        ax[np.isclose(qu[...,0],1.,rtol=0.)] = np.array([0.,0.,1.,0.])
        return ax

    @staticmethod
    def _qu2ro(qu: np.ndarray) -> np.ndarray:
        """Quaternion to Rodrigues–Frank vector."""
        with np.errstate(invalid='ignore',divide='ignore'):
            s  = np.linalg.norm(qu[...,1:4],axis=-1,keepdims=True)
            ro = np.where(np.broadcast_to(np.abs(qu[...,0:1]) < 1.e-12,qu.shape),
                          np.block([qu[...,1:2], qu[...,2:3], qu[...,3:4], np.broadcast_to(np.inf,qu[...,0:1].shape)]),
                          np.block([qu[...,1:2]/s,qu[...,2:3]/s,qu[...,3:4]/s,
                                    np.tan(np.arccos(np.clip(qu[...,0:1],-1.,1.)))
                                   ])
                       )
        ro[np.abs(s).squeeze(-1) < 1.e-12] = np.array([0.,0.,_P,0.])
        return ro

    @staticmethod
    def _qu2ho(qu: np.ndarray) -> np.ndarray:
        """Quaternion to homochoric vector."""
        with np.errstate(invalid='ignore'):
            omega = 2. * np.arccos(np.clip(qu[...,0:1],-1.,1.))
            ho = np.where(np.abs(omega) < 1.e-12,
                          np.zeros(3),
                          qu[...,1:4]/np.linalg.norm(qu[...,1:4],axis=-1,keepdims=True)
                          * (0.75*(omega - np.sin(omega)))**(1./3.))
        return ho

    @staticmethod
    def _qu2cu(qu: np.ndarray) -> np.ndarray:
        """Quaternion to cubochoric vector."""
        return Rotation._ho2cu(Rotation._qu2ho(qu))


    #---------- Rotation matrix ----------
    @staticmethod
    def _om2qu(om: np.ndarray) -> np.ndarray:
        """
        Rotation matrix to quaternion.

        This formulation is from  www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion.
        The original formulation had issues.
        """
        trace = om[...,0,0:1] + om[...,1,1:2] + om[...,2,2:3]

        with np.errstate(invalid='ignore',divide='ignore'):
            s = np.array([
                 0.5 / np.sqrt( 1. + trace),
                 2.  * np.sqrt( 1. + om[...,0,0:1] - om[...,1,1:2] - om[...,2,2:3]),
                 2.  * np.sqrt( 1. + om[...,1,1:2] - om[...,2,2:3] - om[...,0,0:1]),
                 2.  * np.sqrt( 1. + om[...,2,2:3] - om[...,0,0:1] - om[...,1,1:2] )
                ])
            qu = np.where(trace>0,
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
                         )*np.array([1.,_P,_P,_P])
            qu[qu[...,0] < 0.] *= -1.
        return qu

    @staticmethod
    def _om2eu(om: np.ndarray) -> np.ndarray:
        """Rotation matrix to Bunge Euler angles."""
        with np.errstate(invalid='ignore',divide='ignore'):
            zeta = 1./np.sqrt(1.-om[...,2,2:3]**2)
            eu = np.where(np.isclose(np.abs(om[...,2,2:3]),1.,0.),
                          np.block([np.arctan2(om[...,0,1:2],om[...,0,0:1]),
                                    np.pi*0.5*(1.-om[...,2,2:3]),
                                    np.zeros(om.shape[:-2]+(1,)),
                                   ]),
                          np.block([np.arctan2(om[...,2,0:1]*zeta,-om[...,2,1:2]*zeta),
                                    np.arccos( om[...,2,2:3]),
                                    np.arctan2(om[...,0,2:3]*zeta,+om[...,1,2:3]*zeta)
                                   ])
                          )
        eu[np.abs(eu) < 1.e-8] = 0.0
        return np.where(eu < 0., eu%(np.pi*np.array([2.,1.,2.])),eu)

    @staticmethod
    def _om2ax(om: np.ndarray) -> np.ndarray:
        """Rotation matrix to axis–angle pair."""
        diag_delta = -_P*np.block([om[...,1,2:3]-om[...,2,1:2],
                                   om[...,2,0:1]-om[...,0,2:3],
                                   om[...,0,1:2]-om[...,1,0:1]
                                 ])
        t = 0.5*(om.trace(axis2=-2,axis1=-1) -1.).reshape(om.shape[:-2]+(1,))
        w,vr = np.linalg.eig(om)
        # mask duplicated real eigenvalues
        w[np.isclose(w[...,0],1.+0.j),1:] = 0.
        w[np.isclose(w[...,1],1.+0.j),2:] = 0.
        vr = np.swapaxes(vr,-1,-2)
        ax = np.where(np.abs(diag_delta)<1.e-13,
                             np.real(vr[np.isclose(w,1.+0.j)]).reshape(om.shape[:-2]+(3,)),
                      np.abs(np.real(vr[np.isclose(w,1.+0.j)]).reshape(om.shape[:-2]+(3,)))
                      *np.sign(diag_delta))
        ax = np.block([ax,np.arccos(np.clip(t,-1.,1.))])
        ax[np.abs(ax[...,3]) < 1.e-8] = np.array([0.,0.,1.,0.])
        return ax

    @staticmethod
    def _om2ro(om: np.ndarray) -> np.ndarray:
        """Rotation matrix to Rodrigues–Frank vector."""
        return Rotation._eu2ro(Rotation._om2eu(om))

    @staticmethod
    def _om2ho(om: np.ndarray) -> np.ndarray:
        """Rotation matrix to homochoric vector."""
        return Rotation._ax2ho(Rotation._om2ax(om))

    @staticmethod
    def _om2cu(om: np.ndarray) -> np.ndarray:
        """Rotation matrix to cubochoric vector."""
        return Rotation._ho2cu(Rotation._om2ho(om))


    #---------- Bunge Euler angles ----------
    @staticmethod
    def _eu2qu(eu: np.ndarray) -> np.ndarray:
        """Bunge Euler angles to quaternion."""
        ee = 0.5*eu
        cPhi = np.cos(ee[...,1:2])
        sPhi = np.sin(ee[...,1:2])
        qu = np.block([    cPhi*np.cos(ee[...,0:1]+ee[...,2:3]),
                       -_P*sPhi*np.cos(ee[...,0:1]-ee[...,2:3]),
                       -_P*sPhi*np.sin(ee[...,0:1]-ee[...,2:3]),
                       -_P*cPhi*np.sin(ee[...,0:1]+ee[...,2:3])])
        qu[qu[...,0] < 0.] *= -1.
        return qu

    @staticmethod
    def _eu2om(eu: np.ndarray) -> np.ndarray:
        """Bunge Euler angles to rotation matrix."""
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
        om[np.abs(om) < 1.e-12] = 0.
        return om

    @staticmethod
    def _eu2ax(eu: np.ndarray) -> np.ndarray:
        """Bunge Euler angles to axis–angle pair."""
        t = np.tan(eu[...,1:2]*0.5)
        sigma = 0.5*(eu[...,0:1]+eu[...,2:3])
        delta = 0.5*(eu[...,0:1]-eu[...,2:3])
        tau   = np.linalg.norm(np.block([t,np.sin(sigma)]),axis=-1,keepdims=True)
        alpha = np.where(np.abs(np.cos(sigma))<1.e-12,np.pi,2.*np.arctan(tau/np.cos(sigma)))
        with np.errstate(invalid='ignore',divide='ignore'):
            ax = np.where(np.broadcast_to(np.abs(alpha)<1.e-12,eu.shape[:-1]+(4,)),
                          [0.,0.,1.,0.],
                          np.block([-_P/tau*t*np.cos(delta),
                                    -_P/tau*t*np.sin(delta),
                                    -_P/tau*  np.sin(sigma),
                                     alpha
                                    ]))
        ax[(alpha<0.).squeeze()] *= -1.
        return ax

    @staticmethod
    def _eu2ro(eu: np.ndarray) -> np.ndarray:
        """Bunge Euler angles to Rodrigues–Frank vector."""
        ax = Rotation._eu2ax(eu)
        ro = np.block([ax[...,:3],np.tan(ax[...,3:4]*.5)])
        ro[ax[...,3] >= np.pi,3] = np.inf
        ro[np.abs(ax[...,3])<1.e-16] = np.array([0.,0.,_P,0.])
        return ro

    @staticmethod
    def _eu2ho(eu: np.ndarray) -> np.ndarray:
        """Bunge Euler angles to homochoric vector."""
        return Rotation._ax2ho(Rotation._eu2ax(eu))

    @staticmethod
    def _eu2cu(eu: np.ndarray) -> np.ndarray:
        """Bunge Euler angles to cubochoric vector."""
        return Rotation._ho2cu(Rotation._eu2ho(eu))


    #---------- Axis angle pair ----------
    @staticmethod
    def _ax2qu(ax: np.ndarray) -> np.ndarray:
        """Axis–angle pair to quaternion."""
        c = np.cos(ax[...,3:4]*.5)
        s = np.sin(ax[...,3:4]*.5)
        qu = np.where(np.abs(ax[...,3:4]) < 1.e-6,[1.,0.,0.,0.],np.block([c,ax[...,:3]*s]))
        return qu

    @staticmethod
    def _ax2om(ax: np.ndarray) -> np.ndarray:
        """Axis-angle pair to rotation matrix."""
        c = np.cos(ax[...,3:4])
        s = np.sin(ax[...,3:4])
        omc = 1.-c
        om = np.block([c+omc*ax[...,0:1]**2,
                         omc*ax[...,0:1]*ax[...,1:2] + s*ax[...,2:3],
                         omc*ax[...,0:1]*ax[...,2:3] - s*ax[...,1:2],
                         omc*ax[...,0:1]*ax[...,1:2] - s*ax[...,2:3],
                       c+omc*ax[...,1:2]**2,
                         omc*ax[...,1:2]*ax[...,2:3] + s*ax[...,0:1],
                         omc*ax[...,0:1]*ax[...,2:3] + s*ax[...,1:2],
                         omc*ax[...,1:2]*ax[...,2:3] - s*ax[...,0:1],
                       c+omc*ax[...,2:3]**2]).reshape(ax.shape[:-1]+(3,3))
        return om if _P < 0. else np.swapaxes(om,-1,-2)

    @staticmethod
    def _ax2eu(ax: np.ndarray) -> np.ndarray:
        """Rotation matrix to Bunge Euler angles."""
        return Rotation._om2eu(Rotation._ax2om(ax))

    @staticmethod
    def _ax2ro(ax: np.ndarray) -> np.ndarray:
        """Axis–angle pair to Rodrigues–Frank vector."""
        ro = np.block([ax[...,:3],
                       np.where(np.isclose(ax[...,3:4],np.pi,atol=1.e-15,rtol=.0),
                                np.inf,
                                np.tan(ax[...,3:4]*0.5))
                      ])
        ro[np.abs(ax[...,3]) < 1.e-6] = np.array([.0,.0,_P,.0])
        return ro

    @staticmethod
    def _ax2ho(ax: np.ndarray) -> np.ndarray:
        """Axis–angle pair to homochoric vector."""
        f = (0.75 * ( ax[...,3:4] - np.sin(ax[...,3:4]) ))**(1./3.)
        return ax[...,:3] * f

    @staticmethod
    def _ax2cu(ax: np.ndarray) -> np.ndarray:
        """Axis–angle pair to cubochoric vector."""
        return Rotation._ho2cu(Rotation._ax2ho(ax))


    #---------- Rodrigues-Frank vector ----------
    @staticmethod
    def _ro2qu(ro: np.ndarray) -> np.ndarray:
        """Rodrigues–Frank vector to quaternion."""
        return Rotation._ax2qu(Rotation._ro2ax(ro))

    @staticmethod
    def _ro2om(ro: np.ndarray) -> np.ndarray:
        """Rodgrigues–Frank vector to rotation matrix."""
        return Rotation._ax2om(Rotation._ro2ax(ro))

    @staticmethod
    def _ro2eu(ro: np.ndarray) -> np.ndarray:
        """Rodrigues–Frank vector to Bunge Euler angles."""
        return Rotation._om2eu(Rotation._ro2om(ro))

    @staticmethod
    def _ro2ax(ro: np.ndarray) -> np.ndarray:
        """Rodrigues–Frank vector to axis–angle pair."""
        with np.errstate(invalid='ignore',divide='ignore'):
            ax = np.where(np.isfinite(ro[...,3:4]),
                 np.block([ro[...,0:3]*np.linalg.norm(ro[...,0:3],axis=-1,keepdims=True),2.*np.arctan(ro[...,3:4])]),
                 np.block([ro[...,0:3],np.broadcast_to(np.pi,ro[...,3:4].shape)]))
        ax[np.abs(ro[...,3]) < 1.e-8]  = np.array([0.,0.,1.,0.])
        return ax

    @staticmethod
    def _ro2ho(ro: np.ndarray) -> np.ndarray:
        """Rodrigues–Frank vector to homochoric vector."""
        f = np.where(np.isfinite(ro[...,3:4]),2.*np.arctan(ro[...,3:4]) -np.sin(2.*np.arctan(ro[...,3:4])),np.pi)
        return np.where(np.broadcast_to(np.sum(ro[...,0:3]**2,axis=-1,keepdims=True) < 1.e-8,ro[...,0:3].shape),
                        np.zeros(3), ro[...,0:3]* (0.75*f)**(1./3.))

    @staticmethod
    def _ro2cu(ro: np.ndarray) -> np.ndarray:
        """Rodrigues–Frank vector to cubochoric vector."""
        return Rotation._ho2cu(Rotation._ro2ho(ro))


    #---------- Homochoric vector----------
    @staticmethod
    def _ho2qu(ho: np.ndarray) -> np.ndarray:
        """Homochoric vector to quaternion."""
        return Rotation._ax2qu(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2om(ho: np.ndarray) -> np.ndarray:
        """Homochoric vector to rotation matrix."""
        return Rotation._ax2om(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2eu(ho: np.ndarray) -> np.ndarray:
        """Homochoric vector to Bunge Euler angles."""
        return Rotation._ax2eu(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2ax(ho: np.ndarray) -> np.ndarray:
        """Homochoric vector to axis–angle pair."""
        tfit = np.array([+0.9999999999999968,     -0.49999999999986866,     -0.025000000000632055,
                         -0.003928571496460683,   -0.0008164666077062752,   -0.00019411896443261646,
                         -0.00004985822229871769, -0.000014164962366386031, -1.9000248160936107e-6,
                         -5.72184549898506e-6,    +7.772149920658778e-6,    -0.00001053483452909705,
                         +9.528014229335313e-6,   -5.660288876265125e-6,    +1.2844901692764126e-6,
                         +1.1255185726258763e-6,  -1.3834391419956455e-6,   +7.513691751164847e-7,
                         -2.401996891720091e-7,   +4.386887017466388e-8,    -3.5917775353564864e-9])
        hmag_squared = np.sum(ho**2,axis=-1,keepdims=True)
        s = np.sum(tfit*hmag_squared**np.arange(len(tfit)),axis=-1,keepdims=True)
        with np.errstate(invalid='ignore'):
            return np.where(np.broadcast_to(np.abs(hmag_squared)<1.e-8,ho.shape[:-1]+(4,)),
                            [0.,0.,1.,0.],
                            np.block([ho/np.sqrt(hmag_squared),2.*np.arccos(np.clip(s,-1.,1.))]))

    @staticmethod
    def _ho2ro(ho: np.ndarray) -> np.ndarray:
        """Axis–angle pair to Rodrigues–Frank vector."""
        return Rotation._ax2ro(Rotation._ho2ax(ho))

    @staticmethod
    def _ho2cu(ho: np.ndarray) -> np.ndarray:
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
            xyz2 = xyz3[...,0:2] * np.sqrt( 2.*rs/(rs+np.abs(xyz3[...,2:3])) )
            qxy = np.sum(xyz2**2,axis=-1,keepdims=True)

            q2 = qxy + np.max(np.abs(xyz2),axis=-1,keepdims=True)**2
            sq2 = np.sqrt(q2)
            q = (_beta/np.sqrt(2.)/_R1) * np.sqrt(q2*qxy/(q2-np.max(np.abs(xyz2),axis=-1,keepdims=True)*sq2))
            tt = np.clip((np.min(np.abs(xyz2),axis=-1,keepdims=True)**2\
                +np.max(np.abs(xyz2),axis=-1,keepdims=True)*sq2)/np.sqrt(2.)/qxy,-1.,1.)
            T_inv = np.where(np.abs(xyz2[...,1:2]) <= np.abs(xyz2[...,0:1]),
                                np.block([np.ones_like(tt),np.arccos(tt)/np.pi*12.]),
                                np.block([np.arccos(tt)/np.pi*12.,np.ones_like(tt)]))*q
            T_inv[xyz2<0.] *= -1.
            T_inv[np.broadcast_to(np.isclose(qxy,0.,rtol=0.,atol=1.e-12),T_inv.shape)] = 0.
            cu = np.block([T_inv, np.where(xyz3[...,2:3]<0.,-np.ones_like(xyz3[...,2:3]),np.ones_like(xyz3[...,2:3])) \
                                  * rs/np.sqrt(6./np.pi),
                          ])/ _sc

        cu[np.isclose(np.sum(np.abs(ho),axis=-1),0.,rtol=0.,atol=1.e-16)] = 0.
        return np.take_along_axis(cu,Rotation._get_pyramid_order(ho,'backward'),-1)

    #---------- Cubochoric ----------
    @staticmethod
    def _cu2qu(cu: np.ndarray) -> np.ndarray:
        """Cubochoric vector to quaternion."""
        return Rotation._ho2qu(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2om(cu: np.ndarray) -> np.ndarray:
        """Cubochoric vector to rotation matrix."""
        return Rotation._ho2om(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2eu(cu: np.ndarray) -> np.ndarray:
        """Cubochoric vector to Bunge Euler angles."""
        return Rotation._ho2eu(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2ax(cu: np.ndarray) -> np.ndarray:
        """Cubochoric vector to axis–angle pair."""
        return Rotation._ho2ax(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2ro(cu: np.ndarray) -> np.ndarray:
        """Cubochoric vector to Rodrigues–Frank vector."""
        return Rotation._ho2ro(Rotation._cu2ho(cu))

    @staticmethod
    def _cu2ho(cu: np.ndarray) -> np.ndarray:
        """
        Cubochoric vector to homochoric vector.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        with np.errstate(invalid='ignore',divide='ignore'):
            # get pyramid and scale by grid parameter ratio
            XYZ = np.take_along_axis(cu,Rotation._get_pyramid_order(cu,'forward'),-1) * _sc
            order = np.abs(XYZ[...,1:2]) <= np.abs(XYZ[...,0:1])
            q = np.pi/12. * np.where(order,XYZ[...,1:2],XYZ[...,0:1]) \
                           / np.where(order,XYZ[...,0:1],XYZ[...,1:2])
            c = np.cos(q)
            s = np.sin(q)
            q = _R1*2.**0.25/_beta/ np.sqrt(np.sqrt(2.)-c) \
              * np.where(order,XYZ[...,0:1],XYZ[...,1:2])

            T = np.block([(np.sqrt(2.)*c - 1.), np.sqrt(2.) * s]) * q

            # transform to sphere grid (inverse Lambert)
            c = np.sum(T**2,axis=-1,keepdims=True)
            s = c *         np.pi/24. /XYZ[...,2:3]**2
            c = c * np.sqrt(np.pi/24.)/XYZ[...,2:3]
            q = np.sqrt( 1. - s)

            ho = np.where(np.isclose(np.sum(np.abs(XYZ[...,0:2]),axis=-1,keepdims=True),0.,rtol=0.,atol=1.e-16),
                          np.block([np.zeros_like(XYZ[...,0:2]),np.sqrt(6./np.pi)*XYZ[...,2:3]]),
                          np.block([np.where(order,T[...,0:1],T[...,1:2])*q,
                                    np.where(order,T[...,1:2],T[...,0:1])*q,
                                    np.sqrt(6./np.pi) * XYZ[...,2:3] - c])
                          )

        ho[np.isclose(np.sum(np.abs(cu),axis=-1),0.,rtol=0.,atol=1.e-16)] = 0.
        return np.take_along_axis(ho,Rotation._get_pyramid_order(cu,'backward'),-1)


    @staticmethod
    def _get_pyramid_order(xyz: np.ndarray,
                           direction: Literal['forward', 'backward']) -> np.ndarray:
        """
        Get order of the coordinates.

        Depending on the pyramid in which the point is located, the order need to be adjusted.

        Parameters
        ----------
        xyz : numpy.ndarray
           Coordinates of a point on a uniform refinable grid on a ball or
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
