import copy
import warnings
from typing import Optional, Union, TypeVar, Literal, Sequence, NamedTuple, overload

import numpy as np
import numpy.typing as npt

from ._typehints import FloatSequence, IntSequence, CrystalFamily, BravaisLattice, NumpyRngSeed
from . import Rotation
from . import Crystal
from . import util
from . import tensor


MyType = TypeVar('MyType', bound='Orientation')

class DisorientationTuple(NamedTuple):
    disorientation: 'Orientation'
    operators: np.ndarray


class AverageTuple(NamedTuple):
    average: 'Orientation'
    cloud: 'Orientation'


class ToSSTTuple(NamedTuple):
    vector_sst: np.ndarray
    operator: np.ndarray


class Orientation(Rotation,Crystal):
    """
    Representation of crystallographic orientation as combination of rotation and either crystal family or Bravais lattice.

    The crystal family is one of:

    - triclinic
    - monoclinic
    - orthorhombic
    - tetragonal
    - hexagonal
    - cubic

    and enables symmetry-related operations such as
    :func:`equivalent`, :func:`reduced` :func:`disorientation`,
    :func:`IPF_color`, or :func:`to_SST`.

    The Bravais lattice is given in the Pearson notation:

    - triclinic
       - aP : primitive

    - monoclinic
       - mP : primitive
       - mS : base-centered

    - orthorhombic
       - oP : primitive
       - oS : base-centered
       - oI : body-centered
       - oF : face-centered

    - tetragonal
       - tP : primitive
       - tI : body-centered

    - hexagonal
       - hP : primitive

    - cubic
       - cP : primitive
       - cI : body-centered
       - cF : face-centered

    and inherits the corresponding crystal family.
    Specifying a Bravais lattice, compared to just the crystal family,
    extends the functionality of Orientation objects to include operations
    such as :func:`Schmid`, :func:`related`, or :func:`to_frame` that
    require a lattice type and its parameters.

    Examples
    --------
    An array of 3 x 5 random orientations reduced to the fundamental zone of tetragonal symmetry:

    >>> import damask
    >>> o=damask.Orientation.from_random(shape=(3,5),family='tetragonal').reduced
    """

    def __init__(self,
                 rotation: Union[FloatSequence, Rotation] = np.array([1.,0.,0.,0.]),
                 *,
                 family: Optional[CrystalFamily] = None,
                 lattice: Optional[BravaisLattice] = None,
                 a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                 alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                 degrees: bool = False):
        """
        New orientation.

        Parameters
        ----------
        rotation : list, numpy.ndarray, or Rotation, optional
            Unit quaternion in positive real hemisphere.
            Use .from_quaternion to perform a sanity check.
            Defaults to no rotation.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.
        """
        Rotation.__init__(self,rotation)
        Crystal.__init__(self,
                         family=family, lattice=lattice,
                         a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                         degrees=degrees)


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give unambiguous representation.
        """
        return util.srepr([Crystal.__repr__(self),
                           Rotation.__repr__(self)])


    def __copy__(self: MyType,
                 rotation: Optional[Union[FloatSequence, Rotation]] = None) -> MyType:
        """
        Return deepcopy(self).

        Create deep copy.
        """
        dup = copy.deepcopy(self)
        if rotation is not None:
            dup.quaternion = Rotation(rotation).quaternion
        return dup

    copy = __copy__


    def __eq__(self,                                                                                # type: ignore[override]
               other: object) -> npt.NDArray[np.bool_]:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : Orientation
            Orientation to check for equality.

        Returns
        -------
        equal : numpy.ndarray
            Whether both arguments are equal.
        """
        if not isinstance(other, Orientation):
            return NotImplemented
        matching_type = self.family == other.family and \
                        self.lattice == other.lattice and \
                        self.parameters == other.parameters
        return np.logical_and(matching_type,super(self.__class__,self.reduced).__eq__(other.reduced))

    def __ne__(self,                                                                                # type: ignore[override]
               other: object) -> npt.NDArray[np.bool_]:
        """
        Return self!=other.

        Test inequality of other.

        Parameters
        ----------
        other : Orientation
            Orientation to check for equality.
        """
        return np.logical_not(self==other) if isinstance(other, Orientation) else NotImplemented


    def isclose(self: MyType,
                other: MyType,
                rtol: float = 1.e-5,
                atol: float = 1.e-8,
                equal_nan: bool = True) -> npt.NDArray[np.bool_]:
        """
        Report where values are approximately equal to corresponding ones of other Orientation.

        Parameters
        ----------
        other : Orientation
            Orientation to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        mask : numpy.ndarray of bool, shape (self.shape)
            Mask indicating where corresponding orientations are close.
        """
        matching_type = self.family == other.family and \
                        self.lattice == other.lattice and \
                        self.parameters == other.parameters
        return np.logical_and(matching_type,super(self.__class__,self.reduced).isclose(other.reduced))

    def allclose(self: MyType,
                 other: MyType,
                 rtol: float = 1.e-5,
                 atol: float = 1.e-8,
                 equal_nan: bool = True) -> np.bool_:
        """
        Test whether all values are approximately equal to corresponding ones of other Orientation.

        Parameters
        ----------
        other : Orientation
            Orientation to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        allclose : bool
            Test whether all values are approximately equal to corresponding
            ones of other orientation.
        """
        return np.all(self.isclose(other,rtol,atol,equal_nan))


    def __mul__(self: MyType,
                other: Union[Rotation, 'Orientation']) -> MyType:
        """
        Return self*other.

        Compose with other.

        Parameters
        ----------
        other : Rotation or Orientation, shape (self.shape)
            Object for composition.
            Compatible innermost dimensions will blend.

        Returns
        -------
        composition : Orientation
            Compound rotation self*other, i.e. first other then self rotation.
        """
        if not isinstance(other, (Orientation,Rotation)):
            raise TypeError('use "O@b", i.e. matmul, to apply Orientation "O" to object "b"')
        return self.copy(Rotation(self.quaternion)*Rotation(other.quaternion))


    @staticmethod
    def from_quaternion(q: Union[Sequence[FloatSequence], np.ndarray],
                        accept_homomorph: bool = False,
                        normalize: bool = False,
                        P: Literal[1, -1] = -1,
                        *,
                        family: Optional[CrystalFamily] = None,
                        lattice: Optional[BravaisLattice] = None,
                        a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                        alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                        degrees: bool = False) -> 'Orientation':
        """
        Initialize from quaternion.

        Parameters
        ----------
        q : numpy.ndarray, shape (...,4)
            Unit quaternion (q_0, q_1, q_2, q_3) in positive real hemisphere, i.e. ǀqǀ = 1 and q_0 ≥ 0.
        accept_homomorph : bool, optional
            Allow homomorphic variants, i.e. q_0 < 0 (negative real hemisphere).
            Defaults to False.
        normalize : bool, optional
            Allow ǀqǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given quaternion.
        """
        return Orientation(Rotation.from_quaternion(q,accept_homomorph,normalize,P),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_Euler_angles(phi: np.ndarray,
                          degrees: bool = False,
                          *,
                          family: Optional[CrystalFamily] = None,
                          lattice: Optional[BravaisLattice] = None,
                          a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                          alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                          ) -> 'Orientation':
        """
        Initialize from Bunge Euler angles.

        Parameters
        ----------
        phi : numpy.ndarray, shape (...,3)
            Euler angles (φ_1 ∈ [0,2π], ϕ ∈ [0,π], φ_2 ∈ [0,2π])
            or (φ_1 ∈ [0,360], ϕ ∈ [0,180], φ_2 ∈ [0,360]) if degrees == True.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given Bunge Euler angles.

        Notes
        -----
        Bunge Euler angles correspond to a rotation axis sequence of z–x'–z''.
        """
        return Orientation(Rotation.from_Euler_angles(phi,degrees),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_axis_angle(n_omega: np.ndarray,
                        degrees: bool = False,
                        normalize: bool = False,
                        P: Literal[1, -1] = -1,
                        *,
                        family: Optional[CrystalFamily] = None,
                        lattice: Optional[BravaisLattice] = None,
                        a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                        alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                        ) -> 'Orientation':
        """
        Initialize from axis–angle pair.

        Parameters
        ----------
        n_omega : numpy.ndarray, shape (...,4)
            Axis and angle (n_1, n_2, n_3, ω) with ǀnǀ = 1 and ω ∈ [0,π]
            or ω ∈ [0,180] if degrees == True.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.
        normalize : bool, optional
            Allow ǀnǀ ≠ 1. Defaults to False.
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given axis-angle pair.
        """
        return Orientation(Rotation.from_axis_angle(n_omega,degrees,normalize,P),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_basis(basis: np.ndarray,
                   orthonormal: bool = True,
                   reciprocal: bool = False,
                   *,
                   family: Optional[CrystalFamily] = None,
                   lattice: Optional[BravaisLattice] = None,
                   a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                   alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                   degrees: bool = False) -> 'Orientation':
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
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given basis.
        """
        return Orientation(Rotation.from_basis(basis,orthonormal,reciprocal),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_matrix(R: np.ndarray,
                    normalize: bool = False,
                    *,
                    family: Optional[CrystalFamily] = None,
                    lattice: Optional[BravaisLattice] = None,
                    a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                    alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                    degrees: bool = False) -> 'Orientation':
        """
        Initialize from rotation matrix.

        Parameters
        ----------
        R : numpy.ndarray, shape (...,3,3)
            Rotation matrix with det(R) = 1 and R.T ∙ R = I.
        normalize : bool, optional
            Rescales rotation matrix to unit determinant. Defaults to False.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given rotation matrix.
        """
        return Orientation(Rotation.from_matrix(R,normalize),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_parallel(source: np.ndarray,
                      target: np.ndarray,
                      active: bool = False,
                      *,
                      family: Optional[CrystalFamily] = None,
                      lattice: Optional[BravaisLattice] = None,
                      a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                      alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                      degrees: bool = False) -> 'Orientation':
        """
        Initialize from pairs of two orthogonal basis vectors.

        Parameters
        ----------
        source : numpy.ndarray, shape (...,2,3)
            First and second three-dimensional vector of orthogonal source basis.
        target : numpy.ndarray, shape (...,2,3)
            Corresponding three-dimensional vectors of target basis.
        active : bool, optional
            Consider rotations as active, i.e. return (B^-1⋅A) instead of (B⋅A^-1).
            Defaults to False.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given basis vectors.

        Notes
        -----
        If rotations $A = [s_1,s_2,s_1 × s_2]^T$ and B = $[t_1,t_2,t_1 × t_2]^T$
        are considered "active", the resulting rotation will be $B^{-1}⋅A$ instead
        of the default result $B⋅A^{-1}$.
        """
        return Orientation(Rotation.from_parallel(source,target,active),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_Rodrigues_vector(rho: np.ndarray,
                              normalize: bool = False,
                              P: Literal[1, -1] = -1,
                              *,
                              family: Optional[CrystalFamily] = None,
                              lattice: Optional[BravaisLattice] = None,
                              a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                              alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                              degrees: bool = False) -> 'Orientation':
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
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given Rodrigues–Frank vector.
        """
        return Orientation(Rotation.from_Rodrigues_vector(rho,normalize,P),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_homochoric(h: np.ndarray,
                        P: Literal[1, -1] = -1,
                        *,
                        family: Optional[CrystalFamily] = None,
                        lattice: Optional[BravaisLattice] = None,
                        a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                        alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                        degrees: bool = False) -> 'Orientation':
        """
        Initialize from homochoric vector.

        Parameters
        ----------
        h : numpy.ndarray, shape (...,3)
            Homochoric vector (h_1, h_2, h_3) with ǀhǀ < (3π/4)^(1/3).
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given homochoric vector.
        """
        return Orientation(Rotation.from_homochoric(h,P),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_cubochoric(x: np.ndarray,
                        P: Literal[1, -1] = -1,
                        *,
                        family: Optional[CrystalFamily] = None,
                        lattice: Optional[BravaisLattice] = None,
                        a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                        alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                        degrees: bool = False) -> 'Orientation':
        """
        Initialize from cubochoric vector.

        Parameters
        ----------
        x : numpy.ndarray, shape (...,3)
            Cubochoric vector (x_1, x_2, x_3) with max(x_i) < 1/2 π^(2/3).
        P : int ∈ {-1,1}, optional
            Sign convention. Defaults to -1.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the given cubochoric vector.
        """
        return Orientation(Rotation.from_cubochoric(x,P),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_random(shape: Union[None, int, IntSequence] = None,
                    rng_seed: Optional[NumpyRngSeed] = None,
                    *,
                    family: Optional[CrystalFamily] = None,
                    lattice: Optional[BravaisLattice] = None,
                    a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                    alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                    degrees: bool = False) -> 'Orientation':
        """
        Initialize with samples from a uniform distribution.

        Parameters
        ----------
        shape : (sequence of) int, optional
            Output shape. Defaults to None, which gives a scalar.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Random orientation of given shape.
        """
        return Orientation(Rotation.from_random(shape,rng_seed),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_ODF(weights: np.ndarray,
                 phi: np.ndarray,
                 shape: Union[None, int, IntSequence] = None,
                 degrees: bool = False,
                 fractions: bool = True,
                 rng_seed: Optional[NumpyRngSeed] = None,
                 *,
                 family: Optional[CrystalFamily] = None,
                 lattice: Optional[BravaisLattice] = None,
                 a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                 alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                 ) -> 'Orientation':
        """
        Initialize with samples from a binned orientation distribution function (ODF).

        Parameters
        ----------
        weights : numpy.ndarray, shape (n)
            Texture intensity values (probability density or volume fraction) at Euler space grid points.
        phi : numpy.ndarray, shape (n,3)
            Grid coordinates in Euler space at which weights are defined.
        shape : (sequence of) int, optional
            Output shape. Defaults to None, which gives a scalar.
        degrees : bool, optional
            Angles are given in degrees. Defaults to True.
        fractions : bool, optional
            ODF values correspond to volume fractions, not probability densities.
            Defaults to True.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.

        Returns
        -------
        new : damask.Orientation
            Orientation sampled from given ODF.

        Notes
        -----
        Due to the distortion of Euler space in the vicinity of ϕ = 0,
        probability densities, p, defined on grid points with ϕ = 0 will never
        result in reconstructed orientations as their dV/V = p dγ = p × 0.
        Hence, it is recommended to transform any such dataset to a
        cell-centered version, which avoids grid points at ϕ = 0.

        References
        ----------
        P. Eisenlohr and F. Roters, Computational Materials Science 42(4):670-678, 2008
        https://doi.org/10.1016/j.commatsci.2007.09.015
        """
        return Orientation(Rotation.from_ODF(weights,phi,shape,degrees,fractions,rng_seed),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_spherical_component(center: Union[Rotation,'Orientation'],
                                 sigma: float,
                                 shape: Union[None, int, IntSequence] = None,
                                 degrees: bool = False,
                                 rng_seed: Optional[NumpyRngSeed] = None,
                                 *,
                                 family: Optional[CrystalFamily] = None,
                                 lattice: Optional[BravaisLattice] = None,
                                 a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                                 alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                                 ) -> 'Orientation':
        """
        Initialize with samples from a Gaussian distribution around a given center.

        Parameters
        ----------
        center : Rotation or Orientation
            Central rotation.
        sigma : float
            Standard deviation of (Gaussian) misorientation distribution.
        shape : (sequence of) int, optional
            Output shape. Defaults to None, which gives a scalar.
        degrees : bool, optional
            Standard deviation is given in degrees. Defaults to False.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.

        Returns
        -------
        new : damask.Orientation
            Orientation sampled from normal distribution around a center.
        """
        return Orientation(Rotation.from_spherical_component(center,sigma,shape,degrees,rng_seed),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_fiber_component(crystal: IntSequence,
                             sample: IntSequence,
                             sigma: float = 0.,
                             shape: Union[None, int, IntSequence] = None,
                             degrees: bool = False,
                             rng_seed: Optional[NumpyRngSeed] = None,
                             *,
                             family: Optional[CrystalFamily] = None,
                             lattice: Optional[BravaisLattice] = None,
                             a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                             alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                             ) -> 'Orientation':
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
            Output shape. Defaults to None, which gives a scalar.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator.
            Defaults to None, i.e. unpredictable entropy will be pulled from the OS.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.

        Returns
        -------
        new : damask.Orientation
            Orientation sampled from normal distribution around a direction.

        Notes
        -----
        The crystal direction for (θ=0,φ=0) is [0 0 1];
        the sample direction for (θ=0,φ=0) is z.

        Polar coordinates follow the ISO 80000-2:2019 convention
        typically used in physics.
        See https://en.wikipedia.org/wiki/Spherical_coordinate_system.

        Ranges 0≤θ≤π and 0≤φ≤2π give a unique set of coordinates.

        The value for sigma should be selected small enough such that
        sampling from a normal distribution does not result in values
        larger than the maximum disorientation of the crystal family.
        The maximum disorientation angle d_max is 62.80°, 93.84°,
        98.43°, and 120° for cubic, hexagonal, tetragonal, and
        orthorhombic crystal families, respectively.
        A typically safe value is σ = d_max/5.

        References
        ----------
        A. Heinz and P. Neumann, Acta Crystallographica Section A 47:780-789, 1991
        https://doi.org/10.1107/S0108767391006864
        """
        return Orientation(Rotation.from_fiber_component(crystal,sample,sigma,shape,degrees,rng_seed),
                           family=family,lattice=lattice,
                           a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                           degrees=degrees)

    @staticmethod
    def from_directions(uvw: IntSequence,
                        hkl: IntSequence,
                        *,
                        family: Optional[CrystalFamily] = None,
                        lattice: Optional[BravaisLattice] = None,
                        a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                        alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                        degrees: bool = False) -> 'Orientation':
        """
        Initialize orientation object from the crystallographic direction and plane parallel to lab x and z, respectively.

        Parameters
        ----------
        uvw : numpy.ndarray, shape (...,3)
            Lattice direction as Miller indices aligned with lab frame x-direction.
        hkl : numpy.ndarray, shape (...,3)
            Lattice plane normal as Miller indices aligned with lab frame z-direction.
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

        Returns
        -------
        new : damask.Orientation
            Orientation representing the relationship between direction and plane
            and lab x- and z-direction.
        """
        o = Orientation(rotation=[1,0,0,0],
                        family=family,lattice=lattice,
                        a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma,
                        degrees=degrees)
        x = o.to_frame(uvw=uvw)
        z = o.to_frame(hkl=hkl)
        om = np.stack([x,np.cross(z,x),z],axis=-2)
        return o.copy(Rotation.from_matrix(tensor.transpose(om/np.linalg.norm(om,axis=-1,keepdims=True))))


    @property
    def equivalent(self: MyType) -> MyType:
        """
        Orientations that are symmetrically equivalent.

        One dimension (length corresponds to number of symmetrically equivalent orientations)
        is added to the left of the Rotation array.
        """
        sym_ops = self.symmetry_operations
        o = sym_ops.broadcast_to(sym_ops.shape+self.shape,mode='right')
        return self.copy(o*Rotation(self.quaternion).broadcast_to(o.shape,mode='left'))


    @property
    def reduced(self: MyType) -> MyType:
        """Select symmetrically equivalent orientation that falls into fundamental zone according to symmetry."""
        eq   = self.equivalent
        ok   = eq.in_FZ
        ok  &= np.cumsum(ok,axis=0) == 1
        loc  = np.where(ok)
        sort = 0 if len(loc) == 1 else np.lexsort(loc[:0:-1])
        return eq[ok][sort].reshape(self.shape)


    @property
    def in_FZ(self) -> Union[np.bool_, np.ndarray]:
        """
        Check whether orientation falls into fundamental zone of own symmetry.

        Returns
        -------
        in : numpy.ndarray of bool, shape (self.shape)
            Whether Rodrigues-Frank vector falls into fundamental zone.

        Notes
        -----
        Fundamental zones in Rodrigues space are point-symmetric around origin.

        References
        ----------
        A. Heinz and P. Neumann, Acta Crystallographica Section A 47:780-789, 1991
        https://doi.org/10.1107/S0108767391006864
        """
        rho_abs = np.abs(self.as_Rodrigues_vector(compact=True))*(1.-1.e-9)

        with np.errstate(invalid='ignore'):
            # using '*'/prod for 'and'
            if   self.family == 'cubic':
                return (np.prod(np.sqrt(2)-1. >= rho_abs,axis=-1) *
                                   (1. >= np.sum(rho_abs,axis=-1))).astype(bool)
            if self.family == 'hexagonal':
                return (np.prod(1.  >= rho_abs,axis=-1) *
                                (2. >= np.sqrt(3)*rho_abs[...,0] + rho_abs[...,1]) *
                                (2. >= np.sqrt(3)*rho_abs[...,1] + rho_abs[...,0]) *
                                (2. >= np.sqrt(3)                + rho_abs[...,2])).astype(bool)
            if self.family == 'tetragonal':
                return (np.prod(1.  >= rho_abs[...,:2],axis=-1) *
                        (np.sqrt(2) >= rho_abs[...,0] + rho_abs[...,1]) *
                        (np.sqrt(2) >= rho_abs[...,2] + 1.)).astype(bool)
            if self.family == 'orthorhombic':
                return (np.prod(1. >= rho_abs,axis=-1)).astype(bool)
            if self.family == 'monoclinic':
                return np.logical_or(   1. >= rho_abs[...,1],
                                     np.isnan(rho_abs[...,1]))
            if self.family == 'triclinic':
                return np.ones(rho_abs.shape[:-1]).astype(bool)

            raise TypeError(f'unknown symmetry "{self.family}"')


    @property
    def in_disorientation_FZ(self) -> np.ndarray:
        """
        Check whether orientation falls into fundamental zone of disorientations.

        Returns
        -------
        in : numpy.ndarray of bool, shape (self.shape)
            Whether Rodrigues-Frank vector falls into disorientation FZ.

        References
        ----------
        A. Heinz and P. Neumann, Acta Crystallographica Section A 47:780-789, 1991
        https://doi.org/10.1107/S0108767391006864
        """
        def larger_or_equal(v,c):
            return ((np.isclose(c[0],v[...,0]) | (v[...,0] > c[0])) &
                    (np.isclose(c[1],v[...,1]) | (v[...,1] > c[1])) &
                    (np.isclose(c[2],v[...,2]) | (v[...,2] > c[2]))).astype(bool)

        rho = self.as_Rodrigues_vector(compact=True)
        return larger_or_equal(rho,
                                     [rho[...,1],           rho[...,2],0] if self.family == 'cubic'
                                else [rho[...,1]*np.sqrt(3),0,         0] if self.family == 'hexagonal'
                                else [rho[...,1],           0,         0] if self.family == 'tetragonal'
                                else [0,                    0,         0] if self.family == 'orthorhombic'
                                else [-np.inf,              0,         0] if self.family == 'monoclinic'
                                else [-np.inf,        -np.inf,   -np.inf]) & self.in_FZ


    @overload
    def disorientation(self: MyType, other: MyType,
                       return_operators: Literal[False] = False) -> MyType:
        ...
    @overload
    def disorientation(self: MyType, other: MyType,
                       return_operators: Literal[True] = True) -> DisorientationTuple:
        ...
    def disorientation(self: MyType,
                       other: MyType,
                       return_operators: bool = False) -> Union[MyType, DisorientationTuple]:
        """
        Calculate disorientation between self and given other orientation.

        Parameters
        ----------
        other : Orientation
            Orientation to which the disorientation is computed.
            Compatible innermost dimensions will blend.
        return_operators : bool, optional
            Return index pair of symmetrically equivalent orientations
            that result in misorientation falling into disorientation FZ.
            Defaults to False.

        Returns
        -------
        disorientation : Orientation
            Disorientation between self and other.
        operators : numpy.ndarray of int, shape (...,2), conditional
            Index pair of symmetrically equivalent orientations
            that result in misorientation falling into disorientation FZ.

        Notes
        -----
        Requires same crystal family for both orientations.

        Examples
        --------
        Disorientation between two specific orientations of hexagonal symmetry:

        >>> import damask
        >>> a = damask.Orientation.from_Euler_angles(phi=[123,32,21],degrees=True,family='hexagonal')
        >>> b = damask.Orientation.from_Euler_angles(phi=[104,11,87],degrees=True,family='hexagonal')
        >>> a.disorientation(b)
        Crystal family: hexagonal
        array((0.976478  ,     0.18880082,  0.01784483,  0.10259889))

        Plot a sample from the Mackenzie distribution.

        >>> import matplotlib.pyplot as plt
        >>> import damask
        >>> N = 10000
        >>> a = damask.Orientation.from_random(shape=N,family='cubic')
        >>> b = damask.Orientation.from_random(shape=N,family='cubic')
        >>> n,omega = a.disorientation(b).as_axis_angle(degrees=True,pair=True)
        >>> plt.hist(omega,25)
        (...)
        >>> plt.show(block=False)
        """
        # For extension to cases with differing symmetry see
        # https://doi.org/10.1107/S0021889808016373 and https://doi.org/10.1107/S0108767391006864
        if self.family != other.family:
            raise NotImplementedError('disorientation between different crystal families')

        blend = util.shapeblender( self.shape,other.shape)
        s_m   = util.shapeshifter( self.shape,blend,mode='right')
        s_o   = util.shapeshifter(other.shape,blend,mode='left')

        s =  self.broadcast_to(s_m).equivalent
        o = other.broadcast_to(s_o).equivalent

        r_ = s[:,np.newaxis,...].misorientation(o[np.newaxis,:,...])
        _r = ~r_
        shp = r_.shape[2:]

        forward = r_.in_disorientation_FZ
        reverse = _r.in_disorientation_FZ
        ok  = forward | reverse
        ok &= (np.cumsum(ok.reshape((-1,*shp)),axis=0) == 1).reshape(ok.shape)
        r = np.where(np.any((ok&forward)[...,np.newaxis],axis=(0,1),keepdims=True),
                     r_.quaternion,
                     _r.quaternion)
        loc  = np.where(ok)
        sort = 0 if len(loc) == 2 else np.lexsort(loc[:1:-1])

        quat = r[ok][sort].reshape((*shp,4))

        if return_operators:
            operators = (np.vstack(loc[:2]).T)[sort].reshape((*shp, 2))
            return DisorientationTuple(self.copy(rotation=quat), operators)
        else:
            return self.copy(rotation=quat)


    def disorientation_angle(self: MyType,
                             other: MyType) -> np.ndarray:
        """
        Calculate disorientation angle between self and given other orientation.

        Parameters
        ----------
        other : Orientation
            Orientation to which the disorientation angle is computed.
            Compatible innermost dimensions will blend.

        Returns
        -------
        omega : np.ndarray
            Disorientation angle.

        Notes
        -----
        Requires same crystal family for both orientations.

        References
        ----------
        Lionel Germain, personal communication.
        """
        q_abs = np.abs((self*~other).quaternion)

        if   'triclinic' == other.family == self.family:
            trace_max = q_abs[...,0:1]

        elif 'monoclinic' == other.family == self.family:
            trace_max = np.maximum(q_abs[...,0:1],
                                   q_abs[...,2:3])

        elif 'orthorhombic' == other.family == self.family:
            trace_max = np.maximum.reduce([q_abs[...,0:1],
                                           q_abs[...,1:2],
                                           q_abs[...,2:3],
                                           q_abs[...,3:4]])

        elif 'tetragonal' == other.family == self.family:
            m1,m2,m3,m4 = np.split(q_abs,4,axis=-1)

            trace_max = np.maximum.reduce([m1,m2,m3,m4,
                                           (m1+m4)*np.sqrt(2.)/2.,
                                           (m2+m3)*np.sqrt(2.)/2.])

        elif 'hexagonal' == other.family == self.family:
            m1,m2,m3,m4 = np.split(q_abs,4,axis=-1)

            mask = m1 < m4
            m1[mask],m4[mask] = m4[mask],m1[mask]
            mask = m2 < m3
            m2[mask],m3[mask] = m3[mask],m2[mask]

            trace_max = np.maximum.reduce([m1,m2,
                                           m1*np.sqrt(3.)/2.+m4*.5,
                                           m2*np.sqrt(3.)/2.+m3*.5])

        elif 'cubic' == other.family == self.family:
            m1,m2,m3,m4 = np.split(q_abs,4,axis=-1)

            trace_max = np.sum(q_abs,axis=-1,keepdims=True)*.5

            mask = m1 < m2
            m1[mask],m2[mask] = m2[mask],m1[mask]
            mask = m3 < m4
            m3[mask],m4[mask] = m4[mask],m3[mask]

            mask1 = m1 > m3
            mask2 = np.logical_and(mask1,m2<m3)
            mask3 = np.logical_not(mask1)

            m2[mask2] = m3[mask2]
            m2[mask3] = np.where(m4[mask3]<m1[mask3],m1[mask3],m4[mask3])
            m1[mask3] = m3[mask3]

            trace_max = np.maximum.reduce([trace_max,m1,(m1+m2)*np.sqrt(2.)/2.])

        else:
            return self.disorientation(other).as_axis_angle(pair=True)[1]

        return 2.*np.arccos(np.clip(trace_max[...,0].round(15),None,1.))


    def average(self: MyType,                                                                       # type: ignore[override]
                weights: Optional[FloatSequence] = None,
                return_cloud: Optional[bool] = False) -> Union[MyType, AverageTuple]:
        """
        Return orientation average over last dimension.

        Parameters
        ----------
        weights : numpy.ndarray, shape (self.shape), optional
            Relative weights of orientations.
            Defaults to equal weights.
        return_cloud : bool, optional
            Return the specific (symmetrically equivalent) orientations that were averaged.
            Defaults to False.

        Returns
        -------
        average : Orientation
            Weighted average of original Orientation field.
        cloud : Orientations, conditional
            Symmetrically equivalent version of each orientation that were actually used in averaging.

        References
        ----------
        J.C. Glez and J. Driver, Journal of Applied Crystallography 34:280-288, 2001
        https://doi.org/10.1107/S0021889801003077
        """
        eq = self.equivalent
        m  = eq.misorientation_angle(self[...,0].reshape((1,)+self.shape[:-1]+(1,))
                                                .broadcast_to(eq.shape))
        r = Rotation(np.squeeze(np.take_along_axis(eq.quaternion,
                                                   np.argmin(m,axis=0)[np.newaxis,...,np.newaxis],
                                                   axis=0),
                                axis=0))

        if return_cloud:
            return AverageTuple(self.copy(Rotation(r).average(weights)),self.copy(Rotation(r)))
        else:
            return self.copy(Rotation(r).average(weights))

    @overload
    def to_SST(self, vector: FloatSequence, proper: bool = False,
               return_operator: Literal[False] = False,
               return_operators: bool = False) -> np.ndarray:
        ...
    @overload
    def to_SST(self, vector: FloatSequence, proper: bool = False,
               return_operator: Literal[True] = True,
               return_operators: bool = False) -> ToSSTTuple:
        ...
    def to_SST(self,
               vector: FloatSequence,
               proper: bool = False,
               return_operator: bool = False,
               return_operators: bool = False) -> Union[np.ndarray, ToSSTTuple]:
        """
        Rotate lab frame vector to ensure it falls into (improper or proper) standard stereographic triangle of crystal symmetry.

        Parameters
        ----------
        vector : numpy.ndarray, shape (...,3)
            Lab frame vector to align with an SST crystal frame direction.
            Shape of vector blends with shape of own rotation array.
            For example, a rotation array of shape (3,2) and a vector array of shape (2,4) result in (3,2,4) outputs.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs.
            Defaults to False.
        return_operator : bool, optional
            Return the index of the symmetrically equivalent orientation that rotated vector to SST.
            Defaults to False.

        Returns
        -------
        vector_SST : numpy.ndarray, shape (...,3)
            Rotated vector falling into SST.
        operator : numpy.ndarray of int, shape (...), conditional
            Index of the symmetrically equivalent orientation that rotated vector to SST.
        """
        if return_operators:
            warnings.warn('"return_operators" is deprecated, use "return_operator".',DeprecationWarning,stacklevel=2)
        vector_ = np.array(vector,float)
        if vector_.shape[-1] != 3:
            raise ValueError('input is not a field of three-dimensional vectors')

        blend = util.shapeblender(self.shape,vector_.shape[:-1])
        eq    = self.broadcast_to(util.shapeshifter(self.shape,blend,mode='right')).equivalent
        poles = np.atleast_2d(eq @ np.broadcast_to(vector_,(1,)+blend+(3,)))
        ok    = self.in_SST(poles,proper=proper)
        ok   &= np.cumsum(ok,axis=0) == 1
        loc   = np.where(ok)
        sort  = 0 if len(loc) == 1 else np.lexsort(loc[:0:-1])

        if return_operator:
            return ToSSTTuple(poles[ok][sort].reshape(blend+(3,)), (np.vstack(loc[:1]).T)[sort].reshape(blend))
        else:
            return poles[ok][sort].reshape(blend+(3,))


    def in_SST(self,
               vector: FloatSequence,
               proper: bool = False) -> Union[np.bool_, np.ndarray]:
        """
        Check whether given crystal frame vector falls into standard stereographic triangle of own symmetry.

        Parameters
        ----------
        vector : numpy.ndarray, shape (...,3)
            Vector to check.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs.
            Defaults to False.

        Returns
        -------
        in : numpy.ndarray, shape (...)
            Whether vector falls into SST.
        """
        vector_ = np.array(vector,float)
        if vector_.shape[-1] != 3:
            raise ValueError('input is not a field of three-dimensional vectors')

        if self.standard_triangle is None:                                                          # direct exit for no symmetry
            return np.ones_like(vector_[...,0],bool)

        if proper:
            components_proper   = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(self.standard_triangle['proper'],   vector_.shape+(3,)),
                                                      vector_), 12)
            components_improper = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(self.standard_triangle['improper'], vector_.shape+(3,)),
                                                      vector_), 12)
            return   np.all(components_proper   >= 0.0,axis=-1) \
                   | np.all(components_improper >= 0.0,axis=-1)
        else:
            components = np.around(np.einsum('...ji,...i',
                                             np.broadcast_to(self.standard_triangle['improper'], vector_.shape+(3,)),
                                             np.block([vector_[...,:2],np.abs(vector_[...,2:3])])), 12)

            return np.all(components >= 0.0,axis=-1)


    def IPF_color(self,
                  vector: FloatSequence,
                  in_SST: bool = True,
                  proper: bool = False) -> np.ndarray:
        """
        Map lab frame vector to RGB color within standard stereographic triangle of own symmetry.

        Parameters
        ----------
        vector : numpy.ndarray, shape (...,3)
            Lab frame vector to colorize.
            Shape of vector blends with shape of own rotation array.
            For example, a rotation array of shape (3,2) and a vector array of shape (2,4) result in (3,2,4) outputs.
        in_SST : bool, optional
            Consider symmetrically equivalent orientations such that poles are located in SST.
            Defaults to True.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs (with mirrored colors).
            Defaults to False.

        Returns
        -------
        rgb : numpy.ndarray, shape (...,3)
           RGB array of IPF colors.

        Examples
        --------
        Inverse pole figure color of the e_3 lab direction for a
        crystal in "Cube" orientation with cubic symmetry:

        >>> import damask
        >>> o = damask.Orientation(family='cubic')
        >>> o.IPF_color([0,0,1])
        array([1., 0., 0.])

        Sample standard triangle for hexagonal symmetry:

        >>> import damask
        >>> from matplotlib import pyplot as plt
        >>> lab = [0,0,1]
        >>> o = damask.Orientation.from_random(shape=500000,family='hexagonal')
        >>> coord = damask.util.project_equal_area(o.to_SST(lab))
        >>> color = o.IPF_color(lab)
        >>> plt.scatter(coord[:,0],coord[:,1],color=color,s=.06)
        <matplotlib.collections.PathCollection object at ...>
        >>> plt.axis('scaled')
        (...)
        >>> plt.show(block=False)
        """
        if np.array(vector).shape[-1] != 3:
            raise ValueError('input is not a field of three-dimensional vectors')

        vector_:np.ndarray = self.to_SST(vector,proper) if in_SST else \
                             self @ np.broadcast_to(vector,self.shape+(3,))

        if self.standard_triangle is None:                                                          # direct exit for no symmetry
            return np.zeros_like(vector_)

        if proper:
            components_proper   = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(self.standard_triangle['proper'],   vector_.shape+(3,)),
                                                      vector_), 12)
            components_improper = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(self.standard_triangle['improper'], vector_.shape+(3,)),
                                                      vector_), 12)
            in_SST_ = np.logical_or(np.all(components_proper   >= 0.0,axis=-1),
                                    np.all(components_improper >= 0.0,axis=-1))
            components = np.where((np.logical_and(in_SST_,np.all(components_proper >= 0.0,axis=-1)))[...,np.newaxis],
                                  components_proper,components_improper)
        else:
            components = np.around(np.einsum('...ji,...i',
                                             np.broadcast_to(self.standard_triangle['improper'], vector_.shape+(3,)),
                                             np.block([vector_[...,:2],np.abs(vector_[...,2:3])])), 12)

            in_SST_ = np.all(components >= 0.0,axis=-1)

        with np.errstate(invalid='ignore',divide='ignore'):
            rgb = (components/np.linalg.norm(components,axis=-1,keepdims=True))**(1./3.)            # smoothen color ramps
            rgb = np.clip(rgb,0.,1.)                                                                # clip intensity
            rgb /= np.max(rgb,axis=-1,keepdims=True)                                                # normalize to (HS)V = 1
        rgb[np.broadcast_to(~in_SST_[...,np.newaxis],rgb.shape)] = 0.0

        return rgb


####################################################################################################
    # functions that require lattice, not just family

    def to_lattice(self, *,
                   direction: Optional[FloatSequence] = None,
                   plane: Optional[FloatSequence] = None) -> np.ndarray:                            # numpydoc ignore=PR01,PR02
        """
        Calculate lattice vector corresponding to lab frame direction or plane normal.

        Parameters
        ----------
        direction|plane : numpy.ndarray, shape (...,3)
            Real space vector along direction or
            reciprocal space vector along plane normal.
            Shape of vector blends with shape of own rotation array.
            For example, a rotation array of shape (3,2) and a vector
            array of shape (2,4) result in (3,2,4) outputs.

        Returns
        -------
        Miller : numpy.ndarray, shape (...,3)
            Lattice vector of direction or plane.
            Use util.scale_to_coprime to convert to (integer) Miller indices.

        Examples
        --------
        >>> import numpy as np
        >>> import damask
        >>> cubic = damask.Orientation.from_axis_angle(n_omega=[1,0,0,90],degrees=True,lattice='cI')
        >>> np.round(cubic.to_lattice(direction=[1, 0, 0]),3)
        array([1., 0., 0.])
        >>> np.round(cubic.to_lattice(direction=[0, 1, 0]),3)
        array([ 0., 0., -1.])
        >>> np.round(cubic.to_lattice(direction=[0, 0, 1]),3)
        array([-0., 1., 0.])
        >>> tetragonal = damask.Orientation(lattice='tI',c=0.5)
        >>> damask.util.scale_to_coprime(tetragonal.to_lattice(direction=[1,1,1]))
        array([1, 1, 2])
        >>> damask.util.scale_to_coprime(tetragonal.to_lattice(plane=[1,1,1]))
        array([2, 2, 1])
        """
        if (direction is not None) ^ (plane is None):
            raise KeyError('specify either "direction" or "plane"')
        return (super().to_lattice(direction=self@np.asarray(direction)) if plane is None else
                super().to_lattice(plane=self@np.asarray(plane)))


    def to_frame(self, *,
                 uvw: Optional[IntSequence] = None,
                 hkl: Optional[IntSequence] = None,
                 uvtw: Optional[IntSequence] = None,
                 hkil: Optional[IntSequence] = None,
                 with_symmetry: bool = False,
                 normalize: bool = True,
                 ) -> np.ndarray:                                                                   # numpydoc ignore=PR01,PR02
        """
        Calculate lab frame vector along lattice direction [uvw]/[uvtw] or plane normal (hkl)/(hkil).

        Parameters
        ----------
        uvw|hkl|uvtw|hkil : numpy.ndarray, shape (...,3) or shape (...,4)
            Miller(–Bravais) indices of crystallographic direction or plane normal.
            Shape of vector blends with shape of own rotation array.
            For example, a rotation array of shape (3,2) and a vector
            array of shape (2,4) result in (3,2,4) outputs.
        with_symmetry : bool, optional
            Calculate all N symmetrically equivalent vectors.
            Defaults to False.
        normalize : bool, optional
            Normalize output vector.
            Defaults to True.

        Returns
        -------
        vector : numpy.ndarray, shape (...,3) or (N,...,3)
            Lab frame vector (or N vectors if with_symmetry) along
            [uvw]/[uvtw] direction or (hkl)/(hkil) plane normal.
        """
        v = super().to_frame(uvw=uvw,hkl=hkl,uvtw=uvtw,hkil=hkil)
        s_v = v.shape[:-1]
        blend = util.shapeblender(self.shape,s_v)
        if normalize:
            v /= np.linalg.norm(v,axis=-1,keepdims=len(s_v)>0)
        if with_symmetry:
            sym_ops = self.symmetry_operations
            s_v   += sym_ops.shape
            blend += sym_ops.shape
            v = sym_ops.broadcast_to(s_v) @ v[...,np.newaxis,:]

        return np.moveaxis(~(self.broadcast_to(blend)) @ np.broadcast_to(v,blend+(3,)),
                           -2 if with_symmetry else 0,
                           0)


    def Schmid(self, *,
               N_slip: Optional[Union[IntSequence, Literal['*']]] = None,
               N_twin: Optional[Union[IntSequence, Literal['*']]] = None) -> np.ndarray:            # numpydoc ignore=PR01,PR02
        u"""
        Calculate Schmid matrix P = d ⨂ n in the lab frame for selected deformation systems.

        Parameters
        ----------
        N_slip|N_twin : '*' or sequence of int
            Number of deformation systems per family of the deformation system.
            Use '*' to select all.

        Returns
        -------
        P : numpy.ndarray, shape (N,...,3,3)
            Schmid matrix for each of the N deformation systems.

        Examples
        --------
        Schmid matrix (in lab frame) of first octahedral slip system of a face-centered
        cubic crystal in "Goss" orientation.

        >>> import numpy as np
        >>> import damask
        >>> O = damask.Orientation.from_Euler_angles(phi=[0,45,0],degrees=True,lattice='cF')
        >>> np.round(O.Schmid(N_slip=[12])[0],3)
        array([[ 0.   ,  0.   ,  0.   ],
               [ 0.577, -0.   ,  0.816],
               [ 0.   ,  0.   ,  0.   ]])

        Schmid matrix (in lab frame) of first {110}, {112}, and {123} slip systems of
        a body-centered cubic crystal in "rotated cube" orientation.

        >>> import numpy as np
        >>> import damask
        >>> O = damask.Orientation.from_directions([0,0,1],[1,0,0],lattice='cI')
        >>> np.round(O.Schmid(N_slip=[1,1,1]),3)
        array([[[ 0.408, -0.408, -0.   ],
                [ 0.408, -0.408, -0.   ],
                [ 0.408, -0.408, -0.   ]],
        <BLANKLINE>
               [[-0.236, -0.236,  0.471],
                [-0.236, -0.236,  0.471],
                [-0.236, -0.236,  0.471]],
        <BLANKLINE>
               [[ 0.463, -0.309, -0.154],
                [ 0.463, -0.309, -0.154],
                [ 0.463, -0.309, -0.154]]])

        Schmid matrix (in lab frame) of all three prismatic slip systems of a
        hexagonal crystal in "cube" orientation.

        >>> import numpy as np
        >>> import damask
        >>> O = damask.Orientation(lattice='hP')
        >>> np.round(O.Schmid(N_slip=[0,3]),3)
        array([[[ 0.   ,  1.   , -0.   ],
                [ 0.   ,  0.   ,  0.   ],
                [ 0.   ,  0.   ,  0.   ]],
        <BLANKLINE>
               [[ 0.433,  0.25 , -0.   ],
                [-0.75 , -0.433,  0.   ],
                [ 0.   ,  0.   ,  0.   ]],
        <BLANKLINE>
               [[-0.433,  0.25 , -0.   ],
                [-0.75 ,  0.433, -0.   ],
                [ 0.   ,  0.   ,  0.   ]]])
        """
        return np.moveaxis(~self @
                           super().Schmid(N_slip=N_slip,
                                          N_twin=N_twin)[(np.newaxis,)*len(self.shape)],
                           -3,
                           0)


    def related(self: MyType,
                model: str,
                target = None) -> MyType:
        """
        All orientations related to self by given relationship model.

        Parameters
        ----------
        model : str
            Orientation relationship model selected from self.orientation_relationships.
        target : Crystal, optional
            Crystal to transform to.
            Providing this parameter allows specification of non-standard lattice parameters.
            Default is inferred from selected model and uses standard lattice parameters.

        Returns
        -------
        rel : Orientation, shape (:,self.shape)
            Orientations related to self according to the selected
            model for the orientation relationship.

        Examples
        --------
        Face-centered cubic orientations following from a
        body-centered cubic crystal in "Cube" orientation according
        to the Bain orientation relationship (cF -> cI).

        >>> import numpy as np
        >>> import damask
        >>> damask.Orientation(lattice='cF').related('Bain')
        Crystal family: cubic
        Bravais lattice: cI
        a=1 m, b=1 m, c=1 m
        α=90°, β=90°, γ=90°
        array([( 6.53281482e-01,     2.70598050e-01,  6.53281482e-01,  2.70598050e-01),
               ( 2.70598050e-01,    -2.70598050e-01, -6.53281482e-01, -6.53281482e-01),
               ( 9.23879533e-01,    -5.55111512e-17, -2.77555756e-17, -3.82683432e-01)])
        """
        lattice,o = self.relation_operations(model,target)
        target = Crystal(lattice=lattice) if target is None else target
        return Orientation(rotation=o*Rotation(self.quaternion)[np.newaxis,...],                    # type: ignore[return-value]
                           lattice=target.lattice,
                           a=target.a,
                           b=target.b,
                           c=target.c,
                           alpha=target.alpha,
                           beta =target.beta,
                           gamma=target.gamma,
                           )
