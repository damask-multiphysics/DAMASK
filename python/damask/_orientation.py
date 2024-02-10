import copy
from typing import Optional, Union, TypeVar

import numpy as np

from ._typehints import FloatSequence, IntSequence, CrystalFamily, BravaisLattice
from . import Rotation
from . import Crystal
from . import util
from . import tensor

MyType = TypeVar('MyType', bound='Orientation')

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
    "equivalent", "reduced", "disorientation", "IPF_color", or "to_SST".

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
    extends the functionality of Orientation objects to include operations such as
    "Schmid", "related", or "to_pole" that require a lattice type and its parameters.

    Examples
    --------
    An array of 3 x 5 random orientations reduced to the fundamental zone of tetragonal symmetry:

    >>> import damask
    >>> o=damask.Orientation.from_random(shape=(3,5),family='tetragonal').reduced

    """

    @util.extend_docstring(adopted_parameters=Crystal.__init__)
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

        """
        Rotation.__init__(self,rotation)
        Crystal.__init__(self,family=family, lattice=lattice,
                              a=a,b=b,c=c, alpha=alpha,beta=beta,gamma=gamma, degrees=degrees)


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.

        """
        return util.srepr([Crystal.__repr__(self),
                           Rotation.__repr__(self)])


    def __copy__(self: MyType,
                 rotation: Union[None, FloatSequence, Rotation] = None) -> MyType:
        """
        Return deepcopy(self).

        Create deep copy.

        """
        dup = copy.deepcopy(self)
        if rotation is not None:
            dup.quaternion = Rotation(rotation).quaternion
        return dup

    copy = __copy__



    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : Orientation
            Orientation to check for equality.

        """
        if not isinstance(other, Orientation):
            return NotImplemented
        matching_type = self.family == other.family and \
                        self.lattice == other.lattice and \
                        self.parameters == other.parameters
        return np.logical_and(matching_type,super(self.__class__,self.reduced).__eq__(other.reduced))

    def __ne__(self,
               other: object) -> bool:
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
                rtol: float = 1e-5,
                atol: float = 1e-8,
                equal_nan: bool = True) -> bool:
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
                 rtol: float = 1e-5,
                 atol: float = 1e-8,
                 equal_nan: bool = True) -> bool:
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
        answer : bool
            Whether all values are close between both orientations.

        """
        return bool(np.all(self.isclose(other,rtol,atol,equal_nan)))


    def __mul__(self: MyType,
                other: Union[Rotation, 'Orientation']) -> MyType:
        """
        Return self*other.

        Compose with other.

        Parameters
        ----------
        other : Rotation or Orientation, shape (self.shape)
            Object for composition.

        Returns
        -------
        composition : Orientation
            Compound rotation self*other, i.e. first other then self rotation.

        """
        if not isinstance(other, (Orientation,Rotation)):
            raise TypeError('use "O@b", i.e. matmul, to apply Orientation "O" to object "b"')
        return self.copy(Rotation(self.quaternion)*Rotation(other.quaternion))


    @classmethod
    @util.extend_docstring(Rotation.from_quaternion,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_quaternion, wrapped=__init__)
    def from_quaternion(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_Euler_angles,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_Euler_angles, wrapped=__init__)
    def from_Euler_angles(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_axis_angle,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_axis_angle, wrapped=__init__)
    def from_axis_angle(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_basis,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_basis, wrapped=__init__)
    def from_basis(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_matrix,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_matrix, wrapped=__init__)
    def from_matrix(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_parallel,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_parallel, wrapped=__init__)
    def from_parallel(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_Rodrigues_vector,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_Rodrigues_vector, wrapped=__init__)
    def from_Rodrigues_vector(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_homochoric,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_homochoric, wrapped=__init__)
    def from_homochoric(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_cubochoric,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_cubochoric, wrapped=__init__)
    def from_cubochoric(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_random,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_random, wrapped=__init__)
    def from_random(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_ODF,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_ODF, wrapped=__init__)
    def from_ODF(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_spherical_component,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_spherical_component, wrapped=__init__)
    def from_spherical_component(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)

    @classmethod
    @util.extend_docstring(Rotation.from_fiber_component,
                           adopted_parameters=Crystal.__init__)
    @util.pass_on('rotation', Rotation.from_fiber_component, wrapped=__init__)
    def from_fiber_component(cls, **kwargs) -> 'Orientation':
        return cls(**kwargs)


    @classmethod
    @util.extend_docstring(adopted_parameters=Crystal.__init__)
    def from_directions(cls,
                        uvw: FloatSequence,
                        hkl: FloatSequence,
                        **kwargs) -> 'Orientation':
        """
        Initialize orientation object from the crystallographic direction and plane parallel to lab x and z, respectively.

        Parameters
        ----------
        uvw : numpy.ndarray, shape (...,3)
            Lattice direction aligned with lab frame x-direction.
        hkl : numpy.ndarray, shape (...,3)
            Lattice plane normal aligned with lab frame z-direction.

        Returns
        -------
        new : damask.Orientation

        """
        o = cls(**kwargs)
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

        rho  = self.as_Rodrigues_vector(compact=True)
        return larger_or_equal(rho,
                                     [rho[...,1],           rho[...,2],0] if self.family == 'cubic'
                                else [rho[...,1]*np.sqrt(3),0,         0] if self.family == 'hexagonal'
                                else [rho[...,1],           0,         0] if self.family == 'tetragonal'
                                else [0,                    0,         0] if self.family == 'orthorhombic'
                                else [-np.inf,              0,         0] if self.family == 'monoclinic'
                                else [-np.inf,        -np.inf,   -np.inf]) & self.in_FZ


    def disorientation(self,
                       other: 'Orientation',
                       return_operators: bool = False) -> object:
        """
        Calculate disorientation between self and given other orientation.

        Parameters
        ----------
        other : Orientation
            Orientation to calculate disorientation for.
            Shape of other blends with shape of own rotation array.
            For example, shapes of (2,3) for own rotations
            and (3,2) for other's result in (2,3,2) disorientations.
        return_operators : bool, optional
            Return index pair of symmetrically equivalent orientations
            that result in disorientation axis falling into FZ.
            Defaults to False.

        Returns
        -------
        disorientation : Orientation
            Disorientation between self and other.
        operators : numpy.ndarray of int, shape (...,2), conditional
            Index of symmetrically equivalent orientation that rotated vector to the SST.

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
        Crystal family hexagonal
        Quaternion: (real=0.976, imag=<+0.189, +0.018, +0.103>)
        Matrix:
        [[ 0.97831006  0.20710935  0.00389135]
         [-0.19363288  0.90765544  0.37238141]
         [ 0.07359167 -0.36505797  0.92807163]]
        Bunge Eulers / deg: (11.40, 21.86, 0.60)

        Plot a sample from the Mackenzie distribution.

        >>> import matplotlib.pyplot as plt
        >>> import damask
        >>> N = 10000
        >>> a = damask.Orientation.from_random(shape=N,family='cubic')
        >>> b = damask.Orientation.from_random(shape=N,family='cubic')
        >>> n,omega = a.disorientation(b).as_axis_angle(degrees=True,pair=True)
        >>> plt.hist(omega,25)
        >>> plt.show()

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

        r_ = s[:,np.newaxis,...].misorientation(o[np.newaxis,:,...]) # type: ignore[index]
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

        return (
                (self.copy(rotation=quat),
                 (np.vstack(loc[:2]).T)[sort].reshape((*shp,2)))
                if return_operators else
                self.copy(rotation=quat)
               )


    def average(self,
                weights: Optional[FloatSequence] = None,
                return_cloud: bool = False):
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
        m  = eq.misorientation(self[...,0].reshape((1,)+self.shape[:-1]+(1,))             # type: ignore
                                          .broadcast_to(eq.shape)).as_axis_angle()[...,3] # type: ignore
        r = Rotation(np.squeeze(np.take_along_axis(eq.quaternion,
                                                   np.argmin(m,axis=0)[np.newaxis,...,np.newaxis],
                                                   axis=0),
                                axis=0))

        return ((self.copy(Rotation(r).average(weights)),self.copy(Rotation(r))) if return_cloud else
                self.copy(Rotation(r).average(weights))
               )


    def to_SST(self,
               vector: FloatSequence,
               proper: bool = False,
               return_operators: bool = False) -> np.ndarray:
        """
        Rotate lab frame vector to ensure it falls into (improper or proper) standard stereographic triangle of crystal symmetry.

        Parameters
        ----------
        vector : numpy.ndarray, shape (...,3)
            Lab frame vector to align with crystal frame direction.
            Shape of vector blends with shape of own rotation array.
            For example, a rotation array of shape (3,2) and a vector array of shape (2,4) result in (3,2,4) outputs.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs.
            Defaults to False.
        return_operators : bool, optional
            Return the symmetrically equivalent orientation that rotated vector to SST.
            Defaults to False.

        Returns
        -------
        vector_SST : numpy.ndarray, shape (...,3)
            Rotated vector falling into SST.
        operator : numpy.ndarray of int, shape (...), conditional
            Index of symmetrically equivalent orientation that rotated vector to SST.

        """
        vector_ = np.array(vector,float)
        if vector_.shape[-1] != 3:
            raise ValueError('input is not a field of three-dimensional vectors')

        blend = util.shapeblender( self.shape,vector_.shape[:-1])
        eq    = self.broadcast_to(util.shapeshifter( self.shape,blend,mode='right')).equivalent
        poles = np.atleast_2d(eq @ np.broadcast_to(vector_,(1,)+blend+(3,)))
        ok    = self.in_SST(poles,proper=proper)
        ok   &= np.cumsum(ok,axis=0) == 1
        loc   = np.where(ok)
        sort  = 0 if len(loc) == 1 else np.lexsort(loc[:0:-1])

        return (
                (poles[ok][sort].reshape(blend+(3,)), (np.vstack(loc[:1]).T)[sort].reshape(blend))
                if return_operators else
                poles[ok][sort].reshape(blend+(3,))
               )


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
        >>> plt.axis('scaled')
        >>> plt.show()

        """
        if np.array(vector).shape[-1] != 3:
            raise ValueError('input is not a field of three-dimensional vectors')

        vector_ = self.to_SST(vector,proper) if in_SST else \
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
            in_SST_ = np.all(components_proper   >= 0.0,axis=-1) \
                    | np.all(components_improper >= 0.0,axis=-1)
            components = np.where((in_SST_ & np.all(components_proper   >= 0.0,axis=-1))[...,np.newaxis],
                                  components_proper,components_improper)
        else:
            components = np.around(np.einsum('...ji,...i',
                                             np.broadcast_to(self .standard_triangle['improper'], vector_.shape+(3,)),
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

    def to_pole(self, *,
                uvw: Optional[FloatSequence] = None,
                hkl: Optional[FloatSequence] = None,
                with_symmetry: bool = False,
                normalize: bool = True) -> np.ndarray:
        """
        Calculate lab frame vector along lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        uvw|hkl : numpy.ndarray, shape (...,3)
            Miller indices of crystallographic direction or plane normal.
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
        vector : numpy.ndarray, shape (...,3) or (...,N,3)
            Lab frame vector (or vectors if with_symmetry) along
            [uvw] direction or (hkl) plane normal.

        """
        v = self.to_frame(uvw=uvw,hkl=hkl)
        s_v = v.shape[:-1]
        blend = util.shapeblender(self.shape,s_v)
        if normalize:
            v /= np.linalg.norm(v,axis=-1,keepdims=len(s_v)>0)
        if with_symmetry:
            sym_ops = self.symmetry_operations
            s_v   += sym_ops.shape
            blend += sym_ops.shape
            v = sym_ops.broadcast_to(s_v) @ v[...,np.newaxis,:]

        return ~(self.broadcast_to(blend)) @ np.broadcast_to(v,blend+(3,))


    def Schmid(self, *,
               N_slip: Optional[IntSequence] = None,
               N_twin: Optional[IntSequence] = None) -> np.ndarray:
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
        >>> np.set_printoptions(3,suppress=True,floatmode='fixed')
        >>> O = damask.Orientation.from_Euler_angles(phi=[0,45,0],degrees=True,lattice='cF')
        >>> O.Schmid(N_slip=[12])[0]
        array([[ 0.000,  0.000,  0.000],
               [ 0.577, -0.000,  0.816],
               [ 0.000,  0.000,  0.000]])

        """
        if (N_slip is not None) ^ (N_twin is None):
            raise KeyError('specify either "N_slip" or "N_twin"')

        kinematics,active = (self.kinematics('slip'),N_slip) if N_twin is None else \
                            (self.kinematics('twin'),N_twin)
        if active == '*': active = [len(a) for a in kinematics['direction']]

        if not active:
            raise ValueError('Schmid matrix not defined')
        d = self.to_frame(uvw=np.vstack([kinematics['direction'][i][:n] for i,n in enumerate(active)]))
        p = self.to_frame(hkl=np.vstack([kinematics['plane'][i][:n] for i,n in enumerate(active)]))
        P = np.einsum('...i,...j',d/np.linalg.norm(d,axis=1,keepdims=True),
                                  p/np.linalg.norm(p,axis=1,keepdims=True))

        shape = P.shape[0:1]+self.shape+(3,3)

        return ~self.broadcast_to(shape[:-2]) \
               @ np.broadcast_to(P.reshape(util.shapeshifter(P.shape,shape)),shape)


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
        to the Bain orientation relationship (cI -> cF).

        >>> import numpy as np
        >>> import damask
        >>> np.set_printoptions(3,suppress=True,floatmode='fixed')
        >>> damask.Orientation(lattice='cI').related('Bain')
        Crystal family: cubic
        Bravais lattice: cF
        a=1 m, b=1 m, c=1 m
        α=90°, β=90°, γ=90°
        Quaternions of shape (3,)
        [[0.924 0.383 0.000 0.000]
         [0.924 0.000 0.383 0.000]
         [0.924 0.000 0.000 0.383]]

        """
        lattice,o = self.relation_operations(model,target)
        target = Crystal(lattice=lattice) if target is None else target

        return Orientation(rotation=o*Rotation(self.quaternion)[np.newaxis,...],  # type: ignore
                          lattice=lattice,
                          b = self.b if target.ratio['b'] is None else self.a*target.ratio['b'],
                          c = self.c if target.ratio['c'] is None else self.a*target.ratio['c'],
                          alpha = None if 'alpha' in target.immutable else self.alpha,
                          beta  = None if 'beta'  in target.immutable else self.beta,
                          gamma = None if 'gamma' in target.immutable else self.gamma,
                         )
