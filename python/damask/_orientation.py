import inspect

import numpy as np

from . import Rotation
from . import util
from . import tensor

_parameter_doc = \
       """lattice : str
            Either a crystal family  out of [triclinic, monoclinic, orthorhombic, tetragonal, hexagonal, cubic]
            or     a Bravais lattice out of [aP, mP, mS, oP, oS, oI, oF, tP, tI, hP, cP, cI, cF].
            When specifying a Bravais lattice, additional lattice parameters might be required:
        a : float, optional
            Length of lattice parameter "a".
        b : float, optional
            Length of lattice parameter "b".
        c : float, optional
            Length of lattice parameter "c".
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.

       """


class Orientation(Rotation):
    """
    Representation of crystallographic orientation as combination of rotation and either crystal family or Bravais lattice.

    The crystal family is one of Orientation.crystal_families:

    - triclinic
    - monoclinic
    - orthorhombic
    - tetragonal
    - hexagonal
    - cubic

    and enables symmetry-related operations such as
    "equivalent", "reduced", "disorientation", "IPF_color", or "to_SST".

    The Bravais lattice is one of Orientation.lattice_symmetries:

    - triclinic
       - aP : primitive

    - monoclininic
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

    >>> damask.Orientation.from_random(shape=(3,5),lattice='tetragonal').reduced

    """

    crystal_families = ['triclinic',
                        'monoclinic',
                        'orthorhombic',
                        'tetragonal',
                        'hexagonal',
                        'cubic']

    lattice_symmetries = {
                'aP': 'triclinic',

                'mP': 'monoclinic',
                'mS': 'monoclinic',

                'oP': 'orthorhombic',
                'oS': 'orthorhombic',
                'oI': 'orthorhombic',
                'oF': 'orthorhombic',

                'tP': 'tetragonal',
                'tI': 'tetragonal',

                'hP': 'hexagonal',

                'cP': 'cubic',
                'cI': 'cubic',
                'cF': 'cubic',
               }


    @util.extend_docstring(_parameter_doc)
    def __init__(self,
                 rotation = None,
                 lattice = None,
                 a = None,b = None,c = None,
                 alpha = None,beta = None,gamma = None,
                 degrees = False):
        """
        New orientation.

        Parameters
        ----------
        rotation : list, numpy.ndarray, Rotation, optional
            Unit quaternion in positive real hemisphere.
            Use .from_quaternion to perform a sanity check.
            Defaults to no rotation.

        """
        from damask.lattice import kinematics

        Rotation.__init__(self) if rotation is None else Rotation.__init__(self,rotation=rotation)

        if (    lattice not in self.lattice_symmetries
            and lattice not in self.crystal_families):
            raise KeyError(f'Lattice "{lattice}" is unknown')

        self.family     = None
        self.lattice    = None
        self.a          = None
        self.b          = None
        self.c          = None
        self.alpha      = None
        self.beta       = None
        self.gamma      = None
        self.kinematics = None

        if lattice in self.lattice_symmetries:
            self.family  = self.lattice_symmetries[lattice]
            self.lattice = lattice
            self.a = 1 if a is None else a
            self.b = b
            self.c = c
            self.alpha = (np.radians(alpha) if degrees else alpha) if alpha is not None else None
            self.beta  = (np.radians(beta)  if degrees else beta)  if beta  is not None else None
            self.gamma = (np.radians(gamma) if degrees else gamma) if gamma is not None else None

            self.a = float(self.a) if self.a is not None else \
                     (self.b / self.ratio['b'] if self.b is not None and self.ratio['b'] is not None else
                      self.c / self.ratio['c'] if self.c is not None and self.ratio['c'] is not None else None)
            self.b = float(self.b) if self.b is not None else \
                     (self.a * self.ratio['b'] if self.a is not None and self.ratio['b'] is not None else
                      self.c / self.ratio['c'] * self.ratio['b']
                      if self.c is not None and self.ratio['b'] is not None and self.ratio['c'] is not None else None)
            self.c = float(self.c) if self.c is not None else \
                     (self.a * self.ratio['c'] if self.a is not None and self.ratio['c'] is not None else
                      self.b / self.ratio['b'] * self.ratio['c']
                      if self.c is not None and self.ratio['b'] is not None and self.ratio['c'] is not None else None)
            self.alpha = self.alpha if self.alpha is not None else self.immutable['alpha'] if 'alpha' in self.immutable else None
            self.beta  = self.beta  if self.beta  is not None else self.immutable['beta']  if 'beta'  in self.immutable else None
            self.gamma = self.gamma if self.gamma is not None else self.immutable['gamma'] if 'gamma' in self.immutable else None

            if \
                (self.a     is None) \
             or (self.b     is None or ('b'     in self.immutable and self.b     != self.immutable['b'] * self.a)) \
             or (self.c     is None or ('c'     in self.immutable and self.c     != self.immutable['c'] * self.b)) \
             or (self.alpha is None or ('alpha' in self.immutable and self.alpha != self.immutable['alpha'])) \
             or (self.beta  is None or ( 'beta' in self.immutable and self.beta  != self.immutable['beta'])) \
             or (self.gamma is None or ('gamma' in self.immutable and self.gamma != self.immutable['gamma'])):
                raise ValueError (f'Incompatible parameters {self.parameters} for crystal family {self.family}')

            if np.any(np.array([self.alpha,self.beta,self.gamma]) <= 0):
                raise ValueError ('Lattice angles must be positive')
            if np.any([np.roll([self.alpha,self.beta,self.gamma],r)[0]
              > np.sum(np.roll([self.alpha,self.beta,self.gamma],r)[1:]) for r in range(3)]):
                raise ValueError ('Each lattice angle must be less than sum of others')

            if self.lattice in kinematics:
                master = kinematics[self.lattice]
                self.kinematics = {}
                for m in master:
                    self.kinematics[m] = {'direction':master[m][:,0:3],'plane':master[m][:,3:6]} \
                                         if master[m].shape[-1] == 6 else \
                                         {'direction':self.Bravais_to_Miller(uvtw=master[m][:,0:4]),
                                          'plane':    self.Bravais_to_Miller(hkil=master[m][:,4:8])}
        elif lattice in self.crystal_families:
            self.family = lattice


    def __repr__(self):
        """Represent."""
        return '\n'.join(([] if self.lattice is None else [f'Bravais lattice {self.lattice}'])
                       + ([f'Crystal family {self.family}'])
                       + [super().__repr__()])


    def __copy__(self,**kwargs):
        """Create deep copy."""
        return self.__class__(rotation=kwargs['rotation'] if 'rotation' in kwargs else self.quaternion,
                              lattice =kwargs['lattice']  if 'lattice'  in kwargs else self.lattice
                                                                                    if self.lattice is not None else self.family,
                              a       =kwargs['a']        if 'a'        in kwargs else self.a,
                              b       =kwargs['b']        if 'b'        in kwargs else self.b,
                              c       =kwargs['c']        if 'c'        in kwargs else self.c,
                              alpha   =kwargs['alpha']    if 'alpha'    in kwargs else self.alpha,
                              beta    =kwargs['beta']     if 'beta'     in kwargs else self.beta,
                              gamma   =kwargs['gamma']    if 'gamma'    in kwargs else self.gamma,
                              degrees =kwargs['degrees']  if 'degrees'  in kwargs else None,
                             )

    copy = __copy__


    def __eq__(self,other):
        """
        Equal to other.

        Parameters
        ----------
        other : Orientation
            Orientation to check for equality.

        """
        matching_type = all([hasattr(other,attr) and getattr(self,attr) == getattr(other,attr)
                             for attr in ['family','lattice','parameters']])
        return np.logical_and(matching_type,super(self.__class__,self.reduced).__eq__(other.reduced))

    def __ne__(self,other):
        """
        Not equal to other.

        Parameters
        ----------
        other : Orientation
            Orientation to check for equality.

        """
        return np.logical_not(self==other)


    def isclose(self,other,rtol=1e-5,atol=1e-8,equal_nan=True):
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
        mask : numpy.ndarray bool
            Mask indicating where corresponding orientations are close.

        """
        matching_type = all([hasattr(other,attr) and getattr(self,attr) == getattr(other,attr)
                             for attr in ['family','lattice','parameters']])
        return np.logical_and(matching_type,super(self.__class__,self.reduced).isclose(other.reduced))



    def allclose(self,other,rtol=1e-5,atol=1e-8,equal_nan=True):
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
        return np.all(self.isclose(other,rtol,atol,equal_nan))


    def __mul__(self,other):
        """
        Compose this orientation with other.

        Parameters
        ----------
        other : Rotation or Orientation
            Object for composition.

        Returns
        -------
        composition : Orientation
            Compound rotation self*other, i.e. first other then self rotation.

        """
        if isinstance(other,Orientation) or isinstance(other,Rotation):
            return self.copy(rotation=Rotation.__mul__(self,Rotation(other.quaternion)))
        else:
            raise TypeError('Use "O@b", i.e. matmul, to apply Orientation "O" to object "b"')


    @staticmethod
    def _split_kwargs(kwargs,target):
        """
        Separate keyword arguments in 'kwargs' targeted at 'target' from general keyword arguments of Orientation objects.

        Parameters
        ----------
        kwargs : dictionary
            Contains all **kwargs.
        target: method
            Function to scan for kwarg signature.

        Returns
        -------
        rot_kwargs: dictionary
            Valid keyword arguments of 'target' function of Rotation class.
        ori_kwargs: dictionary
            Valid keyword arguments of Orientation object.

        """
        kws = ()
        for t in (target,Orientation.__init__):
            kws += ({key: kwargs[key] for key in set(inspect.signature(t).parameters) & set(kwargs)},)

        invalid_keys = set(kwargs)-(set(kws[0])|set(kws[1]))
        if invalid_keys:
            raise TypeError(f"{inspect.stack()[1][3]}() got an unexpected keyword argument '{invalid_keys.pop()}'")

        return kws


    @classmethod
    @util.extended_docstring(Rotation.from_random,_parameter_doc)
    def from_random(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_random)
        return cls(rotation=Rotation.from_random(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_quaternion,_parameter_doc)
    def from_quaternion(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_quaternion)
        return cls(rotation=Rotation.from_quaternion(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_Euler_angles,_parameter_doc)
    def from_Euler_angles(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_Euler_angles)
        return cls(rotation=Rotation.from_Euler_angles(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_axis_angle,_parameter_doc)
    def from_axis_angle(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_axis_angle)
        return cls(rotation=Rotation.from_axis_angle(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_basis,_parameter_doc)
    def from_basis(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_basis)
        return cls(rotation=Rotation.from_basis(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_matrix,_parameter_doc)
    def from_matrix(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_matrix)
        return cls(rotation=Rotation.from_matrix(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_Rodrigues_vector,_parameter_doc)
    def from_Rodrigues_vector(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_Rodrigues_vector)
        return cls(rotation=Rotation.from_Rodrigues_vector(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_homochoric,_parameter_doc)
    def from_homochoric(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_homochoric)
        return cls(rotation=Rotation.from_homochoric(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_cubochoric,_parameter_doc)
    def from_cubochoric(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_cubochoric)
        return cls(rotation=Rotation.from_cubochoric(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_spherical_component,_parameter_doc)
    def from_spherical_component(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_spherical_component)
        return cls(rotation=Rotation.from_spherical_component(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extended_docstring(Rotation.from_fiber_component,_parameter_doc)
    def from_fiber_component(cls,**kwargs):
        kwargs_rot,kwargs_ori = Orientation._split_kwargs(kwargs,Rotation.from_fiber_component)
        return cls(rotation=Rotation.from_fiber_component(**kwargs_rot),**kwargs_ori)


    @classmethod
    @util.extend_docstring(_parameter_doc)
    def from_directions(cls,uvw,hkl,**kwargs):
        """
        Initialize orientation object from two crystallographic directions.

        Parameters
        ----------
        uvw : list, numpy.ndarray of shape (...,3)
            lattice direction aligned with lab frame x-direction.
        hkl : list, numpy.ndarray of shape (...,3)
            lattice plane normal aligned with lab frame z-direction.

        """
        o = cls(**kwargs)
        x = o.to_frame(uvw=uvw)
        z = o.to_frame(hkl=hkl)
        om = np.stack([x,np.cross(z,x),z],axis=-2)
        return o.copy(rotation=Rotation.from_matrix(tensor.transpose(om/np.linalg.norm(om,axis=-1,keepdims=True))))


    @property
    def symmetry_operations(self):
        """Symmetry operations as Rotations."""
        if self.family == 'cubic':
            sym_quats =  [
                          [ 1.0,            0.0,            0.0,            0.0            ],
                          [ 0.0,            1.0,            0.0,            0.0            ],
                          [ 0.0,            0.0,            1.0,            0.0            ],
                          [ 0.0,            0.0,            0.0,            1.0            ],
                          [ 0.0,            0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2) ],
                          [ 0.0,            0.0,            0.5*np.sqrt(2),-0.5*np.sqrt(2) ],
                          [ 0.0,            0.5*np.sqrt(2), 0.0,            0.5*np.sqrt(2) ],
                          [ 0.0,            0.5*np.sqrt(2), 0.0,           -0.5*np.sqrt(2) ],
                          [ 0.0,            0.5*np.sqrt(2),-0.5*np.sqrt(2), 0.0            ],
                          [ 0.0,           -0.5*np.sqrt(2),-0.5*np.sqrt(2), 0.0            ],
                          [ 0.5,            0.5,            0.5,            0.5            ],
                          [-0.5,            0.5,            0.5,            0.5            ],
                          [-0.5,            0.5,            0.5,           -0.5            ],
                          [-0.5,            0.5,           -0.5,            0.5            ],
                          [-0.5,           -0.5,            0.5,            0.5            ],
                          [-0.5,           -0.5,            0.5,           -0.5            ],
                          [-0.5,           -0.5,           -0.5,            0.5            ],
                          [-0.5,            0.5,           -0.5,           -0.5            ],
                          [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                          [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                          [-0.5*np.sqrt(2), 0.0,            0.5*np.sqrt(2), 0.0            ],
                          [-0.5*np.sqrt(2), 0.0,           -0.5*np.sqrt(2), 0.0            ],
                          [-0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0,            0.0            ],
                          [-0.5*np.sqrt(2),-0.5*np.sqrt(2), 0.0,            0.0            ],
                        ]
        elif self.family == 'hexagonal':
            sym_quats =  [
                          [ 1.0,            0.0,            0.0,            0.0            ],
                          [-0.5*np.sqrt(3), 0.0,            0.0,           -0.5            ],
                          [ 0.5,            0.0,            0.0,            0.5*np.sqrt(3) ],
                          [ 0.0,            0.0,            0.0,            1.0            ],
                          [-0.5,            0.0,            0.0,            0.5*np.sqrt(3) ],
                          [-0.5*np.sqrt(3), 0.0,            0.0,            0.5            ],
                          [ 0.0,            1.0,            0.0,            0.0            ],
                          [ 0.0,           -0.5*np.sqrt(3), 0.5,            0.0            ],
                          [ 0.0,            0.5,           -0.5*np.sqrt(3), 0.0            ],
                          [ 0.0,            0.0,            1.0,            0.0            ],
                          [ 0.0,           -0.5,           -0.5*np.sqrt(3), 0.0            ],
                          [ 0.0,            0.5*np.sqrt(3), 0.5,            0.0            ],
                        ]
        elif self.family == 'tetragonal':
            sym_quats =  [
                          [ 1.0,            0.0,            0.0,            0.0            ],
                          [ 0.0,            1.0,            0.0,            0.0            ],
                          [ 0.0,            0.0,            1.0,            0.0            ],
                          [ 0.0,            0.0,            0.0,            1.0            ],
                          [ 0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                          [ 0.0,           -0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                          [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                          [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                        ]
        elif self.family == 'orthorhombic':
            sym_quats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                          [ 0.0,1.0,0.0,0.0 ],
                          [ 0.0,0.0,1.0,0.0 ],
                          [ 0.0,0.0,0.0,1.0 ],
                        ]
        elif self.family == 'monoclinic':
            sym_quats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                          [ 0.0,0.0,1.0,0.0 ],
                        ]
        elif self.family == 'triclinic':
            sym_quats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                        ]
        else:
            raise KeyError(f'Crystal family "{self.family}" is unknown')

        return Rotation.from_quaternion(sym_quats,accept_homomorph=True)


    @property
    def equivalent(self):
        """
        Orientations that are symmetrically equivalent.

        One dimension (length corresponds to number of symmetrically equivalent orientations)
        is added to the left of the Rotation array.

        """
        o = self.symmetry_operations.broadcast_to(self.symmetry_operations.shape+self.shape,mode='right')
        return self.copy(rotation=o*Rotation(self.quaternion).broadcast_to(o.shape,mode='left'))


    @property
    def reduced(self):
        """Select symmetrically equivalent orientation that falls into fundamental zone according to symmetry."""
        eq   = self.equivalent
        ok   = eq.in_FZ
        ok  &= np.cumsum(ok,axis=0) == 1
        loc  = np.where(ok)
        sort = 0 if len(loc) == 1 else np.lexsort(loc[:0:-1])
        return eq[ok][sort].reshape(self.shape)


    @property
    def in_FZ(self):
        """
        Check whether orientation falls into fundamental zone of own symmetry.

        Returns
        -------
        in : numpy.ndarray of quaternion.shape
           Boolean array indicating whether Rodrigues-Frank vector falls into fundamental zone.

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
            elif self.family == 'hexagonal':
                return (np.prod(1.  >= rho_abs,axis=-1) *
                                (2. >= np.sqrt(3)*rho_abs[...,0] + rho_abs[...,1]) *
                                (2. >= np.sqrt(3)*rho_abs[...,1] + rho_abs[...,0]) *
                                (2. >= np.sqrt(3)                + rho_abs[...,2])).astype(bool)
            elif self.family == 'tetragonal':
                return (np.prod(1.  >= rho_abs[...,:2],axis=-1) *
                        (np.sqrt(2) >= rho_abs[...,0] + rho_abs[...,1]) *
                        (np.sqrt(2) >= rho_abs[...,2] + 1.)).astype(bool)
            elif self.family == 'orthorhombic':
                return (np.prod(1. >= rho_abs,axis=-1)).astype(bool)
            elif self.family == 'monoclinic':
                return (1. >= rho_abs[...,1]).astype(bool)
            else:
                return np.all(np.isfinite(rho_abs),axis=-1)


    @property
    def in_disorientation_FZ(self):
        """
        Check whether orientation falls into fundamental zone of disorientations.

        Returns
        -------
        in : numpy.ndarray of quaternion.shape
           Boolean array indicating whether Rodrigues-Frank vector falls into disorientation FZ.

        References
        ----------
        A. Heinz and P. Neumann, Acta Crystallographica Section A 47:780-789, 1991
        https://doi.org/10.1107/S0108767391006864

        """
        rho = self.as_Rodrigues_vector(compact=True)*(1.0-1.0e-9)

        with np.errstate(invalid='ignore'):
            if   self.family == 'cubic':
                return ((rho[...,0] >= rho[...,1]) &
                        (rho[...,1] >= rho[...,2]) &
                        (rho[...,2] >= 0)).astype(bool)
            elif self.family == 'hexagonal':
                return ((rho[...,0] >= rho[...,1]*np.sqrt(3)) &
                        (rho[...,1] >= 0) &
                        (rho[...,2] >= 0)).astype(bool)
            elif self.family == 'tetragonal':
                return ((rho[...,0] >= rho[...,1]) &
                        (rho[...,1] >= 0) &
                        (rho[...,2] >= 0)).astype(bool)
            elif self.family == 'orthorhombic':
                return ((rho[...,0] >= 0) &
                        (rho[...,1] >= 0) &
                        (rho[...,2] >= 0)).astype(bool)
            elif self.family == 'monoclinic':
                return ((rho[...,1] >= 0) &
                        (rho[...,2] >= 0)).astype(bool)
            else:
                return np.ones_like(rho[...,0],dtype=bool)


    def relation_operations(self,model,return_lattice=False):
        """
        Crystallographic orientation relationships for phase transformations.

        Parameters
        ----------
        model : str
            Name of orientation relationship.
        return_lattice : bool, optional
            Return the target lattice in addition.

        Returns
        -------
        operations : Rotations
            Rotations characterizing the orientation relationship.

        References
        ----------
        S. Morito et al., Journal of Alloys and Compounds 577:s587-s592, 2013
        https://doi.org/10.1016/j.jallcom.2012.02.004

        K. Kitahara et al., Acta Materialia 54(5):1279-1288, 2006
        https://doi.org/10.1016/j.actamat.2005.11.001

        Y. He et al., Journal of Applied Crystallography 39:72-81, 2006
        https://doi.org/10.1107/S0021889805038276

        H. Kitahara et al., Materials Characterization 54(4-5):378-386, 2005
        https://doi.org/10.1016/j.matchar.2004.12.015

        Y. He et al., Acta Materialia 53(4):1179-1190, 2005
        https://doi.org/10.1016/j.actamat.2004.11.021

        """
        from damask.lattice import relations

        if model not in relations:
            raise KeyError(f'Orientation relationship "{model}" is unknown')
        r = relations[model]

        if self.lattice not in r:
            raise KeyError(f'Relationship "{model}" not supported for lattice "{self.lattice}"')

        sl = self.lattice
        ol = (set(r)-{sl}).pop()
        m = r[sl]
        o = r[ol]

        p_,_p = np.zeros(m.shape[:-1]+(3,)),np.zeros(o.shape[:-1]+(3,))
        p_[...,0,:] = m[...,0,:] if m.shape[-1] == 3 else self.Bravais_to_Miller(uvtw=m[...,0,0:4])
        p_[...,1,:] = m[...,1,:] if m.shape[-1] == 3 else self.Bravais_to_Miller(hkil=m[...,1,0:4])
        _p[...,0,:] = o[...,0,:] if o.shape[-1] == 3 else self.Bravais_to_Miller(uvtw=o[...,0,0:4])
        _p[...,1,:] = o[...,1,:] if o.shape[-1] == 3 else self.Bravais_to_Miller(hkil=o[...,1,0:4])

        return (Rotation.from_parallel(p_,_p),ol) \
                if return_lattice else \
                Rotation.from_parallel(p_,_p)


    def related(self,model):
        """
        Orientations derived from the given relationship.

        One dimension (length according to number of related orientations)
        is added to the left of the Rotation array.

        """
        o,lattice = self.relation_operations(model,return_lattice=True)
        target = Orientation(lattice=lattice)
        o = o.broadcast_to(o.shape+self.shape,mode='right')
        return self.copy(rotation=o*Rotation(self.quaternion).broadcast_to(o.shape,mode='left'),
                         lattice=lattice,
                         b = self.b if target.ratio['b'] is None else self.a*target.ratio['b'],
                         c = self.c if target.ratio['c'] is None else self.a*target.ratio['c'],
                         alpha = None if 'alpha' in target.immutable else self.alpha,
                         beta  = None if 'beta'  in target.immutable else self.beta,
                         gamma = None if 'gamma' in target.immutable else self.gamma,
                        )


    @property
    def parameters(self):
        """Return lattice parameters a, b, c, alpha, beta, gamma."""
        return (self.a,self.b,self.c,self.alpha,self.beta,self.gamma)


    @property
    def immutable(self):
        """Return immutable parameters of own lattice."""
        if self.family == 'triclinic':
            return {}
        if self.family == 'monoclinic':
            return {
                     'alpha': np.pi/2.,
                     'gamma': np.pi/2.,
                   }
        if self.family == 'orthorhombic':
            return {
                     'alpha': np.pi/2.,
                     'beta':  np.pi/2.,
                     'gamma': np.pi/2.,
                   }
        if self.family == 'tetragonal':
            return {
                     'b': 1.0,
                     'alpha': np.pi/2.,
                     'beta':  np.pi/2.,
                     'gamma': np.pi/2.,
                   }
        if self.family == 'hexagonal':
            return {
                     'b': 1.0,
                     'alpha': np.pi/2.,
                     'beta':  np.pi/2.,
                     'gamma': 2.*np.pi/3.,
                   }
        if self.family == 'cubic':
            return {
                     'b': 1.0,
                     'c': 1.0,
                     'alpha': np.pi/2.,
                     'beta':  np.pi/2.,
                     'gamma': np.pi/2.,
                   }


    @property
    def ratio(self):
        """Return axes ratios of own lattice."""
        _ratio = { 'hexagonal': {'c': np.sqrt(8./3.)}}

        return dict(b = self.immutable['b']
                        if 'b' in self.immutable else
                        _ratio[self.family]['b'] if self.family in _ratio and 'b' in _ratio[self.family] else None,
                    c = self.immutable['c']
                        if 'c' in self.immutable else
                        _ratio[self.family]['c'] if self.family in _ratio and 'c' in _ratio[self.family] else None,
                   )


    @property
    def basis_real(self):
        """
        Calculate orthogonal real space crystal basis.

        References
        ----------
        C.T. Young and J.L. Lytton, Journal of Applied Physics 43:1408–1417, 1972
        https://doi.org/10.1063/1.1661333

        """
        if None in self.parameters:
            raise KeyError('Missing crystal lattice parameters')
        return np.array([
                          [1,0,0],
                          [np.cos(self.gamma),np.sin(self.gamma),0],
                          [np.cos(self.beta),
                           (np.cos(self.alpha)-np.cos(self.beta)*np.cos(self.gamma))                     /np.sin(self.gamma),
                           np.sqrt(1 - np.cos(self.alpha)**2 - np.cos(self.beta)**2 - np.cos(self.gamma)**2
                                 + 2 * np.cos(self.alpha)    * np.cos(self.beta)    * np.cos(self.gamma))/np.sin(self.gamma)],
                         ],dtype=float).T \
             * np.array([self.a,self.b,self.c])


    @property
    def basis_reciprocal(self):
        """Calculate reciprocal (dual) crystal basis."""
        return np.linalg.inv(self.basis_real.T)


    def in_SST(self,vector,proper=False):
        """
        Check whether given crystal frame vector falls into standard stereographic triangle of own symmetry.

        Parameters
        ----------
        vector : numpy.ndarray of shape (...,3)
            Vector to check.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs.
            Defaults to False.

        Returns
        -------
        in : numpy.ndarray of shape (...)
           Boolean array indicating whether vector falls into SST.

        References
        ----------
        Bases are computed from

        >>> basis = {
        ...    'cubic' :       np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,1.]/np.sqrt(2.),                #              green
        ...                                            [1.,1.,1.]/np.sqrt(3.)]).T),           #              blue
        ...    'hexagonal' :   np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,0.],                            #              green
        ...                                            [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T),  #              blue
        ...    'tetragonal' :  np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,0.],                            #              green
        ...                                            [1.,1.,0.]/np.sqrt(2.)]).T),           #              blue
        ...    'orthorhombic': np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,0.],                            #              green
        ...                                            [0.,1.,0.]]).T),                       #              blue
        ...    }

        """
        if not isinstance(vector,np.ndarray) or vector.shape[-1] != 3:
            raise ValueError('Input is not a field of three-dimensional vectors.')

        if self.family == 'cubic':
            basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                           [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                           [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                       'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                           [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                           [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                    }
        elif self.family == 'hexagonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -np.sqrt(3.)   ,  0. ],
                                           [ 0.            ,  2.            ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  np.sqrt(3.)   ,  0. ],
                                           [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                    }
        elif self.family == 'tetragonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -1.            ,  0. ],
                                           [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  1.            ,  0. ],
                                           [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                    }
        elif self.family == 'orthorhombic':
            basis = {'improper':np.array([ [ 0., 0., 1.],
                                           [ 1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                       'proper':np.array([ [ 0., 0., 1.],
                                           [-1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                    }
        else:                                                                                       # direct exit for unspecified symmetry
            return  np.ones_like(vector[...,0],bool)

        if proper:
            components_proper   = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(basis['proper'], vector.shape+(3,)),
                                                      vector), 12)
            components_improper = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(basis['improper'], vector.shape+(3,)),
                                                      vector), 12)
            return   np.all(components_proper   >= 0.0,axis=-1) \
                   | np.all(components_improper >= 0.0,axis=-1)
        else:
            components = np.around(np.einsum('...ji,...i',
                                             np.broadcast_to(basis['improper'], vector.shape+(3,)),
                                             np.block([vector[...,:2],np.abs(vector[...,2:3])])), 12)

            return np.all(components >= 0.0,axis=-1)


    def IPF_color(self,vector,in_SST=True,proper=False):
        """
        Map vector to RGB color within standard stereographic triangle of own symmetry.

        Parameters
        ----------
        vector : numpy.ndarray of shape (...,3)
            Vector to colorize.
        in_SST : bool, optional
            Consider symmetrically equivalent orientations such that poles are located in SST.
            Defaults to True.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs (with mirrored colors).
            Defaults to False.

        Returns
        -------
        rgb : numpy.ndarray of shape (...,3)
           RGB array of IPF colors.

        Examples
        --------
        Inverse pole figure color of the e_3 direction for a crystal in "Cube" orientation with cubic symmetry:

        >>> o = damask.Orientation(lattice='cubic')
        >>> o.IPF_color([0,0,1])
        array([1., 0., 0.])

        References
        ----------
        Bases are computed from

        >>> basis = {
        ...    'cubic' :       np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,1.]/np.sqrt(2.),                #              green
        ...                                            [1.,1.,1.]/np.sqrt(3.)]).T),           #              blue
        ...    'hexagonal' :   np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,0.],                            #              green
        ...                                            [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T),  #              blue
        ...    'tetragonal' :  np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,0.],                            #              green
        ...                                            [1.,1.,0.]/np.sqrt(2.)]).T),           #              blue
        ...    'orthorhombic': np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                            [1.,0.,0.],                            #              green
        ...                                            [0.,1.,0.]]).T),                       #              blue
        ...    }

        """
        if np.array(vector).shape[-1] != 3:
            raise ValueError('Input is not a field of three-dimensional vectors.')

        vector_ = self.to_SST(vector,proper) if in_SST else \
                  self @ np.broadcast_to(vector,self.shape+(3,))

        if self.family == 'cubic':
            basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                           [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                           [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                       'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                           [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                           [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                    }
        elif self.family == 'hexagonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -np.sqrt(3.)   ,  0. ],
                                           [ 0.            ,  2.            ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  np.sqrt(3.)   ,  0. ],
                                           [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                    }
        elif self.family == 'tetragonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -1.            ,  0. ],
                                           [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  1.            ,  0. ],
                                           [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                    }
        elif self.family == 'orthorhombic':
            basis = {'improper':np.array([ [ 0., 0., 1.],
                                           [ 1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                       'proper':np.array([ [ 0., 0., 1.],
                                           [-1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                    }
        else:                                                                                       # direct exit for unspecified symmetry
            return np.zeros_like(vector_)

        if proper:
            components_proper   = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(basis['proper'], vector_.shape+(3,)),
                                                      vector_), 12)
            components_improper = np.around(np.einsum('...ji,...i',
                                                      np.broadcast_to(basis['improper'], vector_.shape+(3,)),
                                                      vector_), 12)
            in_SST = np.all(components_proper   >= 0.0,axis=-1) \
                   | np.all(components_improper >= 0.0,axis=-1)
            components = np.where((in_SST & np.all(components_proper   >= 0.0,axis=-1))[...,np.newaxis],
                                  components_proper,components_improper)
        else:
            components = np.around(np.einsum('...ji,...i',
                                             np.broadcast_to(basis['improper'], vector_.shape+(3,)),
                                             np.block([vector_[...,:2],np.abs(vector_[...,2:3])])), 12)

            in_SST = np.all(components >= 0.0,axis=-1)

        with np.errstate(invalid='ignore',divide='ignore'):
            rgb = (components/np.linalg.norm(components,axis=-1,keepdims=True))**0.5            # smoothen color ramps
            rgb = np.clip(rgb,0.,1.)                                                            # clip intensity
            rgb /= np.max(rgb,axis=-1,keepdims=True)                                            # normalize to (HS)V = 1
        rgb[np.broadcast_to(~in_SST[...,np.newaxis],rgb.shape)] = 0.0
        return rgb


    def disorientation(self,other,return_operators=False):
        """
        Calculate disorientation between myself and given other orientation.

        Parameters
        ----------
        other : Orientation
            Orientation to calculate disorientation for.
            Shape of other blends with shape of own rotation array.
            For example, shapes of (2,3) for own rotations and (3,2) for other's result in (2,3,2) disorientations.
        return_operators : bool, optional
            Return index pair of symmetrically equivalent orientations that result in disorientation axis falling into FZ.
            Defaults to False.

        Returns
        -------
        disorientation : Orientation
            Disorientation between self and other.
        operators : numpy.ndarray int of shape (...,2), conditional
            Index of symmetrically equivalent orientation that rotated vector to the SST.

        Notes
        -----
        Currently requires same crystal family for both orientations.
        For extension to cases with differing symmetry see  A. Heinz and P. Neumann 1991 and 10.1107/S0021889808016373.

        Examples
        --------
        Disorientation between two specific orientations of hexagonal symmetry:

        >>> import damask
        >>> a = damask.Orientation.from_Eulers(phi=[123,32,21],degrees=True,lattice='hexagonal')
        >>> b = damask.Orientation.from_Eulers(phi=[104,11,87],degrees=True,lattice='hexagonal')
        >>> a.disorientation(b)
        Crystal family hexagonal
        Quaternion: (real=0.976, imag=<+0.189, +0.018, +0.103>)
        Matrix:
        [[ 0.97831006  0.20710935  0.00389135]
         [-0.19363288  0.90765544  0.37238141]
         [ 0.07359167 -0.36505797  0.92807163]]
        Bunge Eulers / deg: (11.40, 21.86, 0.60)

        """
        if self.family != other.family:
            raise NotImplementedError('disorientation between different crystal families')

        blend = util.shapeblender(self.shape,other.shape)
        s =  self.equivalent
        o = other.equivalent

        s_ = s.reshape((s.shape[0],1)+ self.shape).broadcast_to((s.shape[0],o.shape[0])+blend,mode='right')
        o_ = o.reshape((1,o.shape[0])+other.shape).broadcast_to((s.shape[0],o.shape[0])+blend,mode='right')
        r_ = s_.misorientation(o_)
        _r = ~r_

        forward = r_.in_FZ & r_.in_disorientation_FZ
        reverse = _r.in_FZ & _r.in_disorientation_FZ
        ok  = forward | reverse
        ok &= (np.cumsum(ok.reshape((-1,)+ok.shape[2:]),axis=0) == 1).reshape(ok.shape)
        r = np.where(np.any(forward[...,np.newaxis],axis=(0,1),keepdims=True),
                     r_.quaternion,
                     _r.quaternion)
        loc  = np.where(ok)
        sort = 0 if len(loc) == 2 else np.lexsort(loc[:1:-1])
        quat = r[ok][sort].reshape(blend+(4,))

        return (
                (self.copy(rotation=quat),
                 (np.vstack(loc[:2]).T)[sort].reshape(blend+(2,)))
                if return_operators else
                self.copy(rotation=quat)
               )


    def average(self,weights=None,return_cloud=False):
        """
        Return orientation average over last dimension.

        Parameters
        ----------
        weights : numpy.ndarray, optional
            Relative weights of orientations.
        return_cloud : bool, optional
            Return the set of symmetrically equivalent orientations that was used in averaging.
            Defaults to False.

        Returns
        -------
        average : Orientation
            Weighted average of original Orientation field.
        cloud : Orientations, conditional
            Set of symmetrically equivalent orientations that were used in averaging.

        References
        ----------
        J.C. Glez and J. Driver, Journal of Applied Crystallography 34:280-288, 2001
        https://doi.org/10.1107/S0021889801003077

        """
        eq = self.equivalent
        m  = eq.misorientation(self[...,0].reshape((1,)+self.shape[:-1]+(1,))
                                          .broadcast_to(eq.shape))\
               .as_axis_angle()[...,3]
        r = Rotation(np.squeeze(np.take_along_axis(eq.quaternion,
                                                   np.argmin(m,axis=0)[np.newaxis,...,np.newaxis],
                                                   axis=0),
                                axis=0))
        return (
                (self.copy(rotation=Rotation(r).average(weights)),
                 self.copy(rotation=Rotation(r)))
                if return_cloud else
                self.copy(rotation=Rotation(r).average(weights))
               )


    def to_SST(self,vector,proper=False,return_operators=False):
        """
        Rotate vector to ensure it falls into (improper or proper) standard stereographic triangle of crystal symmetry.

        Parameters
        ----------
        vector : numpy.ndarray of shape (...,3)
            Lab frame vector to align with crystal frame direction.
            Shape of other blends with shape of own rotation array.
            For example, a rotation array of shape (3,2) and a (2,4) vector array result in (3,2,4) outputs.
        proper : bool, optional
            Consider only vectors with z >= 0, hence combine two neighboring SSTs.
            Defaults to False.
        return_operators : bool, optional
            Return the symmetrically equivalent orientation that rotated vector to SST.
            Defaults to False.

        Returns
        -------
        vector_SST : numpy.ndarray of shape (...,3)
            Rotated vector falling into SST.
        operators : numpy.ndarray int of shape (...), conditional
            Index of symmetrically equivalent orientation that rotated vector to SST.

        """
        eq  = self.equivalent
        blend = util.shapeblender(eq.shape,np.array(vector).shape[:-1])
        poles = eq.broadcast_to(blend,mode='right') @ np.broadcast_to(np.array(vector),blend+(3,))
        ok    = self.in_SST(poles,proper=proper)
        ok   &= np.cumsum(ok,axis=0) == 1
        loc   = np.where(ok)
        sort  = 0 if len(loc) == 1 else np.lexsort(loc[:0:-1])
        return (
                (poles[ok][sort].reshape(blend[1:]+(3,)), (np.vstack(loc[:1]).T)[sort].reshape(blend[1:]))
                if return_operators else
                poles[ok][sort].reshape(blend[1:]+(3,))
               )


    @classmethod
    def Bravais_to_Miller(cls,*,uvtw=None,hkil=None):
        """
        Transform 4 Miller–Bravais indices to 3 Miller indices of crystal direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        uvtw|hkil : numpy.ndarray of shape (...,4)
            Miller–Bravais indices of crystallographic direction [uvtw] or plane normal (hkil).

        Returns
        -------
        uvw|hkl : numpy.ndarray of shape (...,3)
            Miller indices of [uvw] direction or (hkl) plane normal.

        """
        if (uvtw is not None) ^ (hkil is None):
            raise KeyError('Specify either "uvtw" or "hkil"')
        axis,basis  = (np.array(uvtw),np.array([[1,0,-1,0],
                                                [0,1,-1,0],
                                                [0,0, 0,1]])) \
                      if hkil is None else \
                      (np.array(hkil),np.array([[1,0,0,0],
                                                [0,1,0,0],
                                                [0,0,0,1]]))
        return np.einsum('il,...l',basis,axis)


    @classmethod
    def Miller_to_Bravais(cls,*,uvw=None,hkl=None):
        """
        Transform 3 Miller indices to 4 Miller–Bravais indices of crystal direction [uvtw] or plane normal (hkil).

        Parameters
        ----------
        uvw|hkl : numpy.ndarray of shape (...,3)
            Miller indices of crystallographic direction [uvw] or plane normal (hkl).

        Returns
        -------
        uvtw|hkil : numpy.ndarray of shape (...,4)
            Miller–Bravais indices of [uvtw] direction or (hkil) plane normal.

        """
        if (uvw is not None) ^ (hkl is None):
            raise KeyError('Specify either "uvw" or "hkl"')
        axis,basis  = (np.array(uvw),np.array([[ 2,-1, 0],
                                               [-1, 2, 0],
                                               [-1,-1, 0],
                                               [ 0, 0, 3]])/3) \
                      if hkl is None else \
                      (np.array(hkl),np.array([[ 1, 0, 0],
                                               [ 0, 1, 0],
                                               [-1,-1, 0],
                                               [ 0, 0, 1]]))
        return np.einsum('il,...l',basis,axis)


    def to_lattice(self,*,direction=None,plane=None):
        """
        Calculate lattice vector corresponding to crystal frame direction or plane normal.

        Parameters
        ----------
        direction|normal : numpy.ndarray of shape (...,3)
            Vector along direction or plane normal.

        Returns
        -------
        Miller : numpy.ndarray of shape (...,3)
            lattice vector of direction or plane.
            Use util.scale_to_coprime to convert to (integer) Miller indices.

        """
        if (direction is not None) ^ (plane is None):
            raise KeyError('Specify either "direction" or "plane"')
        axis,basis  = (np.array(direction),self.basis_reciprocal.T) \
                      if plane is None else \
                      (np.array(plane),self.basis_real.T)
        return np.einsum('il,...l',basis,axis)


    def to_frame(self,*,uvw=None,hkl=None,with_symmetry=False):
        """
        Calculate crystal frame vector along lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        uvw|hkl : numpy.ndarray of shape (...,3)
            Miller indices of crystallographic direction or plane normal.
        with_symmetry : bool, optional
            Calculate all N symmetrically equivalent vectors.

        Returns
        -------
        vector : numpy.ndarray of shape (...,3) or (N,...,3)
            Crystal frame vector (or vectors if with_symmetry) along [uvw] direction or (hkl) plane normal.

        """
        if (uvw is not None) ^ (hkl is None):
            raise KeyError('Specify either "uvw" or "hkl"')
        axis,basis  = (np.array(uvw),self.basis_real) \
                      if hkl is None else \
                      (np.array(hkl),self.basis_reciprocal)
        return (self.symmetry_operations.broadcast_to(self.symmetry_operations.shape+axis.shape[:-1],mode='right')
              @ np.broadcast_to(np.einsum('il,...l',basis,axis),self.symmetry_operations.shape+axis.shape)
                if with_symmetry else
                np.einsum('il,...l',basis,axis))


    def to_pole(self,*,uvw=None,hkl=None,with_symmetry=False):
        """
        Calculate lab frame vector along lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        uvw|hkl : numpy.ndarray of shape (...,3)
            Miller indices of crystallographic direction or plane normal.
        with_symmetry : bool, optional
            Calculate all N symmetrically equivalent vectors.

        Returns
        -------
        vector : numpy.ndarray of shape (...,3) or (N,...,3)
            Lab frame vector (or vectors if with_symmetry) along [uvw] direction or (hkl) plane normal.

        """
        v = self.to_frame(uvw=uvw,hkl=hkl,with_symmetry=with_symmetry)
        return ~(self if self.shape+v.shape[:-1] == () else self.broadcast_to(self.shape+v.shape[:-1],mode='right')) \
               @ np.broadcast_to(v,self.shape+v.shape)


    def Schmid(self,mode):
        u"""
        Calculate Schmid matrix P = d ⨂ n in the lab frame for given lattice shear kinematics.

        Parameters
        ----------
        mode : str
            Type of kinematics, i.e. 'slip' or 'twin'.

        Returns
        -------
        P : numpy.ndarray of shape (...,N,3,3)
            Schmid matrix for each of the N deformation systems.

        Examples
        --------
        Schmid matrix (in lab frame) of slip systems of a face-centered
        cubic crystal in "Goss" orientation.

        >>> import damask
        >>> import numpy as np
        >>> np.set_printoptions(3,suppress=True,floatmode='fixed')
        >>> damask.Orientation.from_Eulers(phi=[0,45,0],degrees=True,lattice='cF').Schmid('slip')[0]
        array([[ 0.000,  0.000,  0.000],
               [ 0.577, -0.000,  0.816],
               [ 0.000,  0.000,  0.000]])

        """
        d = self.to_frame(uvw=self.kinematics[mode]['direction'],with_symmetry=False)
        p = self.to_frame(hkl=self.kinematics[mode]['plane']    ,with_symmetry=False)
        P = np.einsum('...i,...j',d/np.linalg.norm(d,axis=-1,keepdims=True),
                                  p/np.linalg.norm(p,axis=-1,keepdims=True))

        return ~self.broadcast_to( self.shape+P.shape[:-2],mode='right') \
               @ np.broadcast_to(P,self.shape+P.shape)
