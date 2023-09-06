from typing import Optional, Union, Dict, List, Tuple

import numpy as np

from ._typehints import FloatSequence, CrystalFamily, BravaisLattice, CrystalKinematics
from . import util
from . import Rotation

lattice_symmetries: Dict[BravaisLattice, CrystalFamily] = {
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

orientation_relationships: Dict[str, Dict[BravaisLattice,np.ndarray]] = {
 'KS': {
   'cF': np.array([
       [[-1, 0, 1],[ 1, 1, 1]],
       [[-1, 0, 1],[ 1, 1, 1]],
       [[ 0, 1,-1],[ 1, 1, 1]],
       [[ 0, 1,-1],[ 1, 1, 1]],
       [[ 1,-1, 0],[ 1, 1, 1]],
       [[ 1,-1, 0],[ 1, 1, 1]],
       [[ 1, 0,-1],[ 1,-1, 1]],
       [[ 1, 0,-1],[ 1,-1, 1]],
       [[-1,-1, 0],[ 1,-1, 1]],
       [[-1,-1, 0],[ 1,-1, 1]],
       [[ 0, 1, 1],[ 1,-1, 1]],
       [[ 0, 1, 1],[ 1,-1, 1]],
       [[ 0,-1, 1],[-1, 1, 1]],
       [[ 0,-1, 1],[-1, 1, 1]],
       [[-1, 0,-1],[-1, 1, 1]],
       [[-1, 0,-1],[-1, 1, 1]],
       [[ 1, 1, 0],[-1, 1, 1]],
       [[ 1, 1, 0],[-1, 1, 1]],
       [[-1, 1, 0],[ 1, 1,-1]],
       [[-1, 1, 0],[ 1, 1,-1]],
       [[ 0,-1,-1],[ 1, 1,-1]],
       [[ 0,-1,-1],[ 1, 1,-1]],
       [[ 1, 0, 1],[ 1, 1,-1]],
       [[ 1, 0, 1],[ 1, 1,-1]],
       ],dtype=float),
   'cI': np.array([
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       [[-1,-1, 1],[ 0, 1, 1]],
       [[-1, 1,-1],[ 0, 1, 1]],
       ],dtype=float),
 },
 'GT': {
   'cF': np.array([
       [[ -5,-12, 17],[  1,  1,  1]],
       [[ 17, -5,-12],[  1,  1,  1]],
       [[-12, 17, -5],[  1,  1,  1]],
       [[  5, 12, 17],[ -1, -1,  1]],
       [[-17,  5,-12],[ -1, -1,  1]],
       [[ 12,-17, -5],[ -1, -1,  1]],
       [[ -5, 12,-17],[ -1,  1,  1]],
       [[ 17,  5, 12],[ -1,  1,  1]],
       [[-12,-17,  5],[ -1,  1,  1]],
       [[  5,-12,-17],[  1, -1,  1]],
       [[-17, -5, 12],[  1, -1,  1]],
       [[ 12, 17,  5],[  1, -1,  1]],
       [[ -5, 17,-12],[  1,  1,  1]],
       [[-12, -5, 17],[  1,  1,  1]],
       [[ 17,-12, -5],[  1,  1,  1]],
       [[  5,-17,-12],[ -1, -1,  1]],
       [[ 12,  5, 17],[ -1, -1,  1]],
       [[-17, 12, -5],[ -1, -1,  1]],
       [[ -5,-17, 12],[ -1,  1,  1]],
       [[-12,  5,-17],[ -1,  1,  1]],
       [[ 17, 12,  5],[ -1,  1,  1]],
       [[  5, 17, 12],[  1, -1,  1]],
       [[ 12, -5,-17],[  1, -1,  1]],
       [[-17,-12,  5],[  1, -1,  1]],
       ],dtype=float),
   'cI': np.array([
       [[-17, -7, 17],[  1,  0,  1]],
       [[ 17,-17, -7],[  1,  1,  0]],
       [[ -7, 17,-17],[  0,  1,  1]],
       [[ 17,  7, 17],[ -1,  0,  1]],
       [[-17, 17, -7],[ -1, -1,  0]],
       [[  7,-17,-17],[  0, -1,  1]],
       [[-17,  7,-17],[ -1,  0,  1]],
       [[ 17, 17,  7],[ -1,  1,  0]],
       [[ -7,-17, 17],[  0,  1,  1]],
       [[ 17, -7,-17],[  1,  0,  1]],
       [[-17,-17,  7],[  1, -1,  0]],
       [[  7, 17, 17],[  0, -1,  1]],
       [[-17, 17, -7],[  1,  1,  0]],
       [[ -7,-17, 17],[  0,  1,  1]],
       [[ 17, -7,-17],[  1,  0,  1]],
       [[ 17,-17, -7],[ -1, -1,  0]],
       [[  7, 17, 17],[  0, -1,  1]],
       [[-17,  7,-17],[ -1,  0,  1]],
       [[-17,-17,  7],[ -1,  1,  0]],
       [[ -7, 17,-17],[  0,  1,  1]],
       [[ 17,  7, 17],[ -1,  0,  1]],
       [[ 17, 17,  7],[  1, -1,  0]],
       [[  7,-17,-17],[  0, -1,  1]],
       [[-17, -7, 17],[  1,  0,  1]],
       ],dtype=float),
 },
 'GT_prime': {
   'cF' : np.array([
       [[  0,  1, -1],[  7, 17, 17]],
       [[ -1,  0,  1],[ 17,  7, 17]],
       [[  1, -1,  0],[ 17, 17,  7]],
       [[  0, -1, -1],[ -7,-17, 17]],
       [[  1,  0,  1],[-17, -7, 17]],
       [[  1, -1,  0],[-17,-17,  7]],
       [[  0,  1, -1],[  7,-17,-17]],
       [[  1,  0,  1],[ 17, -7,-17]],
       [[ -1, -1,  0],[ 17,-17, -7]],
       [[  0, -1, -1],[ -7, 17,-17]],
       [[ -1,  0,  1],[-17,  7,-17]],
       [[ -1, -1,  0],[-17, 17, -7]],
       [[  0, -1,  1],[  7, 17, 17]],
       [[  1,  0, -1],[ 17,  7, 17]],
       [[ -1,  1,  0],[ 17, 17,  7]],
       [[  0,  1,  1],[ -7,-17, 17]],
       [[ -1,  0, -1],[-17, -7, 17]],
       [[ -1,  1,  0],[-17,-17,  7]],
       [[  0, -1,  1],[  7,-17,-17]],
       [[ -1,  0, -1],[ 17, -7,-17]],
       [[  1,  1,  0],[ 17,-17, -7]],
       [[  0,  1,  1],[ -7, 17,-17]],
       [[  1,  0, -1],[-17,  7,-17]],
       [[  1,  1,  0],[-17, 17, -7]],
       ],dtype=float),
   'cI' : np.array([
       [[  1,  1, -1],[ 12,  5, 17]],
       [[ -1,  1,  1],[ 17, 12,  5]],
       [[  1, -1,  1],[  5, 17, 12]],
       [[ -1, -1, -1],[-12, -5, 17]],
       [[  1, -1,  1],[-17,-12,  5]],
       [[  1, -1, -1],[ -5,-17, 12]],
       [[ -1,  1, -1],[ 12, -5,-17]],
       [[  1,  1,  1],[ 17,-12, -5]],
       [[ -1, -1,  1],[  5,-17,-12]],
       [[  1, -1, -1],[-12,  5,-17]],
       [[ -1, -1,  1],[-17, 12, -5]],
       [[ -1, -1, -1],[ -5, 17,-12]],
       [[  1, -1,  1],[ 12, 17,  5]],
       [[  1,  1, -1],[  5, 12, 17]],
       [[ -1,  1,  1],[ 17,  5, 12]],
       [[ -1,  1,  1],[-12,-17,  5]],
       [[ -1, -1, -1],[ -5,-12, 17]],
       [[ -1,  1, -1],[-17, -5, 12]],
       [[ -1, -1,  1],[ 12,-17, -5]],
       [[ -1,  1, -1],[  5,-12,-17]],
       [[  1,  1,  1],[ 17, -5,-12]],
       [[  1,  1,  1],[-12, 17, -5]],
       [[  1, -1, -1],[ -5, 12,-17]],
       [[  1,  1, -1],[-17,  5,-12]],
       ],dtype=float),
 },
 'NW': {
   'cF' : np.array([
       [[  2, -1, -1],[  1,  1,  1]],
       [[ -1,  2, -1],[  1,  1,  1]],
       [[ -1, -1,  2],[  1,  1,  1]],
       [[ -2, -1, -1],[ -1,  1,  1]],
       [[  1,  2, -1],[ -1,  1,  1]],
       [[  1, -1,  2],[ -1,  1,  1]],
       [[  2,  1, -1],[  1, -1,  1]],
       [[ -1, -2, -1],[  1, -1,  1]],
       [[ -1,  1,  2],[  1, -1,  1]],
       [[  2, -1,  1],[ -1, -1,  1]],
       [[ -1,  2,  1],[ -1, -1,  1]],
       [[ -1, -1, -2],[ -1, -1,  1]],
       ],dtype=float),
   'cI' : np.array([
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       [[  0, -1,  1],[  0,  1,  1]],
       ],dtype=float),
 },
 'Pitsch': {
   'cF' : np.array([
       [[  1,  0,  1],[  0,  1,  0]],
       [[  1,  1,  0],[  0,  0,  1]],
       [[  0,  1,  1],[  1,  0,  0]],
       [[  0,  1, -1],[  1,  0,  0]],
       [[ -1,  0,  1],[  0,  1,  0]],
       [[  1, -1,  0],[  0,  0,  1]],
       [[  1,  0, -1],[  0,  1,  0]],
       [[ -1,  1,  0],[  0,  0,  1]],
       [[  0, -1,  1],[  1,  0,  0]],
       [[  0,  1,  1],[  1,  0,  0]],
       [[  1,  0,  1],[  0,  1,  0]],
       [[  1,  1,  0],[  0,  0,  1]],
       ],dtype=float),
   'cI' : np.array([
       [[  1, -1,  1],[ -1,  0,  1]],
       [[  1,  1, -1],[  1, -1,  0]],
       [[ -1,  1,  1],[  0,  1, -1]],
       [[ -1,  1, -1],[  0, -1, -1]],
       [[ -1, -1,  1],[ -1,  0, -1]],
       [[  1, -1, -1],[ -1, -1,  0]],
       [[  1, -1, -1],[ -1,  0, -1]],
       [[ -1,  1, -1],[ -1, -1,  0]],
       [[ -1, -1,  1],[  0, -1, -1]],
       [[ -1,  1,  1],[  0, -1,  1]],
       [[  1, -1,  1],[  1,  0, -1]],
       [[  1,  1, -1],[ -1,  1,  0]],
       ],dtype=float),
 },
 'Bain': {
   'cF' : np.array([
       [[  0,  1,  0],[  1,  0,  0]],
       [[  0,  0,  1],[  0,  1,  0]],
       [[  1,  0,  0],[  0,  0,  1]],
       ],dtype=float),
   'cI' : np.array([
       [[  0,  1,  1],[  1,  0,  0]],
       [[  1,  0,  1],[  0,  1,  0]],
       [[  1,  1,  0],[  0,  0,  1]],
       ],dtype=float),
 },
 'Burgers' : {
   'cI' : np.array([
       [[ -1,  1,  1],[  1,  1,  0]],
       [[ -1,  1, -1],[  1,  1,  0]],
       [[  1,  1,  1],[  1, -1,  0]],
       [[  1,  1, -1],[  1, -1,  0]],

       [[  1,  1, -1],[  1,  0,  1]],
       [[ -1,  1,  1],[  1,  0,  1]],
       [[  1,  1,  1],[ -1,  0,  1]],
       [[  1, -1,  1],[ -1,  0,  1]],

       [[ -1,  1, -1],[  0,  1,  1]],
       [[  1,  1, -1],[  0,  1,  1]],
       [[ -1,  1,  1],[  0, -1,  1]],
       [[  1,  1,  1],[  0, -1,  1]],
     ],dtype=float),
   'hP' : np.array([
       [[  -1,  2,  -1, 0],[  0,  0,  0,  1]],
       [[  -1, -1,   2, 0],[  0,  0,  0,  1]],
       [[  -1,  2,  -1, 0],[  0,  0,  0,  1]],
       [[  -1, -1,   2, 0],[  0,  0,  0,  1]],

       [[  -1,  2,  -1, 0],[  0,  0,  0,  1]],
       [[  -1, -1,   2, 0],[  0,  0,  0,  1]],
       [[  -1,  2,  -1, 0],[  0,  0,  0,  1]],
       [[  -1, -1,   2, 0],[  0,  0,  0,  1]],

       [[  -1,  2,  -1, 0],[  0,  0,  0,  1]],
       [[  -1, -1,   2, 0],[  0,  0,  0,  1]],
       [[  -1,  2,  -1, 0],[  0,  0,  0,  1]],
       [[  -1, -1,   2, 0],[  0,  0,  0,  1]],
     ],dtype=float),
 },
 }

class Crystal():
    """
    Representation of a crystal as (general) crystal family or (more specific) as a scaled Bravais lattice.

    Examples
    --------
    Cubic crystal family:

    >>> import damask
    >>> (cubic := damask.Crystal(family='cubic'))
    Crystal family: cubic

    Body-centered cubic Bravais lattice with parameters of iron:

    >>> import damask
    >>> (Fe := damask.Crystal(lattice='cI', a=287e-12))
    Crystal family: cubic
    Bravais lattice: cI
    a=2.87e-10 m, b=2.87e-10 m, c=2.87e-10 m
    α=90°, β=90°, γ=90°

    """

    def __init__(self, *,
                 family: Optional[CrystalFamily] = None,
                 lattice: Optional[BravaisLattice] = None,
                 a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                 alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                 degrees: bool = False):
        """
        New representation of a crystal.

        Parameters
        ----------
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional.
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional.
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
        if family is not None and family not in list(lattice_symmetries.values()):
            raise KeyError(f'invalid crystal family "{family}"')
        if lattice is not None and family is not None and family != lattice_symmetries[lattice]:
            raise KeyError(f'incompatible family "{family}" for lattice "{lattice}"')

        self.family = lattice_symmetries[lattice] if family is None else family
        self.lattice = lattice

        if self.lattice is not None:
            self.a = 1 if a is None else a
            self.b = b
            self.c = c
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

            self.alpha = np.radians(alpha) if degrees and alpha is not None else alpha
            self.beta  = np.radians(beta)  if degrees and beta  is not None else beta
            self.gamma = np.radians(gamma) if degrees and gamma is not None else gamma
            if self.alpha is None and 'alpha' in self.immutable: self.alpha = self.immutable['alpha']
            if self.beta  is None and 'beta'  in self.immutable: self.beta  = self.immutable['beta']
            if self.gamma is None and 'gamma' in self.immutable: self.gamma = self.immutable['gamma']

            if \
                (self.a     is None) \
             or (self.b     is None or ('b'     in self.immutable and self.b     != self.immutable['b'] * self.a)) \
             or (self.c     is None or ('c'     in self.immutable and self.c     != self.immutable['c'] * self.b)) \
             or (self.alpha is None or ('alpha' in self.immutable and self.alpha != self.immutable['alpha'])) \
             or (self.beta  is None or ('beta'  in self.immutable and self.beta  != self.immutable['beta'])) \
             or (self.gamma is None or ('gamma' in self.immutable and self.gamma != self.immutable['gamma'])):
                raise ValueError (f'incompatible parameters {self.parameters} for crystal family {self.family}')

            if np.any(np.array([self.alpha,self.beta,self.gamma]) <= 0):
                raise ValueError ('lattice angles must be positive')
            if np.any([np.roll([self.alpha,self.beta,self.gamma],r)[0]
              >= np.sum(np.roll([self.alpha,self.beta,self.gamma],r)[1:]) for r in range(3)]):
                raise ValueError ('each lattice angle must be less than sum of others')


    def __repr__(self):
        """
        Return repr(self).

        Give short, human-readable summary.

        """
        family = f'Crystal family: {self.family}'
        return family if self.lattice is None else \
               util.srepr([family,
                           f'Bravais lattice: {self.lattice}',
                           'a={:.5g} m, b={:.5g} m, c={:.5g} m'.format(*self.parameters[:3]),
                           'α={:.5g}°, β={:.5g}°, γ={:.5g}°'.format(*np.degrees(self.parameters[3:]))])


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : Crystal
            Crystal to check for equality.

        """
        return NotImplemented if not isinstance(other, Crystal) else \
               self.lattice == other.lattice and \
               self.parameters == other.parameters and \
               self.family == other.family

    @property
    def parameters(self):
        """Return lattice parameters a, b, c, alpha, beta, gamma."""
        if hasattr(self,'a'): return (self.a,self.b,self.c,self.alpha,self.beta,self.gamma)

    @property
    def immutable(self):
        """Return immutable lattice parameters."""
        _immutable = {
            'cubic': {
                         'b': 1.0,
                         'c': 1.0,
                         'alpha': np.pi/2.,
                         'beta':  np.pi/2.,
                         'gamma': np.pi/2.,
                       },
            'hexagonal': {
                         'b': 1.0,
                         'alpha': np.pi/2.,
                         'beta':  np.pi/2.,
                         'gamma': 2.*np.pi/3.,
                       },
            'tetragonal': {
                         'b': 1.0,
                         'alpha': np.pi/2.,
                         'beta':  np.pi/2.,
                         'gamma': np.pi/2.,
                       },
            'orthorhombic': {
                         'alpha': np.pi/2.,
                         'beta':  np.pi/2.,
                         'gamma': np.pi/2.,
                       },
            'monoclinic': {
                         'alpha': np.pi/2.,
                         'gamma': np.pi/2.,
                       },
            'triclinic': {}
                     }
        return _immutable[self.family]


    @property
    def orientation_relationships(self):
        """Return labels of orientation relationships."""
        return [k for k,v in orientation_relationships.items() if self.lattice in v]


    @property
    def standard_triangle(self) -> Union[Dict[str, np.ndarray], None]:
        """
        Corners of the standard triangle.

        Notes
        -----
        Not yet defined for monoclinic.


        References
        ----------
        Bases are computed from

        >>> basis = {
        ...    'cubic' :       np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,1.]/np.sqrt(2.),               #              green
        ...                                            [1.,1.,1.]/np.sqrt(3.)]).T),          #              blue
        ...    'hexagonal' :   np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,0.],                           #              green
        ...                                            [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T), #              blue
        ...    'tetragonal' :  np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,0.],                           #              green
        ...                                            [1.,1.,0.]/np.sqrt(2.)]).T),          #              blue
        ...    'orthorhombic': np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,0.],                           #              green
        ...                                            [0.,1.,0.]]).T),                      #              blue
        ...    }

        """
        _basis: Dict[CrystalFamily, Dict[str, np.ndarray]]  = {
            'cubic':    {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                               [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                               [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                           'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                               [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                               [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                        },
            'hexagonal':
                        {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                               [ 1.            , -np.sqrt(3.)   ,  0. ],
                                               [ 0.            ,  2.            ,  0. ] ]),
                           'proper':np.array([ [ 0.            ,  0.            ,  1. ],
                                               [-1.            ,  np.sqrt(3.)   ,  0. ],
                                               [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                        },
            'tetragonal':
                        {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                               [ 1.            , -1.            ,  0. ],
                                               [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                           'proper':np.array([ [ 0.            ,  0.            ,  1. ],
                                               [-1.            ,  1.            ,  0. ],
                                               [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                        },
            'orthorhombic':
                        {'improper':np.array([ [ 0., 0., 1.],
                                               [ 1., 0., 0.],
                                               [ 0., 1., 0.] ]),
                           'proper':np.array([ [ 0., 0., 1.],
                                               [-1., 0., 0.],
                                               [ 0., 1., 0.] ]),
                        }}
        return _basis.get(self.family, None)


    @property
    def symmetry_operations(self) -> Rotation:
        """
        Return symmetry operations.

        References
        ----------
        U.F. Kocks et al.,
        Texture and Anisotropy:
        Preferred Orientations in Polycrystals and their Effect on Materials Properties.
        Cambridge University Press 1998. Table II

        """
        _symmetry_operations: Dict[CrystalFamily, List]  = {
            'cubic':         [
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
                            ],
            'hexagonal':    [
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
                            ],
            'tetragonal':   [
                              [ 1.0,            0.0,            0.0,            0.0            ],
                              [ 0.0,            1.0,            0.0,            0.0            ],
                              [ 0.0,            0.0,            1.0,            0.0            ],
                              [ 0.0,            0.0,            0.0,            1.0            ],
                              [ 0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                              [ 0.0,           -0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                              [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                              [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                            ],
            'orthorhombic': [
                              [ 1.0,0.0,0.0,0.0 ],
                              [ 0.0,1.0,0.0,0.0 ],
                              [ 0.0,0.0,1.0,0.0 ],
                              [ 0.0,0.0,0.0,1.0 ],
                            ],
            'monoclinic':   [
                              [ 1.0,0.0,0.0,0.0 ],
                              [ 0.0,0.0,1.0,0.0 ],
                            ],
            'triclinic':    [
                              [ 1.0,0.0,0.0,0.0 ],
                            ]}
        return Rotation.from_quaternion(_symmetry_operations[self.family],accept_homomorph=True)


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
    def basis_real(self) -> np.ndarray:
        """
        Return orthogonal real space crystal basis.

        References
        ----------
        C.T. Young and J.L. Lytton, Journal of Applied Physics 43:1408–1417, 1972
        https://doi.org/10.1063/1.1661333

        """
        if self.parameters is None:
            raise KeyError('missing crystal lattice parameters')
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
    def basis_reciprocal(self) -> np.ndarray:
        """Return reciprocal (dual) crystal basis."""
        return np.linalg.inv(self.basis_real.T)


    @property
    def lattice_points(self):
        """Return lattice points."""
        _lattice_points = {
                'P': [
                     ],
                'S': [
                      [0.5,0.5,0],
                     ],
                'I': [
                      [0.5,0.5,0.5],
                     ],
                'F': [
                      [0.0,0.5,0.5],
                      [0.5,0.0,0.5],
                      [0.5,0.5,0.0],
                     ],
                'hP': [
                       [2./3.,1./3.,0.5],
                     ],
                }

        if self.lattice is None: raise KeyError('no lattice type specified')
        return np.array([[0,0,0]]
                        + _lattice_points.get(self.lattice if self.lattice == 'hP' else \
                                              self.lattice[-1],None),dtype=float)

    def to_lattice(self, *,
                   direction: Optional[FloatSequence] = None,
                   plane: Optional[FloatSequence] = None) -> np.ndarray:
        """
        Calculate lattice vector corresponding to crystal frame direction or plane normal.

        Parameters
        ----------
        direction|plane : numpy.ndarray, shape (...,3)
            Real space vector along direction or
            reciprocal space vector along plane normal.

        Returns
        -------
        Miller : numpy.ndarray, shape (...,3)
            Lattice vector of direction or plane.
            Use util.scale_to_coprime to convert to (integer) Miller indices.

        """
        if (direction is not None) ^ (plane is None):
            raise KeyError('specify either "direction" or "plane"')
        basis,axis = (self.basis_reciprocal,np.array(direction)) \
                     if plane is None else \
                     (self.basis_real,np.array(plane))
        return np.einsum('li,...l',basis,axis)


    def to_frame(self, *,
                 uvw: Optional[FloatSequence] = None,
                 hkl: Optional[FloatSequence] = None) -> np.ndarray:
        """
        Calculate crystal frame vector corresponding to lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        uvw|hkl : numpy.ndarray, shape (...,3)
            Miller indices of crystallographic direction or plane normal.

        Returns
        -------
        vector : numpy.ndarray, shape (...,3)
            Crystal frame vector in real space along [uvw] direction or
            in reciprocal space along (hkl) plane normal.

        Examples
        --------
        Crystal frame vector (real space) of Magnesium corresponding to [1,1,0] direction:

        >>> import damask
        >>> Mg = damask.Crystal(lattice='hP', a=321e-12, c=521e-12)
        >>> Mg.to_frame(uvw=[1, 1, 0])
        array([1.60500000e-10, 2.77994155e-10, 0.00000000e+00])

        Crystal frame vector (reciprocal space) of Titanium along (1,0,0) plane normal:

        >>> import damask
        >>> Ti = damask.Crystal(lattice='hP', a=0.295e-9, c=0.469e-9)
        >>> Ti.to_frame(hkl=(1, 0, 0))
        array([ 3.38983051e+09,  1.95711956e+09, -4.15134508e-07])

        """
        if (uvw is not None) ^ (hkl is None):
            raise KeyError('specify either "uvw" or "hkl"')
        basis,axis = (self.basis_real,np.array(uvw)) \
                     if hkl is None else \
                     (self.basis_reciprocal,np.array(hkl))
        return np.einsum('il,...l',basis,axis)


    def kinematics(self,
                   mode: CrystalKinematics) -> Dict[str, List[np.ndarray]]:
        """
        Return crystal kinematics systems.

        Parameters
        ----------
        mode : {'slip','twin'}
            Deformation mode.

        Returns
        -------
        direction_plane : dictionary
            Directions and planes of deformation mode families.

        """
        _kinematics: Dict[BravaisLattice, Dict[CrystalKinematics, List[np.ndarray]]] = {
            'cF': {
                'slip': [np.array([
                           [+0,+1,-1, +1,+1,+1],
                           [-1,+0,+1, +1,+1,+1],
                           [+1,-1,+0, +1,+1,+1],
                           [+0,-1,-1, -1,-1,+1],
                           [+1,+0,+1, -1,-1,+1],
                           [-1,+1,+0, -1,-1,+1],
                           [+0,-1,+1, +1,-1,-1],
                           [-1,+0,-1, +1,-1,-1],
                           [+1,+1,+0, +1,-1,-1],
                           [+0,+1,+1, -1,+1,-1],
                           [+1,+0,-1, -1,+1,-1],
                           [-1,-1,+0, -1,+1,-1]]),
                         np.array([
                           [+1,+1,+0, +1,-1,+0],
                           [+1,-1,+0, +1,+1,+0],
                           [+1,+0,+1, +1,+0,-1],
                           [+1,+0,-1, +1,+0,+1],
                           [+0,+1,+1, +0,+1,-1],
                           [+0,+1,-1, +0,+1,+1]])],
                'twin': [np.array([
                           [-2, 1, 1,  1, 1, 1],
                           [ 1,-2, 1,  1, 1, 1],
                           [ 1, 1,-2,  1, 1, 1],
                           [ 2,-1, 1, -1,-1, 1],
                           [-1, 2, 1, -1,-1, 1],
                           [-1,-1,-2, -1,-1, 1],
                           [-2,-1,-1,  1,-1,-1],
                           [ 1, 2,-1,  1,-1,-1],
                           [ 1,-1, 2,  1,-1,-1],
                           [ 2, 1,-1, -1, 1,-1],
                           [-1,-2,-1, -1, 1,-1],
                           [-1, 1, 2, -1, 1,-1]])]
            },
            'cI': {
                'slip': [np.array([
                           [+1,-1,+1, +0,+1,+1],
                           [-1,-1,+1, +0,+1,+1],
                           [+1,+1,+1, +0,-1,+1],
                           [-1,+1,+1, +0,-1,+1],
                           [-1,+1,+1, +1,+0,+1],
                           [-1,-1,+1, +1,+0,+1],
                           [+1,+1,+1, -1,+0,+1],
                           [+1,-1,+1, -1,+0,+1],
                           [-1,+1,+1, +1,+1,+0],
                           [-1,+1,-1, +1,+1,+0],
                           [+1,+1,+1, -1,+1,+0],
                           [+1,+1,-1, -1,+1,+0]]),
                         np.array([
                           [-1,+1,+1, +2,+1,+1],
                           [+1,+1,+1, -2,+1,+1],
                           [+1,+1,-1, +2,-1,+1],
                           [+1,-1,+1, +2,+1,-1],
                           [+1,-1,+1, +1,+2,+1],
                           [+1,+1,-1, -1,+2,+1],
                           [+1,+1,+1, +1,-2,+1],
                           [-1,+1,+1, +1,+2,-1],
                           [+1,+1,-1, +1,+1,+2],
                           [+1,-1,+1, -1,+1,+2],
                           [-1,+1,+1, +1,-1,+2],
                           [+1,+1,+1, +1,+1,-2]]),
                         np.array([
                           [+1,+1,-1, +1,+2,+3],
                           [+1,-1,+1, -1,+2,+3],
                           [-1,+1,+1, +1,-2,+3],
                           [+1,+1,+1, +1,+2,-3],
                           [+1,-1,+1, +1,+3,+2],
                           [+1,+1,-1, -1,+3,+2],
                           [+1,+1,+1, +1,-3,+2],
                           [-1,+1,+1, +1,+3,-2],
                           [+1,+1,-1, +2,+1,+3],
                           [+1,-1,+1, -2,+1,+3],
                           [-1,+1,+1, +2,-1,+3],
                           [+1,+1,+1, +2,+1,-3],
                           [+1,-1,+1, +2,+3,+1],
                           [+1,+1,-1, -2,+3,+1],
                           [+1,+1,+1, +2,-3,+1],
                           [-1,+1,+1, +2,+3,-1],
                           [-1,+1,+1, +3,+1,+2],
                           [+1,+1,+1, -3,+1,+2],
                           [+1,+1,-1, +3,-1,+2],
                           [+1,-1,+1, +3,+1,-2],
                           [-1,+1,+1, +3,+2,+1],
                           [+1,+1,+1, -3,+2,+1],
                           [+1,+1,-1, +3,-2,+1],
                           [+1,-1,+1, +3,+2,-1]])],
                'twin': [np.array([
                           [-1, 1, 1,  2, 1, 1],
                           [ 1, 1, 1, -2, 1, 1],
                           [ 1, 1,-1,  2,-1, 1],
                           [ 1,-1, 1,  2, 1,-1],
                           [ 1,-1, 1,  1, 2, 1],
                           [ 1, 1,-1, -1, 2, 1],
                           [ 1, 1, 1,  1,-2, 1],
                           [-1, 1, 1,  1, 2,-1],
                           [ 1, 1,-1,  1, 1, 2],
                           [ 1,-1, 1, -1, 1, 2],
                           [-1, 1, 1,  1,-1, 2],
                           [ 1, 1, 1,  1, 1,-2]])]
            },
            'hP': {
                'slip': [np.array([
                           [+2,-1,-1,+0, +0,+0,+0,+1],
                           [-1,+2,-1,+0, +0,+0,+0,+1],
                           [-1,-1,+2,+0, +0,+0,+0,+1]]),
                         np.array([
                           [+2,-1,-1,+0, +0,+1,-1,+0],
                           [-1,+2,-1,+0, -1,+0,+1,+0],
                           [-1,-1,+2,+0, +1,-1,+0,+0]]),
                         np.array([
                           [-1,+2,-1,+0, +1,+0,-1,+1],
                           [-2,+1,+1,+0, +0,+1,-1,+1],
                           [-1,-1,+2,+0, -1,+1,+0,+1],
                           [+1,-2,+1,+0, -1,+0,+1,+1],
                           [+2,-1,-1,+0, +0,-1,+1,+1],
                           [+1,+1,-2,+0, +1,-1,+0,+1]]),
                         np.array([
                           [-2,+1,+1,+3, +1,+0,-1,+1],
                           [-1,-1,+2,+3, +1,+0,-1,+1],
                           [-1,-1,+2,+3, +0,+1,-1,+1],
                           [+1,-2,+1,+3, +0,+1,-1,+1],
                           [+1,-2,+1,+3, -1,+1,+0,+1],
                           [+2,-1,-1,+3, -1,+1,+0,+1],
                           [+2,-1,-1,+3, -1,+0,+1,+1],
                           [+1,+1,-2,+3, -1,+0,+1,+1],
                           [+1,+1,-2,+3, +0,-1,+1,+1],
                           [-1,+2,-1,+3, +0,-1,+1,+1],
                           [-1,+2,-1,+3, +1,-1,+0,+1],
                           [-2,+1,+1,+3, +1,-1,+0,+1]]),
                         np.array([
                           [-1,-1,+2,+3, +1,+1,-2,+2],
                           [+1,-2,+1,+3, -1,+2,-1,+2],
                           [+2,-1,-1,+3, -2,+1,+1,+2],
                           [+1,+1,-2,+3, -1,-1,+2,+2],
                           [-1,+2,-1,+3, +1,-2,+1,+2],
                           [-2,+1,+1,+3, +2,-1,-1,+2]])],
                'twin': [np.array([
                           [-1, 0, 1, 1,  1, 0,-1, 2],   # shear = (3-(c/a)^2)/(sqrt(3) c/a) <-10.1>{10.2}
                           [ 0,-1, 1, 1,  0, 1,-1, 2],
                           [ 1,-1, 0, 1, -1, 1, 0, 2],
                           [ 1, 0,-1, 1, -1, 0, 1, 2],
                           [ 0, 1,-1, 1,  0,-1, 1, 2],
                           [-1, 1, 0, 1,  1,-1, 0, 2]]),
                         np.array([
                           [-1,-1, 2, 6,  1, 1,-2, 1],   # shear = 1/(c/a) <11.6>{-1-1.1}
                           [ 1,-2, 1, 6, -1, 2,-1, 1],
                           [ 2,-1,-1, 6, -2, 1, 1, 1],
                           [ 1, 1,-2, 6, -1,-1, 2, 1],
                           [-1, 2,-1, 6,  1,-2, 1, 1],
                           [-2, 1, 1, 6,  2,-1,-1, 1]]),
                         np.array([
                           [ 1, 0,-1,-2,  1, 0,-1, 1],   # shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)  <10.-2>{10.1}
                           [ 0, 1,-1,-2,  0, 1,-1, 1],
                           [-1, 1, 0,-2, -1, 1, 0, 1],
                           [-1, 0, 1,-2, -1, 0, 1, 1],
                           [ 0,-1, 1,-2,  0,-1, 1, 1],
                           [ 1,-1, 0,-2,  1,-1, 0, 1]]),
                         np.array([
                           [ 1, 1,-2,-3,  1, 1,-2, 2],   # shear = 2((c/a)^2-2)/(3 c/a)  <11.-3>{11.2}
                           [-1, 2,-1,-3, -1, 2,-1, 2],
                           [-2, 1, 1,-3, -2, 1, 1, 2],
                           [-1,-1, 2,-3, -1,-1, 2, 2],
                           [ 1,-2, 1,-3,  1,-2, 1, 2],
                           [ 2,-1,-1,-3,  2,-1,-1, 2]])]
            },
            'tI': {
                'slip': [np.array([
                           [+0,+0,+1, +1,+0,+0],
                           [+0,+0,+1, +0,+1,+0]]),
                         np.array([
                           [+0,+0,+1, +1,+1,+0],
                           [+0,+0,+1, -1,+1,+0]]),
                         np.array([
                           [+0,+1,+0, +1,+0,+0],
                           [+1,+0,+0, +0,+1,+0]]),
                         np.array([
                           [+1,-1,+1, +1,+1,+0],
                           [+1,-1,-1, +1,+1,+0],
                           [-1,-1,-1, -1,+1,+0],
                           [-1,-1,+1, -1,+1,+0]]),
                         np.array([
                           [+1,-1,+0, +1,+1,+0],
                           [+1,+1,+0, +1,-1,+0]]),
                         np.array([
                           [+0,+1,+1, +1,+0,+0],
                           [+0,-1,+1, +1,+0,+0],
                           [-1,+0,+1, +0,+1,+0],
                           [+1,+0,+1, +0,+1,+0]]),
                         np.array([
                           [+0,+1,+0, +0,+0,+1],
                           [+1,+0,+0, +0,+0,+1]]),
                         np.array([
                           [+1,+1,+0, +0,+0,+1],
                           [-1,+1,+0, +0,+0,+1]]),
                         np.array([
                           [+0,+1,-1, +0,+1,+1],
                           [+0,-1,-1, +0,-1,+1],
                           [-1,+0,-1, -1,+0,+1],
                           [+1,+0,-1, +1,+0,+1]]),
                         np.array([
                           [+1,-1,+1, +0,+1,+1],
                           [+1,+1,-1, +0,+1,+1],
                           [+1,+1,+1, +0,+1,-1],
                           [-1,+1,+1, +0,+1,-1],
                           [+1,-1,-1, +1,+0,+1],
                           [-1,-1,+1, +1,+0,+1],
                           [+1,+1,+1, +1,+0,-1],
                           [+1,-1,+1, +1,+0,-1]]),
                         np.array([
                           [+1,+0,+0, +0,+1,+1],
                           [+1,+0,+0, +0,+1,-1],
                           [+0,+1,+0, +1,+0,+1],
                           [+0,+1,+0, +1,+0,-1]]),
                         np.array([
                           [+0,+1,-1, +2,+1,+1],
                           [+0,-1,-1, +2,-1,+1],
                           [+1,+0,-1, +1,+2,+1],
                           [-1,+0,-1, -1,+2,+1],
                           [+0,+1,-1, -2,+1,+1],
                           [+0,-1,-1, -2,-1,+1],
                           [-1,+0,-1, -1,-2,+1],
                           [+1,+0,-1, +1,-2,+1]]),
                         np.array([
                           [-1,+1,+1, +2,+1,+1],
                           [-1,-1,+1, +2,-1,+1],
                           [+1,-1,+1, +1,+2,+1],
                           [-1,-1,+1, -1,+2,+1],
                           [+1,+1,+1, -2,+1,+1],
                           [+1,-1,+1, -2,-1,+1],
                           [-1,+1,+1, -1,-2,+1],
                           [+1,+1,+1, +1,-2,+1]])]
                }
        }
        master = _kinematics[self.lattice][mode]
        return {'direction':[util.Bravais_to_Miller(uvtw=m[:,0:4]) if self.lattice == 'hP'
                                                    else m[:,0:3] for m in master],
                'plane':    [util.Bravais_to_Miller(hkil=m[:,4:8]) if self.lattice == 'hP'
                                                    else m[:,3:6] for m in master]}


    def relation_operations(self,
                            model: str) -> Tuple[BravaisLattice, Rotation]:
        """
        Crystallographic orientation relationships for phase transformations.

        Parameters
        ----------
        model : str
            Name of orientation relationship.

        Returns
        -------
        operations : (string, damask.Rotation)
            Resulting lattice and rotations characterizing the orientation relationship.

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
        my_relationships = {k:v for k,v in orientation_relationships.items() if self.lattice in v}
        if model not in my_relationships:
            raise KeyError(f'unknown orientation relationship "{model}"')
        r = my_relationships[model]

        sl = self.lattice
        ol = (set(r)-{sl}).pop()
        m = r[sl]
        o = r[ol]

        p_,_p = np.zeros(m.shape[:-1]+(3,)),np.zeros(o.shape[:-1]+(3,))
        p_[...,0,:] = m[...,0,:] if m.shape[-1] == 3 else util.Bravais_to_Miller(uvtw=m[...,0,0:4])
        p_[...,1,:] = m[...,1,:] if m.shape[-1] == 3 else util.Bravais_to_Miller(hkil=m[...,1,0:4])
        _p[...,0,:] = o[...,0,:] if o.shape[-1] == 3 else util.Bravais_to_Miller(uvtw=o[...,0,0:4])
        _p[...,1,:] = o[...,1,:] if o.shape[-1] == 3 else util.Bravais_to_Miller(hkil=o[...,1,0:4])

        return (ol,Rotation.from_parallel(p_,_p))
