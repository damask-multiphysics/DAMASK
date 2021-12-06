from typing import Union, Dict, List, Tuple

import numpy as np

from . import util
from . import Rotation

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


class Crystal():
    """Crystal lattice."""

    def __init__(self,*,
                 family = None,
                 lattice = None,
                 a = None,b = None,c = None,
                 alpha = None,beta = None,gamma = None,
                 degrees = False):
        """
        Representation of crystal in terms of crystal family or Bravais lattice.

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
        if family not in [None] + list(lattice_symmetries.values()):
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
        else:
            self.a = self.b = self.c = None
            self.alpha = self.beta = self.gamma = None


    def __repr__(self):
        """Represent."""
        return '\n'.join([f'Crystal family {self.family}']
                       + ([] if self.lattice is None else [f'Bravais lattice {self.lattice}']+
                                                           list(map(lambda x:f'{x[0]}: {x[1]:.5g}',
                                                                    zip(['a','b','c','α','β','γ',],
                                                                         self.parameters))))
                        )

    def __eq__(self,other):
        """
        Equal to other.

        Parameters
        ----------
        other : Crystal
            Crystal to check for equality.

        """
        return self.lattice == other.lattice and \
               self.parameters == other.parameters and \
               self.family == other.family

    @property
    def parameters(self):
        """Return lattice parameters a, b, c, alpha, beta, gamma."""
        return (self.a,self.b,self.c,self.alpha,self.beta,self.gamma)


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
        _basis  = {
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
        if None in self.parameters:
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

    def to_lattice(self, *, direction: np.ndarray = None, plane: np.ndarray = None) -> np.ndarray:
        """
        Calculate lattice vector corresponding to crystal frame direction or plane normal.

        Parameters
        ----------
        direction|plane : numpy.ndarray of shape (...,3)
            Vector along direction or plane normal.

        Returns
        -------
        Miller : numpy.ndarray of shape (...,3)
            Lattice vector of direction or plane.
            Use util.scale_to_coprime to convert to (integer) Miller indices.

        """
        if (direction is not None) ^ (plane is None):
            raise KeyError('specify either "direction" or "plane"')
        axis,basis  = (np.array(direction),self.basis_reciprocal.T) \
                      if plane is None else \
                      (np.array(plane),self.basis_real.T)
        return np.einsum('il,...l',basis,axis)


    def to_frame(self, *, uvw: np.ndarray = None, hkl: np.ndarray = None) -> np.ndarray:
        """
        Calculate crystal frame vector along lattice direction [uvw] or plane normal (hkl).

        Parameters
        ----------
        uvw|hkl : numpy.ndarray of shape (...,3)
            Miller indices of crystallographic direction or plane normal.

        Returns
        -------
        vector : numpy.ndarray of shape (...,3)
            Crystal frame vector along [uvw] direction or (hkl) plane normal.

        """
        if (uvw is not None) ^ (hkl is None):
            raise KeyError('specify either "uvw" or "hkl"')
        axis,basis  = (np.array(uvw),self.basis_real) \
                      if hkl is None else \
                      (np.array(hkl),self.basis_reciprocal)
        return np.einsum('il,...l',basis,axis)


    def kinematics(self, mode: str) -> Dict[str, List[np.ndarray]]:
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
        _kinematics = {
            'cF': {
                'slip' :[np.array([
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
                'twin' :[np.array([
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
                'slip' :[np.array([
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
                'twin' :[np.array([
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
                'slip' :[np.array([
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
                'twin' :[np.array([
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
        }
        master = _kinematics[self.lattice][mode]
        if self.lattice == 'hP':
            return {'direction':[util.Bravais_to_Miller(uvtw=m[:,0:4]) for m in master],
                    'plane':    [util.Bravais_to_Miller(hkil=m[:,4:8]) for m in master]}
        else:
            return {'direction':[m[:,0:3] for m in master],
                    'plane':    [m[:,3:6] for m in master]}


    def relation_operations(self, model: str) -> Tuple[str, Rotation]:
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
        _orientation_relationships = {
          'KS': {
            'cF' : np.array([
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
            'cI' : np.array([
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
            'cF' : np.array([
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
            'cI' : np.array([
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
        orientation_relationships = {k:v for k,v in _orientation_relationships.items() if self.lattice in v}
        if model not in orientation_relationships:
            raise KeyError(f'unknown orientation relationship "{model}"')
        r = orientation_relationships[model]

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
