import numpy as np

from . import Rotation

class LatticeFamily():

    def __init__(self,family):
        """
        Symmetry-related operations for crystal family.

        Parameters
        ----------
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}
            Name of the crystal family.

        """
        if family not in self._immutable.keys():
            raise KeyError(f'invalid crystal family "{family}"')
        self.family = family


    def __eq__(self,other):
        """
        Equal to other.

        Parameters
        ----------
        other : LatticeFamily
            Lattice family to check for equality.

        """
        return self.family == other.family


    @property
    def symmetry_operations(self):
        """Symmetry operations as Rotations."""
        return Rotation.from_quaternion(self._symmetry_operations[self.family],accept_homomorph=True)


    @property
    def immutable(self):
        """Return immutable lattice parameters."""
        return self._immutable[self.family]


    @property
    def basis(self):
        """
        Corners of the standard triangle.

        Not yet defined for monoclinic.


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
        return self._basis.get(self.family,None)


    _symmetry_operations = {
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
