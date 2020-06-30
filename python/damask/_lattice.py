import numpy as np

from . import Rotation


class Symmetry:
    """
    Symmetry-related operations for crystal systems.

    References
    ----------
    https://en.wikipedia.org/wiki/Crystal_system

    """

    crystal_systems = [None,'orthorhombic','tetragonal','hexagonal','cubic']

    def __init__(self, system = None):
        """
        Symmetry Definition.

        Parameters
        ----------
        system : {None,'orthorhombic','tetragonal','hexagonal','cubic'}, optional
            Name of the crystal system. Defaults to 'None'.

        """
        if system is not None and system.lower() not in self.crystal_systems:
            raise KeyError(f'Crystal system "{system}" is unknown')

        self.system = system.lower() if isinstance(system,str) else system
        self.lattice = self.system # for compatibility


    def __copy__(self):
        """Copy."""
        return self.__class__(self.system)

    copy = __copy__


    def __repr__(self):
        """Readable string."""
        return '{}'.format(self.system)


    def __eq__(self, other):
        """
        Equal to other.

        Parameters
        ----------
        other : Symmetry
            Symmetry to check for equality.

        """
        return self.system == other.system

    def __neq__(self, other):
        """
        Not Equal to other.

        Parameters
        ----------
        other : Symmetry
            Symmetry to check for inequality.

        """
        return not self.__eq__(other)

    def __cmp__(self,other):
        """
        Linear ordering.

        Parameters
        ----------
        other : Symmetry
            Symmetry to check for for order.

        """
        myOrder    = self.crystal_systems.index(self.system)
        otherOrder = self.crystal_systems.index(other.system)
        return (myOrder > otherOrder) - (myOrder < otherOrder)

    def symmetryOperations(self,members=[]):
        """List (or single element) of symmetry operations as rotations."""
        if self.system == 'cubic':
            symQuats =  [
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
        elif self.system == 'hexagonal':
            symQuats =  [
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
        elif self.system == 'tetragonal':
            symQuats =  [
                          [ 1.0,            0.0,            0.0,            0.0            ],
                          [ 0.0,            1.0,            0.0,            0.0            ],
                          [ 0.0,            0.0,            1.0,            0.0            ],
                          [ 0.0,            0.0,            0.0,            1.0            ],
                          [ 0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                          [ 0.0,           -0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                          [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                          [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                        ]
        elif self.system == 'orthorhombic':
            symQuats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                          [ 0.0,1.0,0.0,0.0 ],
                          [ 0.0,0.0,1.0,0.0 ],
                          [ 0.0,0.0,0.0,1.0 ],
                        ]
        else:
            symQuats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                        ]

        symOps = list(map(Rotation,
                      np.array(symQuats)[np.atleast_1d(members) if members != [] else range(len(symQuats))]))
        try:
            iter(members)                                                                           # asking for (even empty) list of members?
        except TypeError:
            return symOps[0]                                                                        # no, return rotation object
        else:
            return symOps                                                                           # yes, return list of rotations

    @property
    def symmetry_operations(self):
        """Symmetry operations as Rotations."""
        if self.system == 'cubic':
            symQuats =  [
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
        elif self.system == 'hexagonal':
            symQuats =  [
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
        elif self.system == 'tetragonal':
            symQuats =  [
                          [ 1.0,            0.0,            0.0,            0.0            ],
                          [ 0.0,            1.0,            0.0,            0.0            ],
                          [ 0.0,            0.0,            1.0,            0.0            ],
                          [ 0.0,            0.0,            0.0,            1.0            ],
                          [ 0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                          [ 0.0,           -0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                          [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                          [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                        ]
        elif self.system == 'orthorhombic':
            symQuats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                          [ 0.0,1.0,0.0,0.0 ],
                          [ 0.0,0.0,1.0,0.0 ],
                          [ 0.0,0.0,0.0,1.0 ],
                        ]
        else:
            symQuats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                        ]
        return np.array(symQuats)


    def inFZ(self,rodrigues):
        """
        Check whether given Rodrigues-Frank vector falls into fundamental zone of own symmetry.

        Fundamental zone in Rodrigues space is point symmetric around origin.
        """
        if (len(rodrigues) != 3):
            raise ValueError('Input is not a Rodrigues-Frank vector.\n')

        if np.any(rodrigues == np.inf): return False # ToDo: MD: not sure if needed

        Rabs = abs(rodrigues)

        if self.system == 'cubic':
            return     np.sqrt(2.0)-1.0 >= Rabs[0] \
                   and np.sqrt(2.0)-1.0 >= Rabs[1] \
                   and np.sqrt(2.0)-1.0 >= Rabs[2] \
                   and 1.0 >= Rabs[0] + Rabs[1] + Rabs[2]
        elif self.system == 'hexagonal':
            return     1.0 >= Rabs[0] and 1.0 >= Rabs[1] and 1.0 >= Rabs[2] \
                   and 2.0 >= np.sqrt(3)*Rabs[0] + Rabs[1] \
                   and 2.0 >= np.sqrt(3)*Rabs[1] + Rabs[0] \
                   and 2.0 >= np.sqrt(3) + Rabs[2]
        elif self.system == 'tetragonal':
            return     1.0 >= Rabs[0] and 1.0 >= Rabs[1] \
                   and np.sqrt(2.0) >= Rabs[0] + Rabs[1] \
                   and np.sqrt(2.0) >= Rabs[2] + 1.0
        elif self.system == 'orthorhombic':
            return     1.0 >= Rabs[0] and 1.0 >= Rabs[1] and 1.0 >= Rabs[2]
        else:
            return True


    def inDisorientationSST(self,rodrigues):
        """
        Check whether given Rodrigues-Frank vector (of misorientation) falls into standard stereographic triangle of own symmetry.

        References
        ----------
        A. Heinz and P. Neumann, Acta Crystallographica Section A 47:780-789, 1991
        https://doi.org/10.1107/S0108767391006864

        """
        if (len(rodrigues) != 3):
            raise ValueError('Input is not a Rodrigues-Frank vector.\n')
        R = rodrigues

        epsilon = 0.0
        if self.system == 'cubic':
            return R[0] >= R[1]+epsilon              and R[1] >= R[2]+epsilon and R[2] >= epsilon
        elif self.system == 'hexagonal':
            return R[0] >= np.sqrt(3)*(R[1]-epsilon) and R[1] >= epsilon      and R[2] >= epsilon
        elif self.system == 'tetragonal':
            return R[0] >= R[1]-epsilon              and R[1] >= epsilon      and R[2] >= epsilon
        elif self.system == 'orthorhombic':
            return R[0] >= epsilon                   and R[1] >= epsilon      and R[2] >= epsilon
        else:
            return True


    def inSST(self,
              vector,
              proper = False,
              color = False):
        """
        Check whether given vector falls into standard stereographic triangle of own symmetry.

        proper considers only vectors with z >= 0, hence uses two neighboring SSTs.
        Return inverse pole figure color if requested.
        Bases are computed from

        >>> basis = {'cubic' :       np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,1.]/np.sqrt(2.),                # direction of green
        ...                                                  [1.,1.,1.]/np.sqrt(3.)]).T),           # direction of blue
        ...          'hexagonal' :   np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,0.],                            # direction of green
        ...                                                  [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T),  # direction of blue
        ...          'tetragonal' :  np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,0.],                            # direction of green
        ...                                                  [1.,1.,0.]/np.sqrt(2.)]).T),           # direction of blue
        ...          'orthorhombic': np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,0.],                            # direction of green
        ...                                                  [0.,1.,0.]]).T),                       # direction of blue
        ...         }

        """
        if self.system == 'cubic':
            basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                           [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                           [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                       'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                           [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                           [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                    }
        elif self.system == 'hexagonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -np.sqrt(3.)   ,  0. ],
                                           [ 0.            ,  2.            ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  np.sqrt(3.)   ,  0. ],
                                           [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                    }
        elif self.system == 'tetragonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -1.            ,  0. ],
                                           [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  1.            ,  0. ],
                                           [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                    }
        elif self.system == 'orthorhombic':
            basis = {'improper':np.array([ [ 0., 0., 1.],
                                           [ 1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                       'proper':np.array([ [ 0., 0., 1.],
                                           [-1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                    }
        else:                                                                                       # direct exit for unspecified symmetry
            if color:
                return (True,np.zeros(3,'d'))
            else:
                return True

        v = np.array(vector,dtype=float)
        if proper:                                                                                  # check both improper ...
            theComponents = np.around(np.dot(basis['improper'],v),12)
            inSST = np.all(theComponents >= 0.0)
            if not inSST:                                                                           # ... and proper SST
                theComponents = np.around(np.dot(basis['proper'],v),12)
                inSST = np.all(theComponents >= 0.0)
        else:
            v[2] = abs(v[2])                                                                        # z component projects identical
            theComponents = np.around(np.dot(basis['improper'],v),12)                               # for positive and negative values
            inSST = np.all(theComponents >= 0.0)

        if color:                                                                                   # have to return color array
            if inSST:
                rgb = np.power(theComponents/np.linalg.norm(theComponents),0.5)                     # smoothen color ramps
                rgb = np.minimum(np.ones(3,dtype=float),rgb)                                        # limit to maximum intensity
                rgb /= max(rgb)                                                                     # normalize to (HS)V = 1
            else:
                rgb = np.zeros(3,dtype=float)
            return (inSST,rgb)
        else:
            return inSST


    def in_SST(self,
               vector,
               proper = False,
               color = False):
        """
        Check whether given vector falls into standard stereographic triangle of own symmetry.

        proper considers only vectors with z >= 0, hence uses two neighboring SSTs.
        Return inverse pole figure color if requested.
        Bases are computed from

        >>> basis = {'cubic' :       np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,1.]/np.sqrt(2.),                # direction of green
        ...                                                  [1.,1.,1.]/np.sqrt(3.)]).T),           # direction of blue
        ...          'hexagonal' :   np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,0.],                            # direction of green
        ...                                                  [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T),  # direction of blue
        ...          'tetragonal' :  np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,0.],                            # direction of green
        ...                                                  [1.,1.,0.]/np.sqrt(2.)]).T),           # direction of blue
        ...          'orthorhombic': np.linalg.inv(np.array([[0.,0.,1.],                            # direction of red
        ...                                                  [1.,0.,0.],                            # direction of green
        ...                                                  [0.,1.,0.]]).T),                       # direction of blue
        ...         }

        """
        if self.system == 'cubic':
            basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                           [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                           [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                       'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                           [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                           [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                    }
        elif self.system == 'hexagonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -np.sqrt(3.)   ,  0. ],
                                           [ 0.            ,  2.            ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  np.sqrt(3.)   ,  0. ],
                                           [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                    }
        elif self.system == 'tetragonal':
            basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                           [ 1.            , -1.            ,  0. ],
                                           [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                     'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                           [-1.            ,  1.            ,  0. ],
                                           [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                    }
        elif self.system == 'orthorhombic':
            basis = {'improper':np.array([ [ 0., 0., 1.],
                                           [ 1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                       'proper':np.array([ [ 0., 0., 1.],
                                           [-1., 0., 0.],
                                           [ 0., 1., 0.] ]),
                    }
        else:                                                                                       # direct exit for unspecified symmetry
            if color:
                return (np.ones_like(vector[...,0],bool),np.zeros_like(vector))
            else:
                return  np.ones_like(vector[...,0],bool)

        b_p = np.broadcast_to(basis['proper'], vector.shape+(3,))
        if proper:
            b_i = np.broadcast_to(basis['improper'],vector.shape+(3,))
            improper = np.all(np.around(np.einsum('...ji,...i',b_i,vector),12)>=0.0,axis=-1,keepdims=True)
            theComponents = np.where(np.broadcast_to(improper,vector.shape),
                                     np.around(np.einsum('...ji,...i',b_i,vector),12),
                                     np.around(np.einsum('...ji,...i',b_p,vector),12))
        else:
            vector_ = np.block([vector[...,0:2],np.abs(vector[...,2:3])])                           # z component projects identical
            theComponents = np.around(np.einsum('...ji,...i',b_p,vector_),12)

        in_SST = np.all(theComponents >= 0.0,axis=-1)

        if color:                                                                                   # have to return color array
            with np.errstate(invalid='ignore',divide='ignore'):
                rgb = (theComponents/np.linalg.norm(theComponents,axis=-1,keepdims=True))**0.5      # smoothen color ramps
                rgb = np.minimum(1.,rgb)                                                            # limit to maximum intensity
                rgb /= np.max(rgb,axis=-1,keepdims=True)                                            # normalize to (HS)V = 1
            rgb[np.broadcast_to(~in_SST.reshape(vector[...,0].shape+(1,)),vector.shape)] = 0.0
            return (in_SST,rgb)
        else:
            return in_SST


# ******************************************************************************************
class Lattice: # ToDo: Make a subclass of Symmetry!
    """
    Lattice system.

    Currently, this contains only a mapping from Bravais lattice to symmetry
    and orientation relationships. It could include twin and slip systems.

    References
    ----------
    https://en.wikipedia.org/wiki/Bravais_lattice

    """

    lattices = {
                'triclinic':{'symmetry':None},
                'bct':{'symmetry':'tetragonal'},
                'hex':{'symmetry':'hexagonal'},
                'fcc':{'symmetry':'cubic','c/a':1.0},
                'bcc':{'symmetry':'cubic','c/a':1.0},
               }


    def __init__(self, lattice):
        """
        New lattice of given type.

        Parameters
        ----------
        lattice : str
            Bravais lattice.

        """
        self.lattice  = lattice
        self.symmetry = Symmetry(self.lattices[lattice]['symmetry'])


    def __repr__(self):
        """Report basic lattice information."""
        return 'Bravais lattice {} ({} symmetry)'.format(self.lattice,self.symmetry)


    # Kurdjomov--Sachs orientation relationship for fcc <-> bcc transformation
    # from S. Morito et al., Journal of Alloys and Compounds 577:s587-s592, 2013
    # also see K. Kitahara et al., Acta Materialia 54:1279-1288, 2006
    KS = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[  1,  1, -1],[  0,  1,  1]],
        [[  1,  1, -1],[  0,  1,  1]],
        [[  1,  1, -1],[  0,  1,  1]],
        [[  1,  1, -1],[  0,  1,  1]],
        [[  1,  1, -1],[  0,  1,  1]],
        [[  1,  1, -1],[  0,  1,  1]]],dtype='float'),
        'directions': np.array([
        [[ -1,  0,  1],[ -1, -1,  1]],
        [[ -1,  0,  1],[ -1,  1, -1]],
        [[  0,  1, -1],[ -1, -1,  1]],
        [[  0,  1, -1],[ -1,  1, -1]],
        [[  1, -1,  0],[ -1, -1,  1]],
        [[  1, -1,  0],[ -1,  1, -1]],
        [[  1,  0, -1],[ -1, -1,  1]],
        [[  1,  0, -1],[ -1,  1, -1]],
        [[ -1, -1,  0],[ -1, -1,  1]],
        [[ -1, -1,  0],[ -1,  1, -1]],
        [[  0,  1,  1],[ -1, -1,  1]],
        [[  0,  1,  1],[ -1,  1, -1]],
        [[  0, -1,  1],[ -1, -1,  1]],
        [[  0, -1,  1],[ -1,  1, -1]],
        [[ -1,  0, -1],[ -1, -1,  1]],
        [[ -1,  0, -1],[ -1,  1, -1]],
        [[  1,  1,  0],[ -1, -1,  1]],
        [[  1,  1,  0],[ -1,  1, -1]],
        [[ -1,  1,  0],[ -1, -1,  1]],
        [[ -1,  1,  0],[ -1,  1, -1]],
        [[  0, -1, -1],[ -1, -1,  1]],
        [[  0, -1, -1],[ -1,  1, -1]],
        [[  1,  0,  1],[ -1, -1,  1]],
        [[  1,  0,  1],[ -1,  1, -1]]],dtype='float')}

    # Greninger--Troiano orientation relationship for fcc <-> bcc transformation
    # from Y. He et al., Journal of Applied Crystallography 39:72-81, 2006
    GT = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  1,  1,  1],[  1,  0,  1]],
        [[  1,  1,  1],[  1,  1,  0]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[ -1, -1,  1],[ -1,  0,  1]],
        [[ -1, -1,  1],[ -1, -1,  0]],
        [[ -1, -1,  1],[  0, -1,  1]],
        [[ -1,  1,  1],[ -1,  0,  1]],
        [[ -1,  1,  1],[ -1,  1,  0]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  1,  0,  1]],
        [[  1, -1,  1],[  1, -1,  0]],
        [[  1, -1,  1],[  0, -1,  1]],
        [[  1,  1,  1],[  1,  1,  0]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  1,  0,  1]],
        [[ -1, -1,  1],[ -1, -1,  0]],
        [[ -1, -1,  1],[  0, -1,  1]],
        [[ -1, -1,  1],[ -1,  0,  1]],
        [[ -1,  1,  1],[ -1,  1,  0]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[ -1,  0,  1]],
        [[  1, -1,  1],[  1, -1,  0]],
        [[  1, -1,  1],[  0, -1,  1]],
        [[  1, -1,  1],[  1,  0,  1]]],dtype='float'),
        'directions': np.array([
        [[ -5,-12, 17],[-17, -7, 17]],
        [[ 17, -5,-12],[ 17,-17, -7]],
        [[-12, 17, -5],[ -7, 17,-17]],
        [[  5, 12, 17],[ 17,  7, 17]],
        [[-17,  5,-12],[-17, 17, -7]],
        [[ 12,-17, -5],[  7,-17,-17]],
        [[ -5, 12,-17],[-17,  7,-17]],
        [[ 17,  5, 12],[ 17, 17,  7]],
        [[-12,-17,  5],[ -7,-17, 17]],
        [[  5,-12,-17],[ 17, -7,-17]],
        [[-17, -5, 12],[-17,-17,  7]],
        [[ 12, 17,  5],[  7, 17, 17]],
        [[ -5, 17,-12],[-17, 17, -7]],
        [[-12, -5, 17],[ -7,-17, 17]],
        [[ 17,-12, -5],[ 17, -7,-17]],
        [[  5,-17,-12],[ 17,-17, -7]],
        [[ 12,  5, 17],[  7, 17, 17]],
        [[-17, 12, -5],[-17,  7,-17]],
        [[ -5,-17, 12],[-17,-17,  7]],
        [[-12,  5,-17],[ -7, 17,-17]],
        [[ 17, 12,  5],[ 17,  7, 17]],
        [[  5, 17, 12],[ 17, 17,  7]],
        [[ 12, -5,-17],[  7,-17,-17]],
        [[-17,-12,  5],[-17,-7, 17]]],dtype='float')}

    # Greninger--Troiano' orientation relationship for fcc <-> bcc transformation
    # from Y. He et al., Journal of Applied Crystallography 39:72-81, 2006
    GTprime = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  7, 17, 17],[ 12,  5, 17]],
        [[ 17,  7, 17],[ 17, 12,  5]],
        [[ 17, 17,  7],[  5, 17, 12]],
        [[ -7,-17, 17],[-12, -5, 17]],
        [[-17, -7, 17],[-17,-12,  5]],
        [[-17,-17,  7],[ -5,-17, 12]],
        [[  7,-17,-17],[ 12, -5,-17]],
        [[ 17, -7,-17],[ 17,-12, -5]],
        [[ 17,-17, -7],[  5,-17,-12]],
        [[ -7, 17,-17],[-12,  5,-17]],
        [[-17,  7,-17],[-17, 12, -5]],
        [[-17, 17, -7],[ -5, 17,-12]],
        [[  7, 17, 17],[ 12, 17,  5]],
        [[ 17,  7, 17],[  5, 12, 17]],
        [[ 17, 17,  7],[ 17,  5, 12]],
        [[ -7,-17, 17],[-12,-17,  5]],
        [[-17, -7, 17],[ -5,-12, 17]],
        [[-17,-17,  7],[-17, -5, 12]],
        [[  7,-17,-17],[ 12,-17, -5]],
        [[ 17, -7,-17],[ 5, -12,-17]],
        [[ 17,-17, -7],[ 17, -5,-12]],
        [[ -7, 17,-17],[-12, 17, -5]],
        [[-17,  7,-17],[ -5, 12,-17]],
        [[-17, 17, -7],[-17,  5,-12]]],dtype='float'),
        'directions': np.array([
        [[  0,  1, -1],[  1,  1, -1]],
        [[ -1,  0,  1],[ -1,  1,  1]],
        [[  1, -1,  0],[  1, -1,  1]],
        [[  0, -1, -1],[ -1, -1, -1]],
        [[  1,  0,  1],[  1, -1,  1]],
        [[  1, -1,  0],[  1, -1, -1]],
        [[  0,  1, -1],[ -1,  1, -1]],
        [[  1,  0,  1],[  1,  1,  1]],
        [[ -1, -1,  0],[ -1, -1,  1]],
        [[  0, -1, -1],[  1, -1, -1]],
        [[ -1,  0,  1],[ -1, -1,  1]],
        [[ -1, -1,  0],[ -1, -1, -1]],
        [[  0, -1,  1],[  1, -1,  1]],
        [[  1,  0, -1],[  1,  1, -1]],
        [[ -1,  1,  0],[ -1,  1,  1]],
        [[  0,  1,  1],[ -1,  1,  1]],
        [[ -1,  0, -1],[ -1, -1, -1]],
        [[ -1,  1,  0],[ -1,  1, -1]],
        [[  0, -1,  1],[ -1, -1,  1]],
        [[ -1,  0, -1],[ -1,  1, -1]],
        [[  1,  1,  0],[  1,  1,  1]],
        [[  0,  1,  1],[  1,  1,  1]],
        [[  1,  0, -1],[  1, -1, -1]],
        [[  1,  1,  0],[  1,  1, -1]]],dtype='float')}

    # Nishiyama--Wassermann orientation relationship for fcc <-> bcc transformation
    # from H. Kitahara et al., Materials Characterization 54:378-386, 2005
    NW = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[  1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[ -1,  1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[  1, -1,  1],[  0,  1,  1]],
        [[ -1, -1,  1],[  0,  1,  1]],
        [[ -1, -1,  1],[  0,  1,  1]],
        [[ -1, -1,  1],[  0,  1,  1]]],dtype='float'),
        'directions': np.array([
        [[  2, -1, -1],[  0, -1,  1]],
        [[ -1,  2, -1],[  0, -1,  1]],
        [[ -1, -1,  2],[  0, -1,  1]],
        [[ -2, -1, -1],[  0, -1,  1]],
        [[  1,  2, -1],[  0, -1,  1]],
        [[  1, -1,  2],[  0, -1,  1]],
        [[  2,  1, -1],[  0, -1,  1]],
        [[ -1, -2, -1],[  0, -1,  1]],
        [[ -1,  1,  2],[  0, -1,  1]],
        [[  2, -1,  1],[  0, -1,  1]], #It is wrong in the paper, but matrix is correct
        [[ -1,  2,  1],[  0, -1,  1]],
        [[ -1, -1, -2],[  0, -1,  1]]],dtype='float')}

    # Pitsch orientation relationship for fcc <-> bcc transformation
    # from Y. He et al., Acta Materialia 53:1179-1190, 2005
    Pitsch = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  0,  1,  0],[ -1,  0,  1]],
        [[  0,  0,  1],[  1, -1,  0]],
        [[  1,  0,  0],[  0,  1, -1]],
        [[  1,  0,  0],[  0, -1, -1]],
        [[  0,  1,  0],[ -1,  0, -1]],
        [[  0,  0,  1],[ -1, -1,  0]],
        [[  0,  1,  0],[ -1,  0, -1]],
        [[  0,  0,  1],[ -1, -1,  0]],
        [[  1,  0,  0],[  0, -1, -1]],
        [[  1,  0,  0],[  0, -1,  1]],
        [[  0,  1,  0],[  1,  0, -1]],
        [[  0,  0,  1],[ -1,  1,  0]]],dtype='float'),
        'directions': np.array([
        [[  1,  0,  1],[  1, -1,  1]],
        [[  1,  1,  0],[  1,  1, -1]],
        [[  0,  1,  1],[ -1,  1,  1]],
        [[  0,  1, -1],[ -1,  1, -1]],
        [[ -1,  0,  1],[ -1, -1,  1]],
        [[  1, -1,  0],[  1, -1, -1]],
        [[  1,  0, -1],[  1, -1, -1]],
        [[ -1,  1,  0],[ -1,  1, -1]],
        [[  0, -1,  1],[ -1, -1,  1]],
        [[  0,  1,  1],[ -1,  1,  1]],
        [[  1,  0,  1],[  1, -1,  1]],
        [[  1,  1,  0],[  1,  1, -1]]],dtype='float')}

    # Bain orientation relationship for fcc <-> bcc transformation
    # from Y. He et al., Journal of Applied Crystallography 39:72-81, 2006
    Bain = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  1,  0,  0],[  1,  0,  0]],
        [[  0,  1,  0],[  0,  1,  0]],
        [[  0,  0,  1],[  0,  0,  1]]],dtype='float'),
        'directions': np.array([
        [[  0,  1,  0],[  0,  1,  1]],
        [[  0,  0,  1],[  1,  0,  1]],
        [[  1,  0,  0],[  1,  1,  0]]],dtype='float')}

    def relationOperations(self,model):
        """
        Crystallographic orientation relationships for phase transformations.

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
        models={'KS':self.KS, 'GT':self.GT, 'GT_prime':self.GTprime,
                'NW':self.NW, 'Pitsch': self.Pitsch, 'Bain':self.Bain}
        try:
            relationship = models[model]
        except KeyError :
            raise KeyError('Orientation relationship "{}" is unknown'.format(model))

        if self.lattice not in relationship['mapping']:
            raise ValueError('Relationship "{}" not supported for lattice "{}"'.format(model,self.lattice))

        r = {'lattice':Lattice((set(relationship['mapping'])-{self.lattice}).pop()),                # target lattice
             'rotations':[] }

        myPlane_id    = relationship['mapping'][self.lattice]
        otherPlane_id = (myPlane_id+1)%2
        myDir_id      = myPlane_id +2
        otherDir_id   = otherPlane_id +2

        for miller in np.hstack((relationship['planes'],relationship['directions'])):
            myPlane     = miller[myPlane_id]/    np.linalg.norm(miller[myPlane_id])
            myDir       = miller[myDir_id]/      np.linalg.norm(miller[myDir_id])
            myMatrix    = np.array([myDir,np.cross(myPlane,myDir),myPlane])

            otherPlane  = miller[otherPlane_id]/ np.linalg.norm(miller[otherPlane_id])
            otherDir    = miller[otherDir_id]/   np.linalg.norm(miller[otherDir_id])
            otherMatrix = np.array([otherDir,np.cross(otherPlane,otherDir),otherPlane])

            r['rotations'].append(Rotation.from_matrix(np.dot(otherMatrix.T,myMatrix)))

        return r
