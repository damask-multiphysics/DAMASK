import numpy as np


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


    def __copy__(self):
        """Copy."""
        return self.__class__(self.system)

    copy = __copy__


    def __repr__(self):
        """Readable string."""
        return f'{self.system}'


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


    @property
    def symmetry_operations(self):
        """Symmetry operations as quaternions."""
        if self.system == 'cubic':
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
        elif self.system == 'hexagonal':
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
        elif self.system == 'tetragonal':
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
        elif self.system == 'orthorhombic':
            sym_quats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                          [ 0.0,1.0,0.0,0.0 ],
                          [ 0.0,0.0,1.0,0.0 ],
                          [ 0.0,0.0,0.0,1.0 ],
                        ]
        else:
            sym_quats =  [
                          [ 1.0,0.0,0.0,0.0 ],
                        ]
        return np.array(sym_quats)


    def in_FZ(self,rho):
        """
        Check whether given Rodrigues-Frank vector falls into fundamental zone.

        Fundamental zone in Rodrigues space is point symmetric around origin.
        """
        if(rho.shape[-1] != 3):
            raise ValueError('Input is not a Rodrigues-Frank vector field.')

        rho_abs = np.abs(rho)

        with np.errstate(invalid='ignore'):
            # using '*'/prod for 'and'
            if self.system == 'cubic':
                return np.where(np.prod(np.sqrt(2)-1. >= rho_abs,axis=-1) * \
                                (1. >= np.sum(rho_abs,axis=-1)),True,False)
            elif self.system == 'hexagonal':
                return np.where(np.prod(1.             >= rho_abs,axis=-1) * \
                                (2. >= np.sqrt(3)*rho_abs[...,0] + rho_abs[...,1]) * \
                                (2. >= np.sqrt(3)*rho_abs[...,1] + rho_abs[...,0]) * \
                                (2. >= np.sqrt(3) + rho_abs[...,2]),True,False)
            elif self.system == 'tetragonal':
                return np.where(np.prod(1.             >= rho_abs[...,:2],axis=-1) * \
                                (np.sqrt(2) >= rho_abs[...,0] + rho_abs[...,1]) * \
                                (np.sqrt(2) >= rho_abs[...,2] + 1.),True,False)
            elif self.system == 'orthorhombic':
                return np.where(np.prod(1.             >= rho_abs,axis=-1),True,False)
            else:
                return np.where(np.all(np.isfinite(rho_abs),axis=-1),True,False)


    def in_disorientation_SST(self,rho):
        """
        Check whether given Rodrigues-Frank vector (of misorientation) falls into standard stereographic triangle.

        References
        ----------
        A. Heinz and P. Neumann, Acta Crystallographica Section A 47:780-789, 1991
        https://doi.org/10.1107/S0108767391006864

        """
        if(rho.shape[-1] != 3):
            raise ValueError('Input is not a Rodrigues-Frank vector field.')

        with np.errstate(invalid='ignore'):
            # using '*' for 'and'
            if self.system == 'cubic':
                return np.where((rho[...,0] >= rho[...,1]) * \
                                (rho[...,1] >= rho[...,2]) * \
                                (rho[...,2] >= 0),True,False)
            elif self.system == 'hexagonal':
                return np.where((rho[...,0] >= rho[...,1]*np.sqrt(3)) * \
                                (rho[...,1] >= 0) * \
                                (rho[...,2] >= 0),True,False)
            elif self.system == 'tetragonal':
                return np.where((rho[...,0] >= rho[...,1]) * \
                                (rho[...,1] >= 0) * \
                                (rho[...,2] >= 0),True,False)
            elif self.system == 'orthorhombic':
                return np.where((rho[...,0] >= 0) * \
                                (rho[...,1] >= 0) * \
                                (rho[...,2] >= 0),True,False)
            else:
                return np.ones_like(rho[...,0],dtype=bool)


    #ToDo: IPF color in separate function
    def in_SST(self,vector,proper=False,color=False):
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
        if(vector.shape[-1] != 3):
            raise ValueError('Input is not a 3D vector field.')

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


        b_i = np.broadcast_to(basis['improper'],vector.shape+(3,))
        if proper:
            b_p = np.broadcast_to(basis['proper'], vector.shape+(3,))
            improper = np.all(np.around(np.einsum('...ji,...i',b_i,vector),12)>=0.0,axis=-1,keepdims=True)
            theComponents = np.where(np.broadcast_to(improper,vector.shape),
                                     np.around(np.einsum('...ji,...i',b_i,vector),12),
                                     np.around(np.einsum('...ji,...i',b_p,vector),12))
        else:
            vector_ = np.block([vector[...,0:2],np.abs(vector[...,2:3])])                           # z component projects identical
            theComponents = np.around(np.einsum('...ji,...i',b_i,vector_),12)

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
    Bravais lattice.

    This contains only a mapping from Bravais lattice to symmetry
    and orientation relationships. It could include twin and slip systems.

    References
    ----------
    https://en.wikipedia.org/wiki/Bravais_lattice

    """

    lattices = {
                'triclinic':{'system':None},
                'bct':      {'system':'tetragonal'},
                'hex':      {'system':'hexagonal'},
                'fcc':      {'system':'cubic','c/a':1.0},
                'bcc':      {'system':'cubic','c/a':1.0},
               }


    def __init__(self,lattice,c_over_a=None):
        """
        New lattice of given type.

        Parameters
        ----------
        lattice : str
            Bravais lattice.

        """
        self.lattice  = lattice
        self.symmetry = Symmetry(self.lattices[lattice]['system'])

        # transition to subclass
        self.system                 = self.symmetry.system
        self.in_SST                 = self.symmetry.in_SST
        self.in_FZ                  = self.symmetry.in_FZ
        self.in_disorientation_SST  = self.symmetry.in_disorientation_SST

    def __repr__(self):
        """Report basic lattice information."""
        return f'Bravais lattice {self.lattice} ({self.symmetry} crystal system)'


    # Kurdjomov--Sachs orientation relationship for fcc <-> bcc transformation
    # from S. Morito et al., Journal of Alloys and Compounds 577:s587-s592, 2013
    # also see K. Kitahara et al., Acta Materialia 54:1279-1288, 2006
    _KS = {'mapping':{'fcc':0,'bcc':1},
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
    _GT = {'mapping':{'fcc':0,'bcc':1},
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
    _GTprime = {'mapping':{'fcc':0,'bcc':1},
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
    _NW = {'mapping':{'fcc':0,'bcc':1},
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
    _Pitsch = {'mapping':{'fcc':0,'bcc':1},
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
    _Bain = {'mapping':{'fcc':0,'bcc':1},
        'planes': np.array([
        [[  1,  0,  0],[  1,  0,  0]],
        [[  0,  1,  0],[  0,  1,  0]],
        [[  0,  0,  1],[  0,  0,  1]]],dtype='float'),
        'directions': np.array([
        [[  0,  1,  0],[  0,  1,  1]],
        [[  0,  0,  1],[  1,  0,  1]],
        [[  1,  0,  0],[  1,  1,  0]]],dtype='float')}


    def relation_operations(self,model):
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
        models={'KS':self._KS, 'GT':self._GT,          'GT_prime':self._GTprime,
                'NW':self._NW, 'Pitsch': self._Pitsch, 'Bain':self._Bain}
        try:
            relationship = models[model]
        except KeyError :
            raise KeyError(f'Orientation relationship "{model}" is unknown')

        if self.lattice not in relationship['mapping']:
            raise ValueError(f'Relationship "{model}" not supported for lattice "{self.lattice}"')

        r = {'lattice':Lattice((set(relationship['mapping'])-{self.lattice}).pop()),                # target lattice
             'rotations':[] }

        myPlane_id    = relationship['mapping'][self.lattice]
        otherPlane_id = (myPlane_id+1)%2
        myDir_id      = myPlane_id    +2
        otherDir_id   = otherPlane_id +2

        for miller in np.hstack((relationship['planes'],relationship['directions'])):
            myPlane     = miller[myPlane_id]/    np.linalg.norm(miller[myPlane_id])
            myDir       = miller[myDir_id]/      np.linalg.norm(miller[myDir_id])
            myMatrix    = np.array([myDir,np.cross(myPlane,myDir),myPlane])

            otherPlane  = miller[otherPlane_id]/ np.linalg.norm(miller[otherPlane_id])
            otherDir    = miller[otherDir_id]/   np.linalg.norm(miller[otherDir_id])
            otherMatrix = np.array([otherDir,np.cross(otherPlane,otherDir),otherPlane])

            r['rotations'].append(np.dot(otherMatrix.T,myMatrix))

        r['rotations'] = np.array(r['rotations'])

        return r
