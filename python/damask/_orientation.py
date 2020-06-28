import numpy as np

from . import Lattice
from . import Rotation

class Orientation: # ToDo: make subclass of lattice and Rotation
    """
    Crystallographic orientation.

    A crystallographic orientation contains a rotation and a lattice.

    """

    __slots__ = ['rotation','lattice']

    def __repr__(self):
        """Report lattice type and orientation."""
        return self.lattice.__repr__()+'\n'+self.rotation.__repr__()

    def __init__(self, rotation, lattice):
        """
        New orientation from rotation and lattice.

        Parameters
        ----------
        rotation : Rotation
            Rotation specifying the lattice orientation.
        lattice : Lattice
            Lattice type of the crystal.

        """
        if isinstance(lattice, Lattice):
            self.lattice = lattice
        else:
            self.lattice = Lattice(lattice)                                                         # assume string

        if isinstance(rotation, Rotation):
            self.rotation = rotation
        else:
            self.rotation = Rotation.from_quaternion(rotation)                                      # assume quaternion


    def disorientation(self,
                       other,
                       SST = True,
                       symmetries = False):
        """
        Disorientation between myself and given other orientation.

        Rotation axis falls into SST if SST == True.

        Currently requires same symmetry for both orientations.
        Look into A. Heinz and P. Neumann 1991 for cases with differing sym.

        """
        if self.lattice.symmetry != other.lattice.symmetry:
            raise NotImplementedError('disorientation between different symmetry classes not supported yet.')

        mySymEqs    =  self.equivalentOrientations() if SST else self.equivalentOrientations([0])   # take all or only first sym operation
        otherSymEqs = other.equivalentOrientations()

        for i,sA in enumerate(mySymEqs):
            aInv = sA.rotation.inversed()
            for j,sB in enumerate(otherSymEqs):
                b = sB.rotation
                r = b*aInv
                for k in range(2):
                    r.inverse()
                    breaker = self.lattice.symmetry.inFZ(r.as_Rodrigues(vector=True)) \
                              and (not SST or other.lattice.symmetry.inDisorientationSST(r.as_Rodrigues(vector=True)))
                    if breaker: break
                if breaker: break
            if breaker: break

        return (Orientation(r,self.lattice), i,j, k == 1) if symmetries else r                      # disorientation ...
                                                                                                    # ... own sym, other sym,
                                                                                                    # self-->other: True, self<--other: False

    def inFZ_vec(self):
        """Check if orientations fall into Fundamental Zone."""
        if not self.rotation.shape:
            return self.lattice.symmetry.inFZ(self.rotation.as_Rodrigues(vector=True))
        else:
            return [self.lattice.symmetry.inFZ(\
                self.rotation.as_Rodrigues(vector=True)[l]) for l in range(self.rotation.shape[0])]


    def inFZ(self):
        return self.lattice.symmetry.inFZ(self.rotation.as_Rodrigues(vector=True))

    @property
    def equivalent_vec(self):
        """
        Return orientations which are symmetrically equivalent.

        One dimension (length according to symmetrically equivalent orientations)
        is added to the left of the rotation array.

        """
        s = self.lattice.symmetry.symmetry_operations #24 lines (sym) x 4 columns (quat)
        s = s.reshape(s.shape[:1]+(1,)*len(self.rotation.shape)+(4,)) #reshape zo (24,1,4)
        s = Rotation(np.broadcast_to(s,s.shape[:1]+self.rotation.quaternion.shape))

        r = np.broadcast_to(self.rotation.quaternion,s.shape[:1]+self.rotation.quaternion.shape) #(24,NumRots,4)
        r = Rotation(r) #(24, NumRot)

        return self.__class__(s@r,self.lattice)


    def equivalentOrientations(self,members=[]):
        """List of orientations which are symmetrically equivalent."""
        try:
            iter(members)                                                                           # asking for (even empty) list of members?
        except TypeError:
            return self.__class__(self.lattice.symmetry.symmetryOperations(members)*self.rotation,self.lattice) # no, return rotation object
        else:
            return [self.__class__(q*self.rotation,self.lattice) \
                                      for q in self.lattice.symmetry.symmetryOperations(members)]   # yes, return list of rotations

    def relatedOrientations_vec(self,model):
        """List of orientations related by the given orientation relationship."""
        h = self.lattice.relationOperations(model)
        rot= h['rotations']
        op=np.array([o.as_quaternion() for o in rot])

        s = op.reshape(op.shape[:1]+(1,)*len(self.rotation.shape)+(4,))
        s = Rotation(np.broadcast_to(s,s.shape[:1]+self.rotation.quaternion.shape))
        
        r = np.broadcast_to(self.rotation.quaternion,s.shape[:1]+self.rotation.quaternion.shape)
        r = Rotation(r)
        
        return self.__class__(s@r,h['lattice'])


    def relatedOrientations(self,model):
        """List of orientations related by the given orientation relationship."""
        r = self.lattice.relationOperations(model)
        return [self.__class__(o*self.rotation,r['lattice']) for o in r['rotations']]

    @property
    def reduced_vec(self):
        """Transform orientation to fall into fundamental zone according to symmetry."""
        equi=self.equivalent_vec.rotation #24 x rot x 3(rodrigues vector)
        r= 1 if not self.rotation.shape else equi.shape[1] #number of rotations
        quat=np.empty( [r , 4])
        for rot in range(r):
            for sym in range(equi.shape[0]):
                if self.lattice.symmetry.inFZ(equi.as_Rodrigues(vector=True)[sym,rot]) == True:
                    quat[rot]=equi.as_quaternion()[sym,rot]
        return self.__class__(quat,self.lattice)



    def reduced(self):
        """Transform orientation to fall into fundamental zone according to symmetry."""
        for me in self.equivalentOrientations():
            if self.lattice.symmetry.inFZ(me.rotation.as_Rodrigues(vector=True)): break

        return self.__class__(me.rotation,self.lattice)


    def inversePole(self,
                    axis,
                    proper = False,
                    SST = True):
        """Axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)."""
        if SST:                                                                                     # pole requested to be within SST
            for i,o in enumerate(self.equivalentOrientations()):                                    # test all symmetric equivalent quaternions
                pole = o.rotation@axis                                                              # align crystal direction to axis
                if self.lattice.symmetry.inSST(pole,proper): break                                  # found SST version
        else:
            pole = self.rotation@axis                                                               # align crystal direction to axis

        return (pole,i if SST else 0)


    def IPFcolor(self,axis):
        """TSL color of inverse pole figure for given axis."""
        color = np.zeros(3,'d')

        for o in self.equivalentOrientations():
            pole = o.rotation@axis                                                                  # align crystal direction to axis
            inSST,color = self.lattice.symmetry.inSST(pole,color=True)
            if inSST: break

        return color


    def IPFcolor_vec(self,axis):
        """TSL color of inverse pole figure for given axis. Not for hex or triclinic lattices"""
        eq = self.equivalent_vec
        pole = eq.rotation @ np.broadcast_to(axis/np.linalg.norm(axis),eq.rotation.shape+(3,))
        in_SST, color = self.lattice.symmetry.in_SST(pole,color=True)

        # ignore duplicates (occur for highly symmetric orientations)
        found = np.zeros_like(in_SST[1],dtype=bool)
        c     = np.empty(color.shape[1:])
        for s in range(in_SST.shape[0]):
            c = np.where(np.expand_dims(np.logical_and(in_SST[s],~found),-1),color[s],c)
            found = np.logical_or(in_SST[s],found)

        return c


    @staticmethod
    def fromAverage(orientations,
                    weights = []):
        """Create orientation from average of list of orientations."""
        # further read: Orientation distribution analysis in deformed grains, https://doi.org/10.1107/S0021889801003077
        if not all(isinstance(item, Orientation) for item in orientations):
            raise TypeError("Only instances of Orientation can be averaged.")

        closest = []
        ref = orientations[0]
        for o in orientations:
            closest.append(o.equivalentOrientations(
                           ref.disorientation(o,
                                              SST = False,                                          # select (o[ther]'s) sym orientation
                                              symmetries = True)[2]).rotation)                      # with lowest misorientation

        return Orientation(Rotation.fromAverage(closest,weights),ref.lattice)


    def average(self,other):
        """Calculate the average rotation."""
        return Orientation.fromAverage([self,other])
