import numpy as np

from . import Lattice
from . import Rotation

class Orientation:
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
            self.rotation = Rotation.fromQuaternion(rotation)                                       # assume quaternion

    def disorientation(self,
                       other,
                       SST = True,
                       symmetries = False):
        """
        Disorientation between myself and given other orientation.

        Rotation axis falls into SST if SST == True.
        (Currently requires same symmetry for both orientations.
         Look into A. Heinz and P. Neumann 1991 for cases with differing sym.)
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
                    breaker = self.lattice.symmetry.inFZ(r.asRodrigues(vector=True)) \
                              and (not SST or other.lattice.symmetry.inDisorientationSST(r.asRodrigues(vector=True)))
                    if breaker: break
                if breaker: break
            if breaker: break

        return (Orientation(r,self.lattice), i,j, k == 1) if symmetries else r                      # disorientation ...
                                                                                                    # ... own sym, other sym,
                                                                                                    # self-->other: True, self<--other: False
    def inFZ(self):
        return self.lattice.symmetry.inFZ(self.rotation.asRodrigues(vector=True))


    def equivalentOrientations(self,members=[]):
        """List of orientations which are symmetrically equivalent."""
        try:
            iter(members)                                                                           # asking for (even empty) list of members?
        except TypeError:
            return self.__class__(self.lattice.symmetry.symmetryOperations(members)*self.rotation,self.lattice) # no, return rotation object
        else:
            return [self.__class__(q*self.rotation,self.lattice) \
                                      for q in self.lattice.symmetry.symmetryOperations(members)]   # yes, return list of rotations

    def relatedOrientations(self,model):
        """List of orientations related by the given orientation relationship."""
        r = self.lattice.relationOperations(model)
        return [self.__class__(o*self.rotation,r['lattice']) for o in r['rotations']]


    def reduced(self):
        """Transform orientation to fall into fundamental zone according to symmetry."""
        for me in self.equivalentOrientations():
            if self.lattice.symmetry.inFZ(me.rotation.asRodrigues(vector=True)): break

        return self.__class__(me.rotation,self.lattice)


    def inversePole(self,
                    axis,
                    proper = False,
                    SST = True):
        """Axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)."""
        if SST:                                                                                     # pole requested to be within SST
            for i,o in enumerate(self.equivalentOrientations()):                                    # test all symmetric equivalent quaternions
                pole = o.rotation*axis                                                              # align crystal direction to axis
                if self.lattice.symmetry.inSST(pole,proper): break                                  # found SST version
        else:
            pole = self.rotation*axis                                                               # align crystal direction to axis

        return (pole,i if SST else 0)


    def IPFcolor(self,axis):
        """TSL color of inverse pole figure for given axis."""
        color = np.zeros(3,'d')

        for o in self.equivalentOrientations():
            pole = o.rotation*axis                                                                  # align crystal direction to axis
            inSST,color = self.lattice.symmetry.inSST(pole,color=True)
            if inSST: break

        return color


    @staticmethod
    def fromAverage(orientations,
                    weights = []):
        """Create orientation from average of list of orientations."""
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
