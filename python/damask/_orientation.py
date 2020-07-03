import numpy as np

from . import Lattice
from . import Rotation

class Orientation: # ToDo: make subclass of lattice and Rotation?
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

    def __getitem__(self,item):
        """Iterate over leading/leftmost dimension of Orientation array."""
        return self.__class__(self.rotation[item],self.lattice)


    # ToDo: Discuss vectorization/calling signature
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

        mySymEqs    =  self.equivalent if SST else self.equivalent[0] #ToDo: This is just me!       # take all or only first sym operation
        otherSymEqs = other.equivalent

        for i,sA in enumerate(mySymEqs):
            aInv = sA.rotation.inversed()
            for j,sB in enumerate(otherSymEqs):
                b = sB.rotation
                r = b*aInv
                for k in range(2):
                    r.inverse()
                    breaker = self.lattice.in_FZ(r.as_Rodrigues(vector=True)) \
                              and (not SST or other.lattice.in_disorientation_SST(r.as_Rodrigues(vector=True)))
                    if breaker: break
                if breaker: break
            if breaker: break

        return (Orientation(r,self.lattice), i,j, k == 1) if symmetries else r                      # disorientation ...
                                                                                                    # ... own sym, other sym,
                                                                                                    # self-->other: True, self<--other: False

    @property
    def in_FZ(self):
        """Check if orientations fall into Fundamental Zone."""
        return self.lattice.in_FZ(self.rotation.as_Rodrigues(vector=True))


    @property
    def equivalent(self):
        """
        Orientations which are symmetrically equivalent.

        One dimension (length according to number of symmetrically equivalent orientations)
        is added to the left of the Rotation array.

        """
        o = self.lattice.symmetry.symmetry_operations
        o = o.reshape(o.shape[:1]+(1,)*len(self.rotation.shape)+(4,))
        o = Rotation(np.broadcast_to(o,o.shape[:1]+self.rotation.quaternion.shape))

        s = np.broadcast_to(self.rotation.quaternion,o.shape[:1]+self.rotation.quaternion.shape)

        return self.__class__(o@Rotation(s),self.lattice)


    def related(self,model):
        """
        Orientations related by the given orientation relationship.

        One dimension (length according to number of related orientations)
        is added to the left of the Rotation array.

        """
        o = Rotation.from_matrix(self.lattice.relation_operations(model)['rotations']).as_quaternion()
        o = o.reshape(o.shape[:1]+(1,)*len(self.rotation.shape)+(4,))
        o = Rotation(np.broadcast_to(o,o.shape[:1]+self.rotation.quaternion.shape))

        s = np.broadcast_to(self.rotation.quaternion,o.shape[:1]+self.rotation.quaternion.shape)

        return self.__class__(o@Rotation(s),self.lattice.relation_operations(model)['lattice'])


    @property
    def reduced(self):
        """Transform orientation to fall into fundamental zone according to symmetry."""
        eq = self.equivalent
        in_FZ = eq.in_FZ

        # remove duplicates (occur for highly symmetric orientations)
        found = np.zeros_like(in_FZ[0],dtype=bool)
        q     = self.rotation.quaternion[0]
        for s in range(in_FZ.shape[0]):
            #something fishy... why does q needs to be initialized?
            q = np.where(np.expand_dims(np.logical_and(in_FZ[s],~found),-1),eq.rotation.quaternion[s],q)
            found = np.logical_or(in_FZ[s],found)

        return self.__class__(q,self.lattice)


    def inverse_pole(self,axis,proper=False,SST=True):
        """Axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)."""
        if SST:
            eq = self.equivalent
            pole = eq.rotation @ np.broadcast_to(axis/np.linalg.norm(axis),eq.rotation.shape+(3,))
            in_SST = self.lattice.in_SST(pole,proper=proper)

            # remove duplicates (occur for highly symmetric orientations)
            found = np.zeros_like(in_SST[0],dtype=bool)
            p     = pole[0]
            for s in range(in_SST.shape[0]):
                p = np.where(np.expand_dims(np.logical_and(in_SST[s],~found),-1),pole[s],p)
                found = np.logical_or(in_SST[s],found)

            return p
        else:
            return self.rotation @ np.broadcast_to(axis/np.linalg.norm(axis),self.rotation.shape+(3,))



    def IPF_color(self,axis): #ToDo axis or direction?
        """TSL color of inverse pole figure for given axis."""
        eq = self.equivalent
        pole = eq.rotation @ np.broadcast_to(axis/np.linalg.norm(axis),eq.rotation.shape+(3,))
        in_SST, color = self.lattice.in_SST(pole,color=True)

        # remove duplicates (occur for highly symmetric orientations)
        found = np.zeros_like(in_SST[0],dtype=bool)
        c     = color[0]
        for s in range(in_SST.shape[0]):
            c = np.where(np.expand_dims(np.logical_and(in_SST[s],~found),-1),color[s],c)
            found = np.logical_or(in_SST[s],found)

        return c


    # ToDo: Discuss vectorization/calling signature
    @staticmethod
    def from_average(orientations,
                    weights = []):
        """Create orientation from average of list of orientations."""
        # further read: Orientation distribution analysis in deformed grains
        # https://doi.org/10.1107/S0021889801003077
        if not all(isinstance(item, Orientation) for item in orientations):
            raise TypeError("Only instances of Orientation can be averaged.")

        closest = []
        ref = orientations[0]
        for o in orientations:
            closest.append(o.equivalent[
                           ref.disorientation(o,
                                              SST = False,                                          # select (o[ther]'s) sym orientation
                                              symmetries = True)[2]].rotation)                      # with lowest misorientation

        return Orientation(Rotation.from_average(closest,weights),ref.lattice)


    # ToDo: Discuss vectorization/calling signature
    def average(self,other):
        """Calculate the average rotation."""
        return Orientation.from_average([self,other])
