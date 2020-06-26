import random

import pytest
import numpy as np

from damask import Symmetry

class TestSymmetry:

    @pytest.mark.parametrize('invalid_symmetry',['fcc','bcc','hello'])
    def test_invalid_symmetry(self,invalid_symmetry):
        with pytest.raises(KeyError):
            s = Symmetry(invalid_symmetry)                                                          # noqa

    def test_equal(self):
        symmetry = random.choice(Symmetry.lattices)
        print(symmetry)
        assert Symmetry(symmetry) == Symmetry(symmetry)

    def test_not_equal(self):
        symmetries = random.sample(Symmetry.lattices,k=2)
        assert Symmetry(symmetries[0]) != Symmetry(symmetries[1])

    @pytest.mark.parametrize('lattice',Symmetry.lattices)
    def test_inFZ(self,lattice):
          assert Symmetry(lattice).inFZ(np.zeros(3))

    @pytest.mark.parametrize('lattice',Symmetry.lattices)
    def test_inDisorientationSST(self,lattice):
          assert Symmetry(lattice).inDisorientationSST(np.zeros(3))

    @pytest.mark.parametrize('lattice',Symmetry.lattices)
    @pytest.mark.parametrize('proper',[True,False])
    def test_inSST(self,lattice,proper):
          assert Symmetry(lattice).inSST(np.zeros(3),proper)

    @pytest.mark.parametrize('function',['inFZ','inDisorientationSST'])
    def test_invalid_argument(self,function):
        s = Symmetry()                                                                              # noqa
        with pytest.raises(ValueError):
            eval(f's.{function}(np.ones(4))')
