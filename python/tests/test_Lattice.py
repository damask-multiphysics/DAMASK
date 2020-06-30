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
        symmetry = random.choice(Symmetry.crystal_systems)
        print(symmetry)
        assert Symmetry(symmetry) == Symmetry(symmetry)

    def test_not_equal(self):
        symmetries = random.sample(Symmetry.crystal_systems,k=2)
        assert Symmetry(symmetries[0]) != Symmetry(symmetries[1])

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_inFZ(self,system):
          assert Symmetry(system).inFZ(np.zeros(3))

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    def test_inDisorientationSST(self,system):
          assert Symmetry(system).inDisorientationSST(np.zeros(3))

    @pytest.mark.parametrize('system',Symmetry.crystal_systems)
    @pytest.mark.parametrize('proper',[True,False])
    def test_inSST(self,system,proper):
          assert Symmetry(system).inSST(np.zeros(3),proper)

    @pytest.mark.parametrize('function',['inFZ','inDisorientationSST'])
    def test_invalid_argument(self,function):
        s = Symmetry()                                                                              # noqa
        with pytest.raises(ValueError):
            eval('s.{}(np.ones(4))'.format(function))
