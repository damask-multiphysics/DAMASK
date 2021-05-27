import pytest
import numpy as np

from damask import lattice

class TestLattice:

    def test_double_Bravais_to_Miller(self):
        with pytest.raises(KeyError):
            lattice.Bravais_to_Miller(uvtw=np.ones(4),hkil=np.ones(4))                              # noqa

    def test_double_Miller_to_Bravais(self):
        with pytest.raises(KeyError):
            lattice.Miller_to_Bravais(uvw=np.ones(4),hkl=np.ones(4))                                # noqa

    @pytest.mark.parametrize('vector',np.array([
                                                [1,0,0],
                                                [1,1,0],
                                                [1,1,1],
                                                [1,0,-2],
                                               ]))
    @pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
    def test_Miller_Bravais_Miller(self,vector,kw_Miller,kw_Bravais):
        assert np.all(vector == lattice.Bravais_to_Miller(**{kw_Bravais:lattice.Miller_to_Bravais(**{kw_Miller:vector})}))

    @pytest.mark.parametrize('vector',np.array([
                                                [1,0,-1,2],
                                                [1,-1,0,3],
                                                [1,1,-2,-3],
                                                [0,0,0,1],
                                               ]))
    @pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
    def test_Bravais_Miller_Bravais(self,vector,kw_Miller,kw_Bravais):
        assert np.all(vector == lattice.Miller_to_Bravais(**{kw_Miller:lattice.Bravais_to_Miller(**{kw_Bravais:vector})}))
