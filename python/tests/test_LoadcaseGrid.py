import numpy as np
from numpy import ma

from damask import LoadcaseGrid

class TestGridConfig:

    def test_dumper(self):
        a = ma.MaskedArray(np.arange(3,dtype=float),mask=[0,1,0])
        assert str(LoadcaseGrid({'a':a}))             == 'a: [0.0, x, 2.0]\n'
        assert str(LoadcaseGrid({'a':a.astype(int)})) == 'a: [0, x, 2]\n'
        assert str(LoadcaseGrid({'a':a.data}))        == 'a: [0.0, 1.0, 2.0]\n'

    def test_init(self):
        assert LoadcaseGrid() \
            == LoadcaseGrid({'solver':{},
                            'loadstep':[],
                            })
        assert LoadcaseGrid(solver={'mechanical':'spectral_basic'}) \
            == LoadcaseGrid({'solver':{'mechanical':'spectral_basic'},
                             'loadstep':[],
                            })
