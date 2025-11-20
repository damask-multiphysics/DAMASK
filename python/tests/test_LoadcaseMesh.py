# SPDX-License-Identifier: AGPL-3.0-or-later
import numpy as np
from numpy import ma

from damask import LoadcaseMesh

def test_dumper():
    a = ma.MaskedArray(np.arange(3,dtype=float),mask=[0,1,0])
    assert str(LoadcaseMesh({'a':a}))             == 'a: [0.0, x, 2.0]\n'
    assert str(LoadcaseMesh({'a':a.astype(int)})) == 'a: [0, x, 2]\n'
    assert str(LoadcaseMesh({'a':a.data}))        == 'a: [0.0, 1.0, 2.0]\n'

def test_init():
    assert LoadcaseMesh() \
        == LoadcaseMesh({
                         'loadstep':[],
                        })
