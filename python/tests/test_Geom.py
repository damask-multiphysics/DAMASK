import copy
import os

import pytest
import numpy as np

from damask import Geom


def geom_equal(a,b):
    return np.all(a.get_microstructure() == b.get_microstructure()) and \
           np.all(a.get_size()           == b.get_size())           and \
           np.all(a.get_grid()           == b.get_grid())

@pytest.fixture
def default():
    """Simple geometry."""
    x=np.concatenate((np.ones(40,dtype=int),
                      np.arange(2,42),
                      np.ones(40,dtype=int)*2,
                      np.arange(1,41))).reshape((8,5,4))
    return Geom(x,[8e-6,5e-6,4e-6])

@pytest.fixture
def reference_dir(reference_dir_base):
    """directory containing reference results."""
    return os.path.join(reference_dir_base,'Geom')


class TestGeom:
    
    def test_update(self,default):
        modified = copy.deepcopy(default)
        modified.update(
                        default.get_microstructure(),
                        default.get_size(),
                        default.get_origin()
                       )
        assert geom_equal(modified,default)


    def test_write_read_str(self,default,tmpdir):
        default.to_file(str(tmpdir.join('default.geom')))
        new = Geom.from_file(str(tmpdir.join('default.geom')))
        assert geom_equal(new,default)

    def test_write_read_file(self,default,tmpdir):
        with open(tmpdir.join('default.geom'),'w') as f:
            default.to_file(f)
        with open(tmpdir.join('default.geom')) as f:
            new = Geom.from_file(f)
        assert geom_equal(new,default)

    def test_mirror(self,default,update,reference_dir):
        modified = copy.deepcopy(default)
        modified.mirror(['x','z']) 
        reference = os.path.join(reference_dir,'mirror.geom')
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))
