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
    """Directory containing reference results."""
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

    @pytest.mark.parametrize('pack',[True,False])
    def test_pack(self,default,tmpdir,pack):
        default.to_file(tmpdir.join('default.geom'),pack=pack)
        new = Geom.from_file(tmpdir.join('default.geom'))
        assert geom_equal(new,default)

    @pytest.mark.parametrize('directions,reflect',[
                                                   (['x'],        False),
                                                   (['x','y','z'],True),
                                                   (['z','x','y'],False),
                                                   (['y','z'],    False)
                                                  ]
                            )
    def test_mirror(self,default,update,reference_dir,directions,reflect):
        modified = copy.deepcopy(default)
        modified.mirror(directions,reflect)
        tag = 'directions={}_reflect={}'.format('-'.join(directions),reflect)
        reference = os.path.join(reference_dir,'mirror_{}.geom'.format(tag))
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))

    @pytest.mark.parametrize('stencil',[(1),(2),(3),(4)])
    def test_clean(self,default,update,reference_dir,stencil):
        modified = copy.deepcopy(default)
        modified.clean(stencil)
        tag = 'stencil={}'.format(stencil)
        reference = os.path.join(reference_dir,'clean_{}.geom'.format(tag))
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))

    @pytest.mark.parametrize('grid',[
                                     ((10,11,10)),
                                     ([10,13,10]),
                                     (np.array((10,10,10))),
                                     (np.array((8, 10,12))),
                                     (np.array((5, 4, 20))),
                                     (np.array((10,20,2)) )
                                    ]
                            )
    def test_scale(self,default,update,reference_dir,grid):
        modified = copy.deepcopy(default)
        modified.scale(grid)
        tag = 'grid={}'.format('-'.join([str(x) for x in grid]))
        reference = os.path.join(reference_dir,'scale_{}.geom'.format(tag))
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))
