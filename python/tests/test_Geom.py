import os
import time

import pytest
import numpy as np

from damask import VTK
from damask import Geom
from damask import Rotation
from damask import util


def geom_equal(a,b):
    return np.all(a.get_microstructure() == b.get_microstructure()) and \
           np.all(a.get_grid()           == b.get_grid()) and \
           np.allclose(a.get_size(), b.get_size())

@pytest.fixture
def default():
    """Simple geometry."""
    x=np.concatenate((np.ones(40,dtype=int),
                      np.arange(2,42),
                      np.ones(40,dtype=int)*2,
                      np.arange(1,41))).reshape(8,5,4,order='F')
    return Geom(x,[8e-6,5e-6,4e-6])

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return reference_dir_base/'Geom'


class TestGeom:

    def test_update(self,default):
        modified = default.copy()
        modified.update(
                        default.get_microstructure(),
                        default.get_size(),
                        default.get_origin()
                       )
        print(modified)
        assert geom_equal(modified,default)

    @pytest.mark.parametrize('masked',[True,False])
    def test_set_microstructure(self,default,masked):
        old = default.get_microstructure()
        new = np.random.randint(200,size=default.grid)
        default.set_microstructure(np.ma.MaskedArray(new,np.full_like(new,masked)))
        assert np.all(default.microstructure==(old if masked else new))


    def test_write_read_str(self,default,tmpdir):
        default.to_file(str(tmpdir/'default.geom'))
        new = Geom.from_file(str(tmpdir/'default.geom'))
        assert geom_equal(new,default)

    def test_write_read_file(self,default,tmpdir):
        with open(tmpdir/'default.geom','w') as f:
            default.to_file(f)
        with open(tmpdir/'default.geom') as f:
            new = Geom.from_file(f)
        assert geom_equal(new,default)

    def test_write_show(self,default,tmpdir):
        with open(tmpdir/'str.geom','w') as f:
            f.write(default.show())
        with open(tmpdir/'str.geom') as f:
            new = Geom.from_file(f)
        assert geom_equal(new,default)

    def test_read_write_vtr(self,default,tmpdir):
        default.to_vtr(tmpdir/'default')
        for _ in range(10):
            time.sleep(.2)
            if os.path.exists(tmpdir/'default.vtr'): break

        new = Geom.from_vtr(tmpdir/'default.vtr')
        assert geom_equal(new,default)

    def test_invalid_vtr(self,tmpdir):
        v = VTK.from_rectilinearGrid(np.random.randint(5,10,3)*2,np.random.random(3) + 1.0)
        v.write(tmpdir/'no_materialpoint.vtr')
        for _ in range(10):
            time.sleep(.2)
            if os.path.exists(tmpdir/'no_materialpoint.vtr'): break
        with pytest.raises(ValueError):
            Geom.from_vtr(tmpdir/'no_materialpoint.vtr')


    @pytest.mark.parametrize('pack',[True,False])
    def test_pack(self,default,tmpdir,pack):
        default.to_file(tmpdir/'default.geom',pack=pack)
        new = Geom.from_file(tmpdir/'default.geom')
        assert geom_equal(new,default)

    def test_invalid_combination(self,default):
        with pytest.raises(ValueError):
            default.update(default.microstructure[1:,1:,1:],size=np.ones(3), rescale=True)

    def test_invalid_size(self,default):
        with pytest.raises(ValueError):
            default.update(default.microstructure[1:,1:,1:],size=np.ones(2))

    def test_invalid_origin(self,default):
        with pytest.raises(ValueError):
            default.update(default.microstructure[1:,1:,1:],origin=np.ones(4))

    def test_invalid_microstructure_size(self,default):
        microstructure = np.ones((3,3))
        with pytest.raises(ValueError):
            default.update(microstructure)

    def test_invalid_microstructure_type(self,default):
        microstructure = np.random.randint(1,300,(3,4,5))==1
        with pytest.raises(TypeError):
            default.update(microstructure)

    def test_invalid_homogenization(self,default):
        with pytest.raises(TypeError):
            default.set_homogenization(homogenization=0)

    @pytest.mark.parametrize('directions,reflect',[
                                                   (['x'],        False),
                                                   (['x','y','z'],True),
                                                   (['z','x','y'],False),
                                                   (['y','z'],    False)
                                                  ]
                            )
    def test_mirror(self,default,update,reference_dir,directions,reflect):
        modified = default.copy()
        modified.mirror(directions,reflect)
        tag = f'directions={"-".join(directions)}_reflect={reflect}'
        reference = reference_dir/f'mirror_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_mirror_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.mirror(directions)

    @pytest.mark.parametrize('stencil',[1,2,3,4])
    @pytest.mark.parametrize('selection',[None,1,2])
    @pytest.mark.parametrize('periodic',[True,False])
    def test_clean(self,update,reference_dir,stencil,selection,periodic):
        current = Geom.from_vtr((reference_dir/'clean').with_suffix('.vtr'))
        current.clean(stencil,None if selection is None else [selection],periodic)
        reference = reference_dir/f'clean_{stencil}_{selection}_{periodic}'
        if update and stencil !=1:
            current.to_vtr(reference)
            for _ in range(10):
                time.sleep(.2)
                if os.path.exists(reference.with_suffix('.vtr')): break
        assert geom_equal(current,Geom.from_vtr(reference if stencil !=1 else reference_dir/'clean'))

    @pytest.mark.parametrize('grid',[
                                     (10,11,10),
                                     [10,13,10],
                                     np.array((10,10,10)),
                                     np.array((8, 10,12)),
                                     np.array((5, 4, 20)),
                                     np.array((10,20,2))
                                    ]
                            )
    def test_scale(self,default,update,reference_dir,grid):
        modified = default.copy()
        modified.scale(grid)
        tag = f'grid={util.srepr(grid,"-")}'
        reference = reference_dir/f'scale_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))

    def test_renumber(self,default):
        modified = default.copy()
        microstructure = modified.get_microstructure()
        for m in np.unique(microstructure):
            microstructure[microstructure==m] = microstructure.max() + np.random.randint(1,30)
        modified.update(microstructure)
        assert not geom_equal(modified,default)
        modified.renumber()
        assert geom_equal(modified,default)

    def test_substitute(self,default):
        modified = default.copy()
        microstructure = modified.get_microstructure()
        offset = np.random.randint(1,500)
        microstructure += offset
        modified.update(microstructure)
        assert not geom_equal(modified,default)
        modified.substitute(np.arange(default.microstructure.max())+1+offset,
                            np.arange(default.microstructure.max())+1)
        assert geom_equal(modified,default)

    @pytest.mark.parametrize('axis_angle',[np.array([1,0,0,86.7]), np.array([0,1,0,90.4]), np.array([0,0,1,90]),
                                           np.array([1,0,0,175]),np.array([0,-1,0,178]),np.array([0,0,1,180])])
    def test_rotate360(self,default,axis_angle):
        modified = default.copy()
        for i in range(np.rint(360/axis_angle[3]).astype(int)):
            modified.rotate(Rotation.from_axis_angle(axis_angle,degrees=True))
        assert geom_equal(modified,default)

    @pytest.mark.parametrize('Eulers',[[32.0,68.0,21.0],
                                       [0.0,32.0,240.0]])
    def test_rotate(self,default,update,reference_dir,Eulers):
        modified = default.copy()
        modified.rotate(Rotation.from_Eulers(Eulers,degrees=True))
        tag = f'Eulers={util.srepr(Eulers,"-")}'
        reference = reference_dir/f'rotate_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(modified,Geom.from_file(reference))

    def test_canvas(self,default):
        grid_add = np.random.randint(0,30,(3))
        modified = default.copy()
        modified.canvas(modified.grid + grid_add)
        e = default.grid
        assert np.all(modified.microstructure[:e[0],:e[1],:e[2]] == default.microstructure)

    @pytest.mark.parametrize('center1,center2',[(np.random.random(3)*.5,np.random.random(3)),
                                                (np.random.randint(4,8,(3)),np.random.randint(9,12,(3)))])
    @pytest.mark.parametrize('diameter',[np.random.random(3)*.5,
                                        np.random.randint(4,10,(3))])
    def test_add_primitive(self,diameter,center1,center2):
        """Same volume fraction for periodic microstructures and different center."""
        o = np.random.random(3)-.5
        g = np.random.randint(8,32,(3))
        s = np.random.random(3)+.5
        G_1 = Geom(np.ones(g,'i'),s,o)
        G_2 = Geom(np.ones(g,'i'),s,o)
        G_1.add_primitive(diameter,center1,1)
        G_2.add_primitive(diameter,center2,1)
        assert np.count_nonzero(G_1.microstructure!=2) == np.count_nonzero(G_2.microstructure!=2)

    @pytest.mark.parametrize('trigger',[[1],[]])
    def test_vicinity_offset(self,trigger):
        offset = np.random.randint(2,4)
        vicinity = np.random.randint(2,4)

        g = np.random.randint(28,40,(3))
        m = np.ones(g,'i')
        x = (g*np.random.permutation(np.array([.5,1,1]))).astype('i')
        m[slice(0,x[0]),slice(0,x[1]),slice(0,x[2])] = 2
        m2 = m.copy()
        for i in [0,1,2]:
            m2[(np.roll(m,+vicinity,i)-m)!=0] += offset
            m2[(np.roll(m,-vicinity,i)-m)!=0] += offset
        if len(trigger) > 0:
            m2[m==1] = 1

        geom = Geom(m,np.random.rand(3))
        geom.vicinity_offset(vicinity,offset,trigger=trigger)

        assert np.all(m2==geom.microstructure)

    @pytest.mark.parametrize('periodic',[True,False])
    def test_vicinity_offset_invariant(self,default,periodic):
        old = default.get_microstructure()
        default.vicinity_offset(trigger=[old.max()+1,old.min()-1])
        assert np.all(old==default.microstructure)

    @pytest.mark.parametrize('periodic',[True,False])
    def test_tessellation_approaches(self,periodic):
        grid   = np.random.randint(10,20,3)
        size   = np.random.random(3) + 1.0
        N_seeds= np.random.randint(10,30)
        seeds  = np.random.rand(N_seeds,3) * np.broadcast_to(size,(N_seeds,3))
        Voronoi  = Geom.from_Voronoi_tessellation( grid,size,seeds,                 periodic)
        Laguerre = Geom.from_Laguerre_tessellation(grid,size,seeds,np.ones(N_seeds),periodic)
        assert geom_equal(Laguerre,Voronoi)

    def test_Laguerre_weights(self):
        grid   = np.random.randint(10,20,3)
        size   = np.random.random(3) + 1.0
        N_seeds= np.random.randint(10,30)
        seeds  = np.random.rand(N_seeds,3) * np.broadcast_to(size,(N_seeds,3))
        weights= np.full((N_seeds),-np.inf)
        ms     = np.random.randint(1, N_seeds+1)
        weights[ms-1] = np.random.random()
        Laguerre = Geom.from_Laguerre_tessellation(grid,size,seeds,weights,np.random.random()>0.5)
        assert np.all(Laguerre.microstructure == ms)

    @pytest.mark.parametrize('approach',['Laguerre','Voronoi'])
    def test_tessellate_bicrystal(self,approach):
        grid  = np.random.randint(5,10,3)*2
        size  = grid.astype(np.float)
        seeds = np.vstack((size*np.array([0.5,0.25,0.5]),size*np.array([0.5,0.75,0.5])))
        microstructure = np.ones(grid)
        microstructure[:,grid[1]//2:,:] = 2
        if   approach == 'Laguerre':
            geom = Geom.from_Laguerre_tessellation(grid,size,seeds,np.ones(2),np.random.random()>0.5)
        elif approach == 'Voronoi':
            geom = Geom.from_Voronoi_tessellation(grid,size,seeds,            np.random.random()>0.5)
        assert np.all(geom.microstructure == microstructure)
