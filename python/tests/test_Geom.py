import os
import time

import pytest
import numpy as np

from damask import VTK
from damask import Geom
from damask import Rotation
from damask import util


def geom_equal(a,b):
    return np.all(a.microstructure == b.microstructure) and \
           np.all(a.grid           == b.grid) and \
           np.allclose(a.size, b.size) and \
           str(a.diff(b)) == str(b.diff(a))

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

    @pytest.mark.parametrize('flavor',['plain','explicit'])
    def test_duplicate(self,default,flavor):
        if flavor == 'plain':
            modified = default.duplicate()
        elif flavor == 'explicit':
            modified = default.duplicate(
                            default.microstructure,
                            default.size,
                            default.origin
                           )
        print(modified)
        assert geom_equal(default,modified)

    def test_diff_equal(self,default):
        assert str(default.diff(default)) == ''

    def test_diff_not_equal(self,default):
        new = Geom(default.microstructure[1:,1:,1:]+1,default.size*.9,np.ones(3)-default.origin,comments=['modified'])
        assert str(default.diff(new)) != ''


    def test_set_inplace_outofplace_homogenization(self,default):
        default.set_homogenization(123,inplace=True)
        outofplace = default.set_homogenization(321,inplace=False)
        assert default.homogenization == 123 and outofplace.homogenization == 321


    def test_set_inplace_outofplace_microstructure(self,default):
        default.set_microstructure(np.arange(72).reshape((2,4,9)),inplace=True)
        outofplace = default.set_microstructure(np.arange(72).reshape((8,3,3)),inplace=False)
        assert np.array_equal(default.grid,[2,4,9]) and np.array_equal(outofplace.grid,[8,3,3])


    def test_set_inplace_outofplace_size(self,default):
        default.set_size(np.array([1,2,3]),inplace=True)
        outofplace = default.set_size(np.array([3,2,1]),inplace=False)
        assert np.array_equal(default.size,[1,2,3]) and np.array_equal(outofplace.size,[3,2,1])


    def test_set_inplace_outofplace_comments(self,default):
        default.set_comments(['a','and','b'],inplace=True)
        outofplace = default.set_comments(['b','or','a'],inplace=False)
        assert default.comments == ['a','and','b'] and outofplace.comments == ['b','or','a']


    @pytest.mark.parametrize('masked',[True,False])
    def test_set_microstructure(self,default,masked):
        old = default.microstructure
        new = np.random.randint(200,size=default.grid)
        default.set_microstructure(np.ma.MaskedArray(new,np.full_like(new,masked)),inplace=True)
        assert np.all(default.microstructure==(old if masked else new))


    def test_write_read_str(self,default,tmpdir):
        default.to_file(str(tmpdir/'default.geom'),format='ASCII')
        new = Geom.from_file(str(tmpdir/'default.geom'))
        assert geom_equal(default,new)

    def test_write_read_file(self,default,tmpdir):
        with open(tmpdir/'default.geom','w') as f:
            default.to_file(f,format='ASCII',pack=True)
        with open(tmpdir/'default.geom') as f:
            new = Geom.from_file(f)
        assert geom_equal(default,new)

    def test_write_as_ASCII(self,default,tmpdir):
        with open(tmpdir/'str.geom','w') as f:
            f.write(default.as_ASCII())
        with open(tmpdir/'str.geom') as f:
            new = Geom.from_file(f)
        assert geom_equal(default,new)

    def test_read_write_vtr(self,default,tmpdir):
        default.to_file(tmpdir/'default',format='vtr')
        for _ in range(10):
            time.sleep(.2)
            if os.path.exists(tmpdir/'default.vtr'): break

        new = Geom.from_vtr(tmpdir/'default.vtr')
        assert geom_equal(new,default)

    def test_invalid_geom(self,tmpdir):
        with open('invalid_file','w') as f:
            f.write('this is not a valid header')
        with open('invalid_file','r') as f:
            with pytest.raises(TypeError):
                Geom.from_file(f)

    def test_invalid_vtr(self,tmpdir):
        v = VTK.from_rectilinearGrid(np.random.randint(5,10,3)*2,np.random.random(3) + 1.0)
        v.to_file(tmpdir/'no_materialpoint.vtr')
        for _ in range(10):
            time.sleep(.2)
            if os.path.exists(tmpdir/'no_materialpoint.vtr'): break
        with pytest.raises(ValueError):
            Geom.from_vtr(tmpdir/'no_materialpoint.vtr')


    @pytest.mark.parametrize('pack',[True,False])
    def test_pack(self,default,tmpdir,pack):
        default.to_file(tmpdir/'default.geom',format='ASCII',pack=pack)
        new = Geom.from_file(tmpdir/'default.geom')
        assert geom_equal(new,default)

    def test_invalid_combination(self,default):
        with pytest.raises(ValueError):
            default.duplicate(default.microstructure[1:,1:,1:],size=np.ones(3), autosize=True)

    def test_invalid_size(self,default):
        with pytest.raises(ValueError):
            default.duplicate(default.microstructure[1:,1:,1:],size=np.ones(2))

    def test_invalid_origin(self,default):
        with pytest.raises(ValueError):
            default.duplicate(default.microstructure[1:,1:,1:],origin=np.ones(4))

    def test_invalid_microstructure_size(self,default):
        microstructure = np.ones((3,3))
        with pytest.raises(ValueError):
            default.duplicate(microstructure)

    def test_invalid_microstructure_type(self,default):
        microstructure = np.random.randint(1,300,(3,4,5))==1
        with pytest.raises(TypeError):
            default.duplicate(microstructure)

    def test_invalid_homogenization(self,default):
        with pytest.raises(TypeError):
            default.set_homogenization(homogenization=0)

    def test_invalid_write_format(self,default):
        with pytest.raises(TypeError):
            default.to_file(format='invalid')

    @pytest.mark.parametrize('directions,reflect',[
                                                   (['x'],        False),
                                                   (['x','y','z'],True),
                                                   (['z','x','y'],False),
                                                   (['y','z'],    False)
                                                  ]
                            )
    def test_mirror(self,default,update,reference_dir,directions,reflect):
        modified = default.mirror(directions,reflect)
        tag = f'directions={"-".join(directions)}_reflect={reflect}'
        reference = reference_dir/f'mirror_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(Geom.from_file(reference),
                          modified)

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_mirror_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.mirror(directions)


    @pytest.mark.parametrize('directions',[
                                           ['x'],
                                           ['x','y','z'],
                                           ['z','x','y'],
                                           ['y','z'],
                                          ]
                            )
    def test_flip(self,default,update,reference_dir,directions):
        modified = default.flip(directions)
        tag = f'directions={"-".join(directions)}'
        reference = reference_dir/f'flip_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(Geom.from_file(reference),
                          modified)

    def test_flip_invariant(self,default):
        assert geom_equal(default,default.flip([]))

    @pytest.mark.parametrize('direction',[['x'],['x','y']])
    def test_flip_double(self,default,direction):
        assert geom_equal(default,default.flip(direction).flip(direction))

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_flip_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.flip(directions)


    @pytest.mark.parametrize('stencil',[1,2,3,4])
    @pytest.mark.parametrize('selection',[None,[1],[1,2,3]])
    @pytest.mark.parametrize('periodic',[True,False])
    def test_clean(self,default,update,reference_dir,stencil,selection,periodic):
        current = default.clean(stencil,selection,periodic)
        reference = reference_dir/f'clean_{stencil}_{"+".join(map(str,[None] if selection is None else selection))}_{periodic}'
        if update and stencil > 1:
            current.to_file(reference,format='vtr')
            for _ in range(10):
                time.sleep(.2)
                if os.path.exists(reference.with_suffix('.vtr')): break
        assert geom_equal(Geom.from_vtr(reference) if stencil > 1 else default,
                          current
                         )

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
        modified = default.scale(grid)
        tag = f'grid={util.srepr(grid,"-")}'
        reference = reference_dir/f'scale_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(Geom.from_file(reference),
                          modified)

    def test_renumber(self,default):
        microstructure = default.microstructure.copy()
        for m in np.unique(microstructure):
            microstructure[microstructure==m] = microstructure.max() + np.random.randint(1,30)
        modified = default.duplicate(microstructure)
        assert not geom_equal(modified,default)
        assert geom_equal(default,
                          modified.renumber())

    def test_substitute(self,default):
        offset = np.random.randint(1,500)
        modified = default.duplicate(default.microstructure + offset)
        assert not geom_equal(modified,default)
        assert geom_equal(default,
                          modified.substitute(np.arange(default.microstructure.max())+1+offset,
                                              np.arange(default.microstructure.max())+1))

    @pytest.mark.parametrize('axis_angle',[np.array([1,0,0,86.7]), np.array([0,1,0,90.4]), np.array([0,0,1,90]),
                                           np.array([1,0,0,175]),np.array([0,-1,0,178]),np.array([0,0,1,180])])
    def test_rotate360(self,default,axis_angle):
        modified = default.copy()
        for i in range(np.rint(360/axis_angle[3]).astype(int)):
            modified.rotate(Rotation.from_axis_angle(axis_angle,degrees=True))
        assert geom_equal(default,modified)

    @pytest.mark.parametrize('Eulers',[[32.0,68.0,21.0],
                                       [0.0,32.0,240.0]])
    def test_rotate(self,default,update,reference_dir,Eulers):
        modified = default.rotate(Rotation.from_Eulers(Eulers,degrees=True))
        tag = f'Eulers={util.srepr(Eulers,"-")}'
        reference = reference_dir/f'rotate_{tag}.geom'
        if update: modified.to_file(reference)
        assert geom_equal(Geom.from_file(reference),
                          modified)

    def test_canvas(self,default):
        grid = default.grid
        grid_add = np.random.randint(0,30,(3))
        modified = default.canvas(grid + grid_add)
        assert np.all(modified.microstructure[:grid[0],:grid[1],:grid[2]] == default.microstructure)

    @pytest.mark.parametrize('center1,center2',[(np.random.random(3)*.5,np.random.random()*8),
                                                (np.random.randint(4,8,(3)),np.random.randint(9,12,(3)))])
    @pytest.mark.parametrize('diameter',[np.random.random(3)*.5,
                                         np.random.randint(4,10,(3)),
                                         np.random.rand(),
                                         np.random.randint(30)])
    @pytest.mark.parametrize('exponent',[np.random.random(3)*.5,
                                         np.random.randint(4,10,(3)),
                                         np.random.rand()*4,
                                         np.random.randint(20)])
    def test_add_primitive_shift(self,center1,center2,diameter,exponent):
        """Same volume fraction for periodic microstructures and different center."""
        o = np.random.random(3)-.5
        g = np.random.randint(8,32,(3))
        s = np.random.random(3)+.5
        G_1 = Geom(np.ones(g,'i'),s,o).add_primitive(diameter,center1,exponent)
        G_2 = Geom(np.ones(g,'i'),s,o).add_primitive(diameter,center2,exponent)
        assert np.count_nonzero(G_1.microstructure!=2) == np.count_nonzero(G_2.microstructure!=2)

    @pytest.mark.parametrize('center',[np.random.randint(4,10,(3)),
                                       np.random.randint(2,10),
                                       np.random.rand()*4,
                                       np.random.rand(3)*10])
    @pytest.mark.parametrize('inverse',[True,False])
    @pytest.mark.parametrize('periodic',[True,False])
    def test_add_primitive_rotation(self,center,inverse,periodic):
        """Rotation should not change result for sphere (except for discretization errors)."""
        g = np.array([32,32,32])
        fill = np.random.randint(10)+2
        eu=np.array([np.random.randint(4),np.random.randint(2),np.random.randint(4)])*.5*np.pi
        G_1 = Geom(np.ones(g,'i'),[1.,1.,1.]).add_primitive(.3,center,1,fill,inverse=inverse,periodic=periodic)
        G_2 = Geom(np.ones(g,'i'),[1.,1.,1.]).add_primitive(.3,center,1,fill,Rotation.from_Eulers(eu),inverse,periodic=periodic)
        assert geom_equal(G_1,G_2)

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

        geom = Geom(m,np.random.rand(3)).vicinity_offset(vicinity,offset,trigger=trigger)

        assert np.all(m2==geom.microstructure)

    @pytest.mark.parametrize('periodic',[True,False])
    def test_vicinity_offset_invariant(self,default,periodic):
        offset = default.vicinity_offset(trigger=[default.microstructure.max()+1,
                                                  default.microstructure.min()-1])
        assert np.all(offset.microstructure==default.microstructure)

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
