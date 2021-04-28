import pytest
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk as np_to_vtk

from damask import VTK
from damask import Grid
from damask import Table
from damask import Rotation
from damask import util
from damask import seeds
from damask import grid_filters


def grid_equal(a,b):
    return np.all(a.material == b.material) and \
           np.all(a.cells    == b.cells) and \
           np.allclose(a.size, b.size)

@pytest.fixture
def default():
    """Simple geometry."""
    x=np.concatenate((np.ones(40,dtype=int),
                      np.arange(2,42),
                      np.ones(40,dtype=int)*2,
                      np.arange(1,41))).reshape(8,5,4,order='F')
    return Grid(x,[8e-6,5e-6,4e-6])

@pytest.fixture
def ref_path(ref_path_base):
    """Directory containing reference results."""
    return ref_path_base/'Grid'


class TestGrid:

    @pytest.fixture(autouse=True)
    def _patch_execution_stamp(self, patch_execution_stamp):
        print('patched damask.util.execution_stamp')

    @pytest.fixture(autouse=True)
    def _patch_datetime_now(self, patch_datetime_now):
        print('patched datetime.datetime.now')


    def test_equal(self,default):
        assert default == default

    def test_repr(self,default):
        print(default)

    def test_read_write_vtr(self,default,tmp_path):
        default.save(tmp_path/'default')
        new = Grid.load(tmp_path/'default.vtr')
        assert grid_equal(new,default)

    def test_invalid_no_material(self,tmp_path):
        v = VTK.from_rectilinear_grid(np.random.randint(5,10,3)*2,np.random.random(3) + 1.0)
        v.save(tmp_path/'no_materialpoint.vtr',parallel=False)
        with pytest.raises(ValueError):
            Grid.load(tmp_path/'no_materialpoint.vtr')

    def test_invalid_spacing(self,tmp_path,default):
        default.save(tmp_path/'spacing_ok.vtr')
        vtk = VTK.load(tmp_path/'spacing_ok.vtr')
        vtk.vtk_data.SetXCoordinates(np_to_vtk(np.sort(np.random.random(default.cells[0]))))
        vtk.save(tmp_path/'invalid_spacing.vtr',parallel=False)
        with pytest.raises(ValueError):
            Grid.load(tmp_path/'invalid_spacing.vtr')

    def test_invalid_material_type(self):
        with pytest.raises(TypeError):
            Grid(np.zeros((3,3,3),dtype='complex'),np.ones(3))

    def test_cast_to_int(self):
        g = Grid(np.zeros((3,3,3)),np.ones(3))
        assert g.material.dtype in np.sctypes['int']

    def test_invalid_size(self,default):
        with pytest.raises(ValueError):
            Grid(default.material[1:,1:,1:],
                 size=np.ones(2))

    def test_save_load_ASCII(self,default,tmp_path):
        default.save_ASCII(tmp_path/'ASCII')
        default.material -= 1
        assert grid_equal(Grid.load_ASCII(tmp_path/'ASCII'),default)

    def test_invalid_origin(self,default):
        with pytest.raises(ValueError):
            Grid(default.material[1:,1:,1:],
                 size=np.ones(3),
                 origin=np.ones(4))


    def test_invalid_materials_shape(self,default):
        material = np.ones((3,3))
        with pytest.raises(ValueError):
            Grid(material,
                 size=np.ones(3))


    def test_invalid_materials_type(self,default):
        material = np.random.randint(1,300,(3,4,5))==1
        with pytest.raises(TypeError):
            Grid(material)


    @pytest.mark.parametrize('directions,reflect',[
                                                   (['x'],        False),
                                                   (['x','y','z'],True),
                                                   (['z','x','y'],False),
                                                   (['y','z'],    False)
                                                  ]
                            )
    def test_mirror(self,default,update,ref_path,directions,reflect):
        modified = default.mirror(directions,reflect)
        tag = f'directions_{"-".join(directions)}+reflect_{reflect}'
        reference = ref_path/f'mirror_{tag}.vtr'
        if update: modified.save(reference)
        assert grid_equal(Grid.load(reference),
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
    def test_flip(self,default,update,ref_path,directions):
        modified = default.flip(directions)
        tag = f'directions_{"-".join(directions)}'
        reference = ref_path/f'flip_{tag}.vtr'
        if update: modified.save(reference)
        assert grid_equal(Grid.load(reference),
                          modified)


    def test_flip_invariant(self,default):
        assert grid_equal(default,default.flip([]))


    @pytest.mark.parametrize('direction',[['x'],['x','y']])
    def test_flip_double(self,default,direction):
        assert grid_equal(default,default.flip(direction).flip(direction))


    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_flip_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.flip(directions)


    @pytest.mark.parametrize('stencil',[1,2,3,4])
    @pytest.mark.parametrize('selection',[None,[1],[1,2,3]])
    @pytest.mark.parametrize('periodic',[True,False])
    def test_clean(self,default,update,ref_path,stencil,selection,periodic):
        current = default.clean(stencil,selection,periodic)
        reference = ref_path/f'clean_{stencil}_{"+".join(map(str,[None] if selection is None else selection))}_{periodic}'
        if update and stencil > 1:
            current.save(reference)
        assert grid_equal(Grid.load(reference) if stencil > 1 else default,
                          current
                         )


    @pytest.mark.parametrize('cells',[
                                     (10,11,10),
                                     [10,13,10],
                                     np.array((10,10,10)),
                                     np.array((8, 10,12)),
                                     np.array((5, 4, 20)),
                                     np.array((10,20,2))
                                    ]
                            )
    def test_scale(self,default,update,ref_path,cells):
        modified = default.scale(cells)
        tag = f'grid_{util.srepr(cells,"-")}'
        reference = ref_path/f'scale_{tag}.vtr'
        if update: modified.save(reference)
        assert grid_equal(Grid.load(reference),
                          modified)


    def test_renumber(self,default):
        material = default.material.copy()
        for m in np.unique(material):
            material[material==m] = material.max() + np.random.randint(1,30)
        default.material -= 1
        modified = Grid(material,
                        default.size,
                        default.origin)
        assert not grid_equal(modified,default)
        assert grid_equal(default,
                          modified.renumber())


    def test_substitute(self,default):
        offset = np.random.randint(1,500)
        modified = Grid(default.material + offset,
                        default.size,
                        default.origin)
        assert not grid_equal(modified,default)
        assert grid_equal(default,
                          modified.substitute(np.arange(default.material.max())+1+offset,
                                              np.arange(default.material.max())+1))

    def test_substitute_invariant(self,default):
        f = np.unique(default.material.flatten())[:np.random.randint(1,default.material.max())]
        t = np.random.permutation(f)
        modified = default.substitute(f,t)
        assert np.array_equiv(t,f) or (not grid_equal(modified,default))
        assert grid_equal(default, modified.substitute(t,f))

    def test_sort(self):
        cells = np.random.randint(5,20,3)
        m = Grid(np.random.randint(1,20,cells)*3,np.ones(3)).sort().material.flatten(order='F')
        for i,v in enumerate(m):
            assert i==0 or v > m[:i].max() or v in m[:i]

    @pytest.mark.parametrize('axis_angle',[np.array([1,0,0,86.7]), np.array([0,1,0,90.4]), np.array([0,0,1,90]),
                                           np.array([1,0,0,175]),np.array([0,-1,0,178]),np.array([0,0,1,180])])
    def test_rotate360(self,default,axis_angle):
        modified = default.copy()
        for i in range(np.rint(360/axis_angle[3]).astype(int)):
            modified.rotate(Rotation.from_axis_angle(axis_angle,degrees=True))
        assert grid_equal(default,modified)


    @pytest.mark.parametrize('Eulers',[[32.0,68.0,21.0],
                                       [0.0,32.0,240.0]])
    def test_rotate(self,default,update,ref_path,Eulers):
        modified = default.rotate(Rotation.from_Euler_angles(Eulers,degrees=True))
        tag = f'Eulers_{util.srepr(Eulers,"-")}'
        reference = ref_path/f'rotate_{tag}.vtr'
        if update: modified.save(reference)
        assert grid_equal(Grid.load(reference),
                          modified)


    def test_canvas(self,default):
        cells = default.cells
        grid_add = np.random.randint(0,30,(3))
        modified = default.canvas(cells + grid_add)
        assert np.all(modified.material[:cells[0],:cells[1],:cells[2]] == default.material)


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
        """Same volume fraction for periodic geometries and different center."""
        o = np.random.random(3)-.5
        g = np.random.randint(8,32,(3))
        s = np.random.random(3)+.5
        G_1 = Grid(np.ones(g,'i'),s,o).add_primitive(diameter,center1,exponent)
        G_2 = Grid(np.ones(g,'i'),s,o).add_primitive(diameter,center2,exponent)
        assert np.count_nonzero(G_1.material!=2) == np.count_nonzero(G_2.material!=2)


    @pytest.mark.parametrize('center',[np.random.randint(4,10,(3)),
                                       np.random.randint(2,10),
                                       np.random.rand()*4,
                                       np.random.rand(3)*10])
    @pytest.mark.parametrize('inverse',[True,False])
    @pytest.mark.parametrize('periodic',[True,False])
    def test_add_primitive_rotation(self,center,inverse,periodic):
        """Rotation should not change result for sphere."""
        g = np.random.randint(8,32,(3))
        s = np.random.random(3)+.5
        fill = np.random.randint(10)+2
        G_1 = Grid(np.ones(g,'i'),s).add_primitive(.3,center,1,fill,inverse=inverse,periodic=periodic)
        G_2 = Grid(np.ones(g,'i'),s).add_primitive(.3,center,1,fill,Rotation.from_random(),inverse,periodic=periodic)
        assert grid_equal(G_1,G_2)


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

        grid = Grid(m,np.random.rand(3)).vicinity_offset(vicinity,offset,trigger=trigger)

        assert np.all(m2==grid.material)


    @pytest.mark.parametrize('periodic',[True,False])
    def test_vicinity_offset_invariant(self,default,periodic):
        offset = default.vicinity_offset(trigger=[default.material.max()+1,
                                                  default.material.min()-1])
        assert np.all(offset.material==default.material)


    @pytest.mark.parametrize('periodic',[True,False])
    def test_tessellation_approaches(self,periodic):
        cells  = np.random.randint(10,20,3)
        size   = np.random.random(3) + 1.0
        N_seeds= np.random.randint(10,30)
        seeds  = np.random.rand(N_seeds,3) * np.broadcast_to(size,(N_seeds,3))
        Voronoi  = Grid.from_Voronoi_tessellation( cells,size,seeds,                 np.arange(N_seeds)+5,periodic)
        Laguerre = Grid.from_Laguerre_tessellation(cells,size,seeds,np.ones(N_seeds),np.arange(N_seeds)+5,periodic)
        assert grid_equal(Laguerre,Voronoi)


    def test_Laguerre_weights(self):
        cells  = np.random.randint(10,20,3)
        size   = np.random.random(3) + 1.0
        N_seeds= np.random.randint(10,30)
        seeds  = np.random.rand(N_seeds,3) * np.broadcast_to(size,(N_seeds,3))
        weights= np.full((N_seeds),-np.inf)
        ms     = np.random.randint(N_seeds)
        weights[ms] = np.random.random()
        Laguerre = Grid.from_Laguerre_tessellation(cells,size,seeds,weights,periodic=np.random.random()>0.5)
        assert np.all(Laguerre.material == ms)


    @pytest.mark.parametrize('approach',['Laguerre','Voronoi'])
    def test_tessellate_bicrystal(self,approach):
        cells = np.random.randint(5,10,3)*2
        size  = cells.astype(float)
        seeds = np.vstack((size*np.array([0.5,0.25,0.5]),size*np.array([0.5,0.75,0.5])))
        material = np.zeros(cells)
        material[:,cells[1]//2:,:] = 1
        if   approach == 'Laguerre':
            grid = Grid.from_Laguerre_tessellation(cells,size,seeds,np.ones(2),periodic=np.random.random()>0.5)
        elif approach == 'Voronoi':
            grid = Grid.from_Voronoi_tessellation(cells,size,seeds,            periodic=np.random.random()>0.5)
        assert np.all(grid.material == material)


    @pytest.mark.parametrize('surface',['Schwarz P',
                                        'Double Primitive',
                                        'Schwarz D',
                                        'Complementary D',
                                        'Double Diamond',
                                        'Dprime',
                                        'Gyroid',
                                        'Gprime',
                                        'Karcher K',
                                        'Lidinoid',
                                        'Neovius',
                                        'Fisher-Koch S',
                                        ])
    def test_minimal_surface_basic_properties(self,surface):
        cells = np.random.randint(60,100,3)
        size  = np.ones(3)+np.random.rand(3)
        threshold = 2*np.random.rand()-1.
        periods = np.random.randint(2)+1
        materials = np.random.randint(0,40,2)
        grid = Grid.from_minimal_surface(cells,size,surface,threshold,periods,materials)
        assert set(grid.material.flatten()) | set(materials) == set(materials) \
               and (grid.size == size).all() and (grid.cells == cells).all()

    @pytest.mark.parametrize('surface,threshold',[('Schwarz P',0),
                                        ('Double Primitive',-1./6.),
                                        ('Schwarz D',0),
                                        ('Complementary D',0),
                                        ('Double Diamond',-0.133),
                                        ('Dprime',-0.0395),
                                        ('Gyroid',0),
                                        ('Gprime',0.22913),
                                        ('Karcher K',0.17045),
                                        ('Lidinoid',0.14455),
                                        ('Neovius',0),
                                        ('Fisher-Koch S',0),
                                        ])
    def test_minimal_surface_volume(self,surface,threshold):
        cells = np.ones(3,dtype=int)*64
        grid = Grid.from_minimal_surface(cells,np.ones(3),surface,threshold)
        assert np.isclose(np.count_nonzero(grid.material==1)/np.prod(grid.cells),.5,rtol=1e-3)


    def test_from_table(self):
        cells = np.random.randint(60,100,3)
        size = np.ones(3)+np.random.rand(3)
        coords = grid_filters.coordinates0_point(cells,size).reshape(-1,3,order='F')
        z=np.ones(cells.prod())
        z[cells[:2].prod()*int(cells[2]/2):]=0
        t = Table(np.column_stack((coords,z)),{'coords':3,'z':1})
        t = t.add('indicator',t.get('coords')[:,0])
        g = Grid.from_table(t,'coords',['indicator','z'])
        assert g.N_materials == g.cells[0]*2 and (g.material[:,:,-1]-g.material[:,:,0] == cells[0]).all()


    def test_from_table_recover(self,tmp_path):
        cells = np.random.randint(60,100,3)
        size = np.ones(3)+np.random.rand(3)
        s = seeds.from_random(size,np.random.randint(60,100))
        grid = Grid.from_Voronoi_tessellation(cells,size,s)
        coords = grid_filters.coordinates0_point(cells,size)
        t = Table(np.column_stack((coords.reshape(-1,3,order='F'),grid.material.flatten(order='F'))),{'c':3,'m':1})
        assert grid_equal(grid.sort().renumber(),Grid.from_table(t,'c',['m']))


    @pytest.mark.parametrize('periodic',[True,False])
    @pytest.mark.parametrize('direction',['x','y','z',['x','y'],'zy','xz',['x','y','z']])
    def test_get_grain_boundaries(self,update,ref_path,periodic,direction):
        grid = Grid.load(ref_path/'get_grain_boundaries_8g12x15x20.vtr')
        current = grid.get_grain_boundaries(periodic,direction)
        if update:
            current.save(ref_path/f'get_grain_boundaries_8g12x15x20_{direction}_{periodic}.vtu',parallel=False)
        reference = VTK.load(ref_path/f'get_grain_boundaries_8g12x15x20_{"".join(direction)}_{periodic}.vtu')
        assert current.__repr__() == reference.__repr__()

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_get_grain_boundaries_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.get_grain_boundaries(directions=directions)

    def test_load_DREAM3D(self,ref_path):
        grain = Grid.load_DREAM3D(ref_path/'2phase_irregularGrid.dream3d','FeatureIds')
        point = Grid.load_DREAM3D(ref_path/'2phase_irregularGrid.dream3d')

        assert np.allclose(grain.origin,point.origin) and \
               np.allclose(grain.size,point.size) and \
               (grain.sort().material == point.material+1).all()


    def test_load_DREAM3D_reference(self,ref_path,update):
        current   = Grid.load_DREAM3D(ref_path/'measured.dream3d')
        reference = Grid.load(ref_path/'measured')
        if update:
            current.save(ref_path/'measured.vtr')

        assert grid_equal(current,reference)
