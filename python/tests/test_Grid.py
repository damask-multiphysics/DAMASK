import sys

import pytest
import numpy as np
from vtkmodules.vtkCommonCore import vtkVersion

from damask import VTK
from damask import Grid
from damask import Table
from damask import Rotation
from damask import Colormap
from damask import util
from damask import seeds
from damask import grid_filters


@pytest.fixture
def default():
    """Simple geometry."""
    g = np.array([8,5,4])
    l = np.prod(g[:2])
    return Grid(np.concatenate((np.ones  (l,dtype=int),
                                np.arange(l,dtype=int)+2,
                                np.ones  (l,dtype=int)*2,
                                np.arange(l,dtype=int)+1)).reshape(g,order='F'),
                g*1e-6)

@pytest.fixture
def random():
    """Simple geometry."""
    size = (1+np.random.rand(3))*1e-5
    cells = np.random.randint(10,20,3)
    s = seeds.from_random(size,np.random.randint(5,25),cells)
    return Grid.from_Voronoi_tessellation(cells,size,s)

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Grid'


class TestGrid:

    @pytest.fixture(autouse=True)
    def _patch_execution_stamp(self, patch_execution_stamp):
        print('patched damask.util.execution_stamp')

    @pytest.fixture(autouse=True)
    def _patch_datetime_now(self, patch_datetime_now):
        print('patched datetime.datetime.now')


    @pytest.mark.parametrize('cmap',[Colormap.from_predefined('stress'),'viridis'])
    @pytest.mark.skipif(sys.platform == 'win32', reason='DISPLAY has no effect on windows')
    def test_show(sef,default,cmap,monkeypatch):
        monkeypatch.delenv('DISPLAY',raising=False)
        default.show(cmap)

    def test_equal(self,default):
        assert default == default
        assert not default == 42

    def test_repr(self,default):
        print(default)

    def test_read_write_vti(self,default,tmp_path):
        default.save(tmp_path/'default')
        new = Grid.load(tmp_path/'default.vti')
        assert new == default

    def test_invalid_no_material(self,tmp_path):
        v = VTK.from_image_data(np.random.randint(5,10,3)*2,np.random.random(3) + 1.0)
        v.save(tmp_path/'no_materialpoint.vti',parallel=False)
        with pytest.raises(KeyError):
            Grid.load(tmp_path/'no_materialpoint.vti')

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
        assert Grid.load_ASCII(tmp_path/'ASCII') == default

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
    def test_mirror(self,default,update,res_path,directions,reflect):
        modified = default.mirror(directions,reflect)
        tag = f'directions_{"-".join(directions)}+reflect_{reflect}'
        reference = res_path/f'mirror_{tag}.vti'
        if update: modified.save(reference)
        assert Grid.load(reference) == modified

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_mirror_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.mirror(directions)

    @pytest.mark.parametrize('reflect',[True,False])
    def test_mirror_order_invariant(self,default,reflect):
        direction = np.array(['x','y','z'])
        assert default.mirror(np.random.permutation(direction),reflect=reflect) \
            == default.mirror(np.random.permutation(direction),reflect=reflect)

    @pytest.mark.parametrize('directions',[
                                           ['x'],
                                           ['x','y','z'],
                                           ['z','x','y'],
                                           ['y','z'],
                                          ]
                            )
    def test_flip(self,default,update,res_path,directions):
        modified = default.flip(directions)
        tag = f'directions_{"-".join(directions)}'
        reference = res_path/f'flip_{tag}.vti'
        if update: modified.save(reference)
        assert Grid.load(reference) == modified

    def test_flip_order_invariant(self,default):
        direction = np.array(['x','y','z'])
        assert default.flip(np.random.permutation(direction)) \
            == default.flip(np.random.permutation(direction))

    def test_flip_mirrored_invariant(self,default):
        direction = np.random.permutation(['x','y','z'])
        assert default.mirror(direction,True) == default.mirror(direction,True).flip(direction)

    def test_flip_equal_halfspin(self,default):
        direction = ['x','y','z']
        i = np.random.choice(3)
        assert default.rotate(Rotation.from_axis_angle(np.hstack((np.identity(3)[i],180)),degrees=True)) \
            == default.flip(direction[:i]+direction[i+1:])

    @pytest.mark.parametrize('direction',[['x'],['x','y']])
    def test_flip_double(self,default,direction):
        assert default == default.flip(direction).flip(direction)

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_flip_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.flip(directions)

    @pytest.mark.parametrize('distance',[1.,np.sqrt(3)])
    @pytest.mark.parametrize('selection',[None,1,[1],[1,2,3]])
    @pytest.mark.parametrize('periodic',[True,False])
    def test_clean_reference(self,default,update,res_path,distance,selection,periodic):
        current = default.clean(distance,selection,periodic=periodic,rng_seed=0)
        reference = res_path/f'clean_{distance}_{util.srepr(selection,"+")}_{periodic}.vti'
        if update:
            current.save(reference)
        assert Grid.load(reference) == current

    @pytest.mark.parametrize('selection',[list(np.random.randint(1,20,6)),np.random.randint(1,20,6)])
    @pytest.mark.parametrize('invert',[True,False])
    def test_clean_invert(self,default,selection,invert):
        selection_inverse = np.setdiff1d(default.material,selection)
        assert default.clean(selection=selection,invert_selection=invert,rng_seed=0) == \
               default.clean(selection=selection_inverse,invert_selection=not invert,rng_seed=0)

    def test_clean_selection_empty(self,random):
        assert random.clean(selection=None,invert_selection=True,rng_seed=0) == random.clean(rng_seed=0) and \
               random.clean(selection=None,invert_selection=False,rng_seed=0) == random.clean(rng_seed=0)

    @pytest.mark.parametrize('cells',[
                                     (10,11,10),
                                     [10,13,10],
                                     np.array((10,10,10)),
                                     np.array((8, 10,12)),
                                     np.array((5, 4, 20)),
                                     np.array((10,20,2))
                                    ]
                            )
    def test_scale(self,default,update,res_path,cells):
        modified = default.scale(cells)
        tag = f'grid_{util.srepr(cells,"-")}'
        reference = res_path/f'scale_{tag}.vti'
        if update: modified.save(reference)
        assert Grid.load(reference) == modified

    def test_renumber(self,default):
        material = default.material.copy()
        for m in np.unique(material):
            material[material==m] = material.max() + np.random.randint(1,30)
        default.material -= 1
        modified = Grid(material,
                        default.size,
                        default.origin)
        assert not default == modified
        assert     default == modified.renumber()

    def test_assemble(self):
        cells = np.random.randint(8,16,3)
        N = cells.prod()
        g = Grid(np.arange(N).reshape(cells),np.ones(3))
        idx = np.random.randint(0,N,N).reshape(cells)
        assert (idx == g.assemble(idx).material).all

    def test_substitute(self,default):
        offset = np.random.randint(1,500)
        modified = Grid(default.material + offset,
                        default.size,
                        default.origin)
        assert not default == modified
        assert     default == modified.substitute(np.arange(default.material.max())+1+offset,
                                                  np.arange(default.material.max())+1)

    def test_substitute_integer_list(self,random):
        f = np.random.randint(30)
        t = np.random.randint(30)
        assert random.substitute(f,t) == random.substitute([f],[t])

    def test_substitute_invariant(self,default):
        f = np.unique(default.material.flatten())[:np.random.randint(1,default.material.max())]
        t = np.random.permutation(f)
        modified = default.substitute(f,t)
        assert np.array_equiv(t,f) or modified != default
        assert default == modified.substitute(t,f)

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
        assert default == modified

    @pytest.mark.parametrize('Eulers',[[32.0,68.0,21.0],
                                       [0.0,32.0,240.0]])
    def test_rotate(self,default,update,res_path,Eulers):
        modified = default.rotate(Rotation.from_Euler_angles(Eulers,degrees=True))
        tag = f'Eulers_{util.srepr(Eulers,"-")}'
        reference = res_path/f'rotate_{tag}.vti'
        if update: modified.save(reference)
        assert Grid.load(reference) == modified

    def test_canvas_extend(self,default):
        cells = default.cells
        cells_add = np.random.randint(0,30,(3))
        modified = default.canvas(cells + cells_add)
        assert np.all(modified.material[:cells[0],:cells[1],:cells[2]] == default.material)

    @pytest.mark.parametrize('sign',[+1,-1])
    @pytest.mark.parametrize('extra_offset',[0,-1])
    def test_canvas_move_out(self,sign,extra_offset):
        g = Grid(np.zeros(np.random.randint(3,30,(3)),int),np.ones(3))
        o = sign*np.ones(3)*g.cells.min() +extra_offset*sign
        if extra_offset == 0:
            assert np.all(g.canvas(offset=o).material == 1)
        else:
            assert np.all(np.unique(g.canvas(offset=o).material) == (0,1))

    def test_canvas_cells(self,default):
        g = Grid(np.zeros(np.random.randint(3,30,(3)),int),np.ones(3))
        cells = np.random.randint(1,30,(3))
        offset = np.random.randint(-30,30,(3))
        assert np.all(g.canvas(cells,offset).cells == cells)

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
        assert G_1 == G_2

    @pytest.mark.parametrize('exponent',[1,np.inf,np.random.random(3)*2.])
    def test_add_primitive_shape_symmetry(self,exponent):
        """Shapes defined in the center should always produce a grid with reflection symmetry along the coordinate axis."""
        o = np.random.random(3)-.5
        s = np.random.random(3)*5.
        grid = Grid(np.zeros(np.random.randint(8,32,3),'i'),s,o).add_primitive(np.random.random(3)*3.,o+s/2.,exponent)
        for axis in [0,1,2]:
            assert np.all(grid.material==np.flip(grid.material,axis=axis))

    @pytest.mark.parametrize('selection',[1,None])
    def test_vicinity_offset(self,selection):
        offset = np.random.randint(2,4)
        distance = np.random.randint(2,4)

        g = np.random.randint(28,40,(3))
        m = np.ones(g,'i')
        x = (g*np.random.permutation(np.array([.5,1,1]))).astype(int)
        m[slice(0,x[0]),slice(0,x[1]),slice(0,x[2])] = 2
        m2 = m.copy()
        for i in [0,1,2]:
            m2[(np.roll(m,+distance,i)-m)!=0] += offset
            m2[(np.roll(m,-distance,i)-m)!=0] += offset
        if selection == 1:
            m2[m==1] = 1

        grid = Grid(m,np.random.rand(3)).vicinity_offset(distance,offset,selection=selection)

        assert np.all(m2==grid.material)

    @pytest.mark.parametrize('selection',[list(np.random.randint(1,20,6)),np.random.randint(1,20,6)])
    @pytest.mark.parametrize('invert',[True,False])
    def test_vicinity_offset_invert(self,random,selection,invert):
        selection_inverse = np.setdiff1d(random.material,selection)
        assert random.vicinity_offset(selection=selection        ,invert_selection=not invert) == \
               random.vicinity_offset(selection=selection_inverse,invert_selection=    invert)

    def test_vicinity_offset_selection_empty(self,random):
        assert random.vicinity_offset(selection=None,invert_selection=False) == random.vicinity_offset() and \
               random.vicinity_offset(selection=None,invert_selection=True ) == random.vicinity_offset()

    @pytest.mark.parametrize('periodic',[True,False])
    def test_vicinity_offset_invariant(self,default,periodic):
        offset = default.vicinity_offset(selection=[default.material.max()+1,
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
        assert Laguerre == Voronoi

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
        z = np.ones(cells.prod())
        z[cells[:2].prod()*int(cells[2]/2):] = 0
        t = Table({'coords':3,'z':1},np.column_stack((coords,z)))
        t = t.set('indicator',t.get('coords')[:,0])
        g = Grid.from_table(t,'coords',['indicator','z'])
        assert g.N_materials == g.cells[0]*2 and (g.material[:,:,-1]-g.material[:,:,0] == cells[0]).all()

    def test_from_table_recover(self,tmp_path):
        cells = np.random.randint(60,100,3)
        size = np.ones(3)+np.random.rand(3)
        s = seeds.from_random(size,np.random.randint(60,100))
        grid = Grid.from_Voronoi_tessellation(cells,size,s)
        coords = grid_filters.coordinates0_point(cells,size)
        t = Table({'c':3,'m':1},np.column_stack((coords.reshape(-1,3,order='F'),grid.material.flatten(order='F'))))
        assert grid.sort().renumber() == Grid.from_table(t,'c',['m'])

    @pytest.mark.parametrize('periodic',[True,False])
    @pytest.mark.parametrize('direction',['x','y','z',['x','y'],'zy','xz',['x','y','z']])
    @pytest.mark.xfail(vtkVersion.GetVTKMajorVersion()<8, reason='missing METADATA')
    def test_get_grain_boundaries(self,update,res_path,periodic,direction):
        grid = Grid.load(res_path/'get_grain_boundaries_8g12x15x20.vti')
        current = grid.get_grain_boundaries(periodic,direction)
        if update:
            current.save(res_path/f'get_grain_boundaries_8g12x15x20_{direction}_{periodic}.vtu',parallel=False)
        reference = VTK.load(res_path/f'get_grain_boundaries_8g12x15x20_{"".join(direction)}_{periodic}.vtu')
        assert current.__repr__() == reference.__repr__()

    @pytest.mark.parametrize('directions',[(1,2,'y'),('a','b','x'),[1]])
    def test_get_grain_boundaries_invalid(self,default,directions):
        with pytest.raises(ValueError):
            default.get_grain_boundaries(directions=directions)

    def test_load_DREAM3D(self,res_path):
        grain = Grid.load_DREAM3D(res_path/'2phase_irregularGrid.dream3d','FeatureIds')
        point = Grid.load_DREAM3D(res_path/'2phase_irregularGrid.dream3d')

        assert np.allclose(grain.origin,point.origin) and \
               np.allclose(grain.size,point.size) and \
               (grain.sort().material == point.material+1).all()

    def test_load_DREAM3D_reference(self,res_path,update):
        current   = Grid.load_DREAM3D(res_path/'measured.dream3d')
        reference = Grid.load(res_path/'measured.vti')
        if update:
            current.save(res_path/'measured.vti')

        assert current == reference

    def test_load_Neper_reference(self,res_path,update):
        current   = Grid.load_Neper(res_path/'n10-id1_scaled.vtk').renumber()
        reference = Grid.load(res_path/'n10-id1_scaled.vti')
        if update:
            current.save(res_path/'n10-id1_scaled.vti')

        assert current == reference
