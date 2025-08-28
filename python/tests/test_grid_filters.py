import pytest
import numpy as np

from damask import grid_filters
from damask import mechanics
from damask import seeds
from damask import GeomGrid


def test_coordinates0_point(np_rng):
    size = np_rng.random(3)
    cells = np_rng.integers(8,32,(3))
    coord = grid_filters.coordinates0_point(cells,size)
    assert np.allclose(coord[0,0,0],size/cells*.5) and coord.shape == tuple(cells) + (3,)

def test_coordinates0_node(np_rng):
    size = np_rng.random(3)
    cells = np_rng.integers(8,32,(3))
    coord = grid_filters.coordinates0_node(cells,size)
    assert np.allclose(coord[-1,-1,-1],size) and coord.shape == tuple(cells+1) + (3,)

def test_coord0(np_rng):
    size = np_rng.random(3)
    cells = np_rng.integers(8,32,(3))
    c = grid_filters.coordinates0_point(cells+1,size+size/cells)
    n = grid_filters.coordinates0_node(cells,size) + size/cells*.5
    assert np.allclose(c,n)


@pytest.mark.parametrize('mode',['point','node'])
@pytest.mark.parametrize('atol',[0.0,1.0e-5])
def test_grid_DNA(np_rng,mode,atol):
    """Ensure that cellsSizeOrigin_coordinates0_xx is the inverse of coordinates0_xx."""
    cells  = np_rng.integers(8,32,(3))
    size   = np_rng.random(3)
    origin = np_rng.random(3)
    coord0 = eval(f'grid_filters.coordinates0_{mode}(cells,size,origin)')                         # noqa
    coord0 += (np_rng.random(coord0.shape)-0.5) * atol                                            # add noise
    _cells,_size,_origin = eval(f'grid_filters.cellsSizeOrigin_coordinates0_{mode}(coord0.reshape(-1,3,order="F"),atol={atol})')
    assert np.allclose(cells,_cells) and np.allclose(size,_size,atol=atol) and np.allclose(origin,_origin,atol=atol)


def test_interpolation_to_node(np_rng):
    size  = np_rng.random(3)
    cells = np_rng.integers(8,32,(3))
    F     = np_rng.random(tuple(cells)+(3,3))
    assert np.allclose(grid_filters.coordinates_node(size,F) [1:-1,1:-1,1:-1],
                       grid_filters.point_to_node(grid_filters.coordinates_point(size,F))[1:-1,1:-1,1:-1])

def test_interpolation_to_cell(np_rng):
    cells = np_rng.integers(1,30,(3))

    coordinates_node_x = np.linspace(0,np.pi*2,num=cells[0]+1)
    node_field_x = np.cos(coordinates_node_x)
    node_field   = np.broadcast_to(node_field_x.reshape(-1,1,1),cells+1)

    coordinates0_point_x = coordinates_node_x[:-1]+coordinates_node_x[1]*.5
    cell_field_x = np.interp(coordinates0_point_x,coordinates_node_x,node_field_x,period=np.pi*2.)
    cell_field   = np.broadcast_to(cell_field_x.reshape(-1,1,1),cells)

    assert np.allclose(cell_field,grid_filters.node_to_point(node_field))


@pytest.mark.parametrize('mode',['point','node'])
def test_coordinates0_origin(np_rng,mode):
    origin= np_rng.random(3)
    size  = np_rng.random(3)
    cells  = np_rng.integers(8,32,(3))
    shifted   = eval(f'grid_filters.coordinates0_{mode}(cells,size,origin)')
    unshifted = eval(f'grid_filters.coordinates0_{mode}(cells,size)')
    if   mode == 'cell':
        assert  np.allclose(shifted,unshifted+np.broadcast_to(origin,tuple(cells)  +(3,)))
    elif mode == 'node':
        assert  np.allclose(shifted,unshifted+np.broadcast_to(origin,tuple(cells+1)+(3,)))

@pytest.mark.parametrize('function',[grid_filters.displacement_avg_point,
                                     grid_filters.displacement_avg_node])
def test_displacement_avg_vanishes(np_rng,function):
    """Ensure that random fluctuations in F do not result in average displacement."""
    size = np_rng.random(3) + 1.0
    cells = np_rng.integers(8,32,(3))
    F  = np_rng.random(tuple(cells)+(3,3))
    F += np.eye(3) - np.average(F,axis=(0,1,2))
    assert np.allclose(function(size,F),0.0)

@pytest.mark.parametrize('function',[grid_filters.displacement_avg_point,
                                     grid_filters.displacement_avg_node])
def test_displacement_avg_vanishes_simple(np_rng,function):
    F = np.eye(3)
    size = np_rng.random(3) + 1.0
    F_c = F.copy()
    F_t = F.copy()

    F_c[0,0] = 0.8
    F_t[0,0] = 1.2

    F_no_avg = np.concatenate([np.broadcast_to(_,(10,20,20,3,3)) for _ in [F_t,F_c]])
    assert np.allclose(function(size,F_no_avg),0.0)

@pytest.mark.parametrize('function',[grid_filters.displacement_fluct_point,
                                     grid_filters.displacement_fluct_node])
def test_displacement_fluct_vanishes_avg(np_rng,function):
    """Ensure that constant F does not result in fluctuating displacement."""
    size = np_rng.random(3)
    cells = np_rng.integers(8,32,(3))
    F    = np.broadcast_to(np_rng.random((3,3)), tuple(cells)+(3,3))
    assert np.allclose(function(size,F),0.0)

displacement_fluct_test_data = [
(['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0',
  '0.0',                                  '0.0', '0.0',
  '0.0',                                  '0.0', '0.0'],
 ['-np.cos(np.pi*2*nodes[...,0]/size[0])/np.pi/2*size[0]', '0.0', '0.0']
),
(['np.cos(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0',
  '0.0',                                  '0.0', '0.0',
  '0.0',                                  '0.0', 'np.cos(np.pi*2*nodes[...,2]/size[2])'],
 ['np.sin(np.pi*2*nodes[...,0]/size[0])/np.pi/2*size[0]',
  '0.0',
  'np.sin(np.pi*2*nodes[...,2]/size[2])/np.pi/2*size[2]']
)]
@pytest.mark.parametrize('F_def,u_def',displacement_fluct_test_data)
def test_displacment_fluct_analytic(np_rng,F_def,u_def):
    size = np_rng.random(3)+1.0
    cells = np_rng.integers(8,32,(3))

    nodes = grid_filters.coordinates0_point(cells,size)
    my_locals = locals()                                                                            # needed for list comprehension

    F = np.stack([np.broadcast_to(eval(F,globals(),my_locals),cells) for F in F_def],axis=-1).reshape(tuple(cells) + (3,3))
    u = np.stack([np.broadcast_to(eval(u,globals(),my_locals),cells) for u in u_def],axis=-1).reshape(tuple(cells) + (3,))

    assert np.allclose(u,grid_filters.displacement_fluct_point(size,F))

def test_displacement_fluct_periodic(np_rng):
    """Ensure that fluctuations are periodic."""
    size  = np_rng.random(3)
    cells = np_rng.integers(8,32,(3))
    F     = np_rng.random(tuple(cells)+(3,3))
    assert np.allclose(grid_filters.displacement_fluct_node(size,F),
                       grid_filters.point_to_node(grid_filters.displacement_fluct_point(size,F)))


def test_coordinates(np_rng):
    cells = np.array([np_rng.integers(40,100)*2,2,2])
    size = (np_rng.random(3)+0.8)*cells
    F = np.broadcast_to(np.eye(3),tuple(cells)+(3,3)).copy()
    F[...,0,0] += np.expand_dims(0.1*np.sin(np.linspace(0,2*np.pi,cells[0],False))+
                                 np_rng.random(cells[0])*0.05,(-1,-2))

    c_n = grid_filters.coordinates_node(size,F)[:,0,0,0]
    l_0 = (size/cells)[0]
    l = c_n[1:] - c_n[:-1]
    epsilon_reconstructed = (l-l_0)/l_0
    epsilon_direct = mechanics.strain(F,'V',1)[:,0,0,0,0]
    assert np.corrcoef(epsilon_reconstructed,epsilon_direct)[0,1] > 0.99

@pytest.mark.parametrize('function',[grid_filters.cellsSizeOrigin_coordinates0_point,
                                     grid_filters.cellsSizeOrigin_coordinates0_node])
def test_invalid_coordinates(np_rng,function):
    invalid_coordinates = np_rng.random((np_rng.integers(12,52),3))
    with pytest.raises(ValueError):
        function(invalid_coordinates)

@pytest.mark.parametrize('function',[grid_filters.coordinates0_point,
                                     grid_filters.coordinates0_node])
@pytest.mark.parametrize('atol',[0.0,1.0e-5])
def test_valid_coordinates_check(np_rng,function,atol):
    valid_coordinates = function(np_rng.integers(4,10,(3)),np_rng.random(3))
    valid_coordinates += (np_rng.random(valid_coordinates.shape)-0.5) * atol                        # add noise
    assert grid_filters.coordinates0_valid(valid_coordinates.reshape(-1,3,order='F'),atol=atol)

def test_invalid_coordinates_check(np_rng):
    invalid_coordinates = np_rng.random((np_rng.integers(12,52),3))
    assert not grid_filters.coordinates0_valid(invalid_coordinates)

@pytest.mark.parametrize('function',[grid_filters.cellsSizeOrigin_coordinates0_node,
                                     grid_filters.cellsSizeOrigin_coordinates0_point])
def test_uneven_spaced_coordinates(np_rng,function):
    start = np_rng.random(3)
    end   = np_rng.random(3)*10. + start
    cells = np_rng.integers(8,32,(3))
    uneven = np.stack(np.meshgrid(np.logspace(start[0],end[0],cells[0]),
                                  np.logspace(start[1],end[1],cells[1]),
                                  np.logspace(start[2],end[2],cells[2]),indexing = 'ij'),
                        axis = -1).reshape((cells.prod(),3),order='F')
    with pytest.raises(ValueError):
        function(uneven)

@pytest.mark.parametrize('mode',[True,False])
@pytest.mark.parametrize('function',[grid_filters.cellsSizeOrigin_coordinates0_node,
                                     grid_filters.cellsSizeOrigin_coordinates0_point])
def test_unordered_coordinates(np_rng,function,mode):
    origin = np_rng.random(3)
    size   = np_rng.random(3)*10.+origin
    cells  = np_rng.integers(8,32,(3))
    unordered = grid_filters.coordinates0_node(cells,size,origin).reshape(-1,3)
    if mode:
        with pytest.raises(ValueError):
            function(unordered,mode)
    else:
        function(unordered,mode)


def test_regrid_identity(np_rng):
        size = np_rng.random(3)
        cells = np_rng.integers(8,32,(3))
        F = np.broadcast_to(np.eye(3), (*cells,3,3))
        assert (grid_filters.regrid(size,F,cells).flatten(order='F') == np.arange(cells.prod())).all()

@pytest.mark.parametrize('factor',[2,3,4])
def test_regrid_fractional_shear(np_rng,factor):
        size = np.ones(3)
        scaling = np.array([1,1,factor])
        res = np_rng.integers(4,12)
        cells = factor*np.array([res]*3)
        F = np.broadcast_to([
            [1.0,0.0,1/factor],
            [0.0,1.0,1/factor],
            [0.0,0.0,1.0],
        ], (*cells,3,3))
        idx,box = grid_filters.regrid(size,
                                      F,
                                      cells*scaling,
                                      max_coeff=5,
                                      max_candidates=1000,
                                      return_size=True)
        counts = np.unique(idx, return_counts=True)[1]
        assert (counts == factor).all(), f'number of mismatches {np.count_nonzero(counts!=factor)}'
        assert np.allclose(box,size*scaling)

@pytest.mark.parametrize('factor',[1,2,3])
@pytest.mark.parametrize('noise',[1.1,1.0,0.9])
def test_regrid_stretch_with_shear(np_rng,factor,noise):
        size = np_rng.random(3)
        cells = np_rng.integers(8,32,(3))
        stretch = np_rng.random()
        F = np.broadcast_to([
            [noise*stretch,0,factor*stretch*size[0]/size[2]],
            [0,noise*stretch,factor*stretch*size[1]/size[2]],
            [0,0,1],
            ], (*cells,3,3))
        if noise == 1.0:
            assert set(grid_filters.regrid(size,F,cells).flatten()).issubset(np.arange(np.prod(cells)))
        else:
            with pytest.raises(ValueError):
                grid_filters.regrid(size,F,cells)

def test_regrid_double_cells(np_rng):
        size = np_rng.random(3)
        cells = np_rng.integers(8,32,(3))
        g = GeomGrid.from_Voronoi_tessellation(cells,size,seeds.from_random(size,10,rng_seed=np_rng))
        F = np.broadcast_to(np.eye(3), (*cells,3,3))
        assert g.scale(cells*2) == g.assemble(grid_filters.regrid(size,F,cells*2))


@pytest.mark.parametrize('differential_operator',[grid_filters.curl,
                                                  grid_filters.divergence,
                                                  grid_filters.gradient])
def test_differential_operator_constant(np_rng,differential_operator):
    size = np_rng.random(3)+1.0
    cells = np_rng.integers(8,32,(3))
    shapes = {
              grid_filters.curl:      [(3,),(3,3)],
              grid_filters.divergence:[(3,),(3,3)],
              grid_filters.gradient:  [(1,),(3,)]
             }
    for shape in shapes[differential_operator]:
        field = np.ones(tuple(cells)+shape)*np_rng.random()*1.0e5
        assert np.allclose(differential_operator(size,field),0.0)


grad_test_data = [
(['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
 ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0',
  '0.0',                                                  '0.0', '0.0',
  '0.0',                                                  '0.0', '0.0']
),
(['0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0' ],
 ['0.0', '0.0',                                                   '0.0',
  '0.0', '-np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])', '0.0',
  '0.0', '0.0',                                                   '0.0']
),
(['1.0', '0.0', '2.0*np.cos(np.pi*2*nodes[...,2]/size[2])'],
 ['0.0', '0.0', '0.0',
  '0.0', '0.0', '0.0',
  '0.0', '0.0', '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']
),
(['np.cos(np.pi*2*nodes[...,2]/size[2])', '3.0', 'np.sin(np.pi*2*nodes[...,2]/size[2])'],
 ['0.0', '0.0', '-np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
  '0.0', '0.0', '0.0',
  '0.0', '0.0', ' np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']
),
(['np.sin(np.pi*2*nodes[...,0]/size[0])',
  'np.sin(np.pi*2*nodes[...,1]/size[1])',
  'np.sin(np.pi*2*nodes[...,2]/size[2])'],
 ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0',
  '0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0',
  '0.0', '0.0', 'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']
),
(['np.sin(np.pi*2*nodes[...,0]/size[0])'],
 ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0']
),
(['8.0'],
 ['0.0', '0.0', '0.0' ]
)]
@pytest.mark.parametrize('field_def,grad_def',grad_test_data)
def test_grad(np_rng,field_def,grad_def):
    size = np_rng.random(3)+1.0
    cells = np_rng.integers(8,32,(3))

    nodes = grid_filters.coordinates0_point(cells,size)
    my_locals = locals()                                                                            # needed for list comprehension

    field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),cells) for f in field_def],axis=-1)
    field = field.reshape(tuple(cells) + ((3,) if len(field_def)==3 else (1,)))
    grad = np.stack([np.broadcast_to(eval(c,globals(),my_locals),cells) for c in grad_def], axis=-1)
    grad = grad.reshape(tuple(cells) + ((3,3) if len(grad_def)==9 else (3,)))

    assert np.allclose(grad,grid_filters.gradient(size,field))


curl_test_data = [
(['np.sin(np.pi*2*nodes[...,2]/size[2])', '0.0', '0.0',
  '0.0',                                  '0.0', '0.0',
  '0.0',                                  '0.0', '0.0'],
 ['0.0'                                                 , '0.0', '0.0',
  'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]', '0.0', '0.0',
  '0.0',                                                  '0.0', '0.0']
),
(['np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0', '0.0',
  '0.0',                                  '0.0', '0.0',
  'np.cos(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
 ['0.0',                                                  '0.0', '0.0',
  '0.0',                                                  '0.0', '0.0',
  'np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0', '0.0']
),
(['np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])',
  'np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])',
  'np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])'],
 ['0.0', '0.0', '0.0',
  '0.0', '0.0', '0.0',
  '0.0', '0.0', '0.0']
),
(['5.0', '0.0', '0.0',
  '0.0', '0.0', '0.0',
  '0.0', '0.0', '2*np.cos(np.pi*2*nodes[...,1]/size[1])'],
 ['0.0', '0.0', '-2*np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])',
  '0.0', '0.0', '0.0',
  '0.0', '0.0', '0.0']
),
([ '4*np.sin(np.pi*2*nodes[...,2]/size[2])',
   '8*np.sin(np.pi*2*nodes[...,0]/size[0])',
   '16*np.sin(np.pi*2*nodes[...,1]/size[1])'],
 ['16*np.pi*2/size[1]*np.cos(np.pi*2*nodes[...,1]/size[1])',
  '4*np.pi*2/size[2]*np.cos(np.pi*2*nodes[...,2]/size[2])',
  '8*np.pi*2/size[0]*np.cos(np.pi*2*nodes[...,0]/size[0])']
),
(['0.0',
  'np.cos(np.pi*2*nodes[...,0]/size[0])+5*np.cos(np.pi*2*nodes[...,2]/size[2])',
  '0.0'],
 ['5*np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
  '0.0',
  '-np.sin(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]']
)]
@pytest.mark.parametrize('field_def,curl_def',curl_test_data)
def test_curl(np_rng,field_def,curl_def):
    size = np_rng.random(3)+1.0
    cells = np_rng.integers(8,32,(3))

    nodes = grid_filters.coordinates0_point(cells,size)
    my_locals = locals()                                                                            # needed for list comprehension

    field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),cells) for f in field_def],axis=-1)
    field = field.reshape(tuple(cells) + ((3,3) if len(field_def)==9 else (3,)))
    curl = np.stack([np.broadcast_to(eval(c,globals(),my_locals),cells) for c in curl_def], axis=-1)
    curl = curl.reshape(tuple(cells) + ((3,3) if len(curl_def)==9 else (3,)))

    assert np.allclose(curl,grid_filters.curl(size,field))


div_test_data =[
(['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0',
  '0.0'                                 , '0.0', '0.0',
  '0.0'                                 , '0.0', '0.0'],
 ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]','0.0', '0.0']
),
(['0.0', '0.0',                                  '0.0',
  '0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0',
  '0.0', '0.0',                                  '0.0'],
 ['0.0', '-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0']
),
(['1.0', '0.0', '0.0',
  '0.0', '0.0', '0.0',
  '0.0', '0.0', '2*np.cos(np.pi*2*nodes[...,2]/size[2])' ],
 ['0.0', '0.0', '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']
),
([ '23.0', '0.0',   'np.sin(np.pi*2*nodes[...,2]/size[2])',
   '0.0',  '100.0', 'np.sin(np.pi*2*nodes[...,2]/size[2])',
   '0.0',  '0.0',   'np.sin(np.pi*2*nodes[...,2]/size[2])'],
 ['np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
  'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
  'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']
 ),
(['400.0',                                '0.0',                                  '0.0',
  'np.sin(np.pi*2*nodes[...,0]/size[0])', 'np.sin(np.pi*2*nodes[...,1]/size[1])', 'np.sin(np.pi*2*nodes[...,2]/size[2])',
  '0.0',                                  '10.0',                                 '6.0'],
 ['0.0','np.sum(np.cos(np.pi*2*nodes/size)*np.pi*2/size,axis=-1)', '0.0' ]
),
(['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
 ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]',]
),
(['0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0' ],
 ['-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]']
)]
@pytest.mark.parametrize('field_def,div_def',div_test_data)
def test_div(np_rng,field_def,div_def):
    size = np_rng.random(3)+1.0
    cells = np_rng.integers(8,32,(3))

    nodes = grid_filters.coordinates0_point(cells,size)
    my_locals = locals()                                                                            # needed for list comprehension

    field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),cells) for f in field_def],axis=-1)
    field = field.reshape(tuple(cells) + ((3,3) if len(field_def)==9 else (3,)))
    div = np.stack([np.broadcast_to(eval(c,globals(),my_locals),cells) for c in div_def], axis=-1)
    div = div.reshape(tuple(cells) + ((3,) if len(div_def)==3 else ()))

    assert np.allclose(div,grid_filters.divergence(size,field))


def test_ravel_index(np_rng):
    cells = np_rng.integers(8,32,(3))

    indices = np.block(list(np.meshgrid(np.arange(cells[0]),
                                        np.arange(cells[1]),
                                        np.arange(cells[2]),indexing='ij'))).reshape(tuple(cells)+(3,),order='F')
    x,y,z = map(np_rng.integers,cells)
    assert grid_filters.ravel_index(indices)[x,y,z] == np.arange(0,np.prod(cells)).reshape(cells,order='F')[x,y,z]

def test_unravel_index(np_rng):
    cells = np_rng.integers(8,32,(3))
    indices = np.arange(np.prod(cells)).reshape(cells,order='F')
    x,y,z = map(np_rng.integers,cells)
    assert np.all(grid_filters.unravel_index(indices)[x,y,z] == [x,y,z])

def test_unravel_ravel_index(np_rng):
    cells = np_rng.integers(8,32,(3))
    indices = np_rng.integers(0,np.prod(cells),cells).reshape(cells)
    assert np.all(indices==grid_filters.ravel_index(grid_filters.unravel_index(indices)))

def test_ravel_unravel_index(np_rng):
    cells = np.hstack([np_rng.integers(8,32,(3)),1])
    indices = np.block([np_rng.integers(0,cells[0],cells),
                        np_rng.integers(0,cells[1],cells),
                        np_rng.integers(0,cells[2],cells)])
    assert np.all(indices==grid_filters.unravel_index(grid_filters.ravel_index(indices)))


def test_unravel_ravel(np_rng):
    shape = tuple(np_rng.integers(1,32,(np_rng.integers(3,6))))
    f = np_rng.random(shape).reshape((np.prod(shape[:3]),)+shape[3:])
    assert np.all(f == grid_filters.ravel(grid_filters.unravel(f,shape[:3])))

def test_ravel_unravel(np_rng):
    shape = np_rng.integers(1,32,(np_rng.integers(3,6)))
    f = np_rng.random(shape)
    assert np.all(f == grid_filters.unravel(grid_filters.ravel(f),shape[:3]))

def test_ravel_unravel_consistency(np_rng):
    cells = tuple(np_rng.integers(8,32,(3)))
    data_3D = np_rng.random(cells+tuple([np_rng.integers(2,4)] * np_rng.integers(3)))
    data_1D = grid_filters.ravel(data_3D)
    indices_1D = np_rng.integers(0,np.prod(cells),cells).reshape(cells)
    indices_3D = grid_filters.unravel_index(indices_1D)
    random_cell = (np_rng.integers(cells[0]),np_rng.integers(cells[1]),np_rng.integers(cells[2]))
    assert np.all(data_1D[indices_1D[random_cell]]==data_3D[tuple(indices_3D[random_cell])])
