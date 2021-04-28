import pytest
import numpy as np

from damask import grid_filters

class TestGridFilters:

    def test_coordinates0_point(self):
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         coord = grid_filters.coordinates0_point(cells,size)
         assert np.allclose(coord[0,0,0],size/cells*.5) and coord.shape == tuple(cells) + (3,)

    def test_coordinates0_node(self):
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         coord = grid_filters.coordinates0_node(cells,size)
         assert np.allclose(coord[-1,-1,-1],size) and coord.shape == tuple(cells+1) + (3,)

    def test_coord0(self):
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         c = grid_filters.coordinates0_point(cells+1,size+size/cells)
         n = grid_filters.coordinates0_node(cells,size) + size/cells*.5
         assert np.allclose(c,n)

    @pytest.mark.parametrize('mode',['point','node'])
    def test_grid_DNA(self,mode):
         """Ensure that cellsSizeOrigin_coordinates0_xx is the inverse of coordinates0_xx."""
         cells   = np.random.randint(8,32,(3))
         size   = np.random.random(3)
         origin = np.random.random(3)
         coord0 = eval(f'grid_filters.coordinates0_{mode}(cells,size,origin)')                     # noqa
         _cells,_size,_origin = eval(f'grid_filters.cellsSizeOrigin_coordinates0_{mode}(coord0.reshape(-1,3,order="F"))')
         assert np.allclose(cells,_cells) and np.allclose(size,_size) and np.allclose(origin,_origin)

    def test_displacement_fluct_equivalence(self):
         """Ensure that fluctuations are periodic."""
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(cells)+(3,3))
         assert np.allclose(grid_filters.displacement_fluct_node(size,F),
                            grid_filters.point_to_node(grid_filters.displacement_fluct_point(size,F)))

    def test_interpolation_to_node(self):
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(cells)+(3,3))
         assert np.allclose(grid_filters.coordinates_node(size,F) [1:-1,1:-1,1:-1],
                            grid_filters.point_to_node(grid_filters.coordinates_point(size,F))[1:-1,1:-1,1:-1])

    def test_interpolation_to_cell(self):
         cells = np.random.randint(1,30,(3))

         coordinates_node_x = np.linspace(0,np.pi*2,num=cells[0]+1)
         node_field_x = np.cos(coordinates_node_x)
         node_field   = np.broadcast_to(node_field_x.reshape(-1,1,1),cells+1)

         coordinates0_point_x = coordinates_node_x[:-1]+coordinates_node_x[1]*.5
         cell_field_x = np.interp(coordinates0_point_x,coordinates_node_x,node_field_x,period=np.pi*2.)
         cell_field   = np.broadcast_to(cell_field_x.reshape(-1,1,1),cells)

         assert np.allclose(cell_field,grid_filters.node_to_point(node_field))

    @pytest.mark.parametrize('mode',['point','node'])
    def test_coordinates0_origin(self,mode):
         origin= np.random.random(3)
         size  = np.random.random(3)                                                                # noqa
         cells  = np.random.randint(8,32,(3))
         shifted   = eval(f'grid_filters.coordinates0_{mode}(cells,size,origin)')
         unshifted = eval(f'grid_filters.coordinates0_{mode}(cells,size)')
         if   mode == 'cell':
            assert  np.allclose(shifted,unshifted+np.broadcast_to(origin,tuple(cells)  +(3,)))
         elif mode == 'node':
            assert  np.allclose(shifted,unshifted+np.broadcast_to(origin,tuple(cells+1)+(3,)))

    @pytest.mark.parametrize('function',[grid_filters.displacement_avg_point,
                                         grid_filters.displacement_avg_node])
    def test_displacement_avg_vanishes(self,function):
         """Ensure that random fluctuations in F do not result in average displacement."""
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(cells)+(3,3))
         F   += np.eye(3) - np.average(F,axis=(0,1,2))
         assert np.allclose(function(size,F),0.0)

    @pytest.mark.parametrize('function',[grid_filters.displacement_fluct_point,
                                         grid_filters.displacement_fluct_node])
    def test_displacement_fluct_vanishes(self,function):
         """Ensure that constant F does not result in fluctuating displacement."""
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         F    = np.broadcast_to(np.random.random((3,3)), tuple(cells)+(3,3))
         assert np.allclose(function(size,F),0.0)

    @pytest.mark.parametrize('function',[grid_filters.cellsSizeOrigin_coordinates0_point,
                                         grid_filters.cellsSizeOrigin_coordinates0_node])
    def test_invalid_coordinates(self,function):
        invalid_coordinates = np.random.random((np.random.randint(12,52),3))
        with pytest.raises(ValueError):
            function(invalid_coordinates)

    @pytest.mark.parametrize('function',[grid_filters.coordinates0_point,
                                         grid_filters.coordinates0_node])
    def test_valid_coordinates_check(self,function):
        valid_coordinates = function(np.random.randint(4,10,(3)),np.random.rand(3))
        assert grid_filters.coordinates0_valid(valid_coordinates.reshape(-1,3,order='F'))

    def test_invalid_coordinates_check(self):
        invalid_coordinates = np.random.random((np.random.randint(12,52),3))
        assert not grid_filters.coordinates0_valid(invalid_coordinates)


    @pytest.mark.parametrize('function',[grid_filters.cellsSizeOrigin_coordinates0_node,
                                         grid_filters.cellsSizeOrigin_coordinates0_point])
    def test_uneven_spaced_coordinates(self,function):
        start = np.random.random(3)
        end   = np.random.random(3)*10. + start
        cells  = np.random.randint(8,32,(3))
        uneven = np.stack(np.meshgrid(np.logspace(start[0],end[0],cells[0]),
                                      np.logspace(start[1],end[1],cells[1]),
                                      np.logspace(start[2],end[2],cells[2]),indexing = 'ij'),
                           axis = -1).reshape((cells.prod(),3),order='F')
        with pytest.raises(ValueError):
            function(uneven)


    @pytest.mark.parametrize('mode',[True,False])
    @pytest.mark.parametrize('function',[grid_filters.cellsSizeOrigin_coordinates0_node,
                                         grid_filters.cellsSizeOrigin_coordinates0_point])
    def test_unordered_coordinates(self,function,mode):
        origin = np.random.random(3)
        size   = np.random.random(3)*10.+origin
        cells  = np.random.randint(8,32,(3))
        unordered = grid_filters.coordinates0_node(cells,size,origin).reshape(-1,3)
        if mode:
            with pytest.raises(ValueError):
                function(unordered,mode)
        else:
            function(unordered,mode)

    def test_regrid(self):
         size = np.random.random(3)
         cells = np.random.randint(8,32,(3))
         F    = np.broadcast_to(np.eye(3), tuple(cells)+(3,3))
         assert all(grid_filters.regrid(size,F,cells) == np.arange(cells.prod()))


    @pytest.mark.parametrize('differential_operator',[grid_filters.curl,
                                                      grid_filters.divergence,
                                                      grid_filters.gradient])
    def test_differential_operator_constant(self,differential_operator):
        size = np.random.random(3)+1.0
        cells = np.random.randint(8,32,(3))
        shapes = {
                  grid_filters.curl:      [(3,),(3,3)],
                  grid_filters.divergence:[(3,),(3,3)],
                  grid_filters.gradient:  [(1,),(3,)]
                 }
        for shape in shapes[differential_operator]:
            field = np.ones(tuple(cells)+shape)*np.random.random()*1.0e5
            assert np.allclose(differential_operator(size,field),0.0)


    grad_test_data = [
    (['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0']),

    (['0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0' ],
     ['0.0', '0.0',                                                   '0.0',
      '0.0', '-np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])', '0.0',
      '0.0', '0.0',                                                   '0.0' ]),

    (['1.0', '0.0', '2.0*np.cos(np.pi*2*nodes[...,2]/size[2])'],
     ['0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']),

    (['np.cos(np.pi*2*nodes[...,2]/size[2])', '3.0', 'np.sin(np.pi*2*nodes[...,2]/size[2])'],
     ['0.0', '0.0', '-np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', ' np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])',
      'np.sin(np.pi*2*nodes[...,1]/size[1])',
      'np.sin(np.pi*2*nodes[...,2]/size[2])'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0',
      '0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0',
      '0.0', '0.0', 'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0']),

    (['8.0'],
     ['0.0', '0.0', '0.0' ])
                    ]

    @pytest.mark.parametrize('field_def,grad_def',grad_test_data)
    def test_grad(self,field_def,grad_def):
        size = np.random.random(3)+1.0
        cells = np.random.randint(8,32,(3))

        nodes = grid_filters.coordinates0_point(cells,size)
        my_locals = locals()                                                                        # needed for list comprehension

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
      '0.0',                                                  '0.0', '0.0']),

    (['np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0', '0.0',
      '0.0',                                  '0.0', '0.0',
      'np.cos(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
     ['0.0',                                                  '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0',
      'np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0', '0.0']),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])',
      'np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])',
      'np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])'],
     ['0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0']),

    (['5.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '2*np.cos(np.pi*2*nodes[...,1]/size[1])'],
     ['0.0', '0.0', '-2*np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0']),

    ([ '4*np.sin(np.pi*2*nodes[...,2]/size[2])',
       '8*np.sin(np.pi*2*nodes[...,0]/size[0])',
      '16*np.sin(np.pi*2*nodes[...,1]/size[1])'],
     ['16*np.pi*2/size[1]*np.cos(np.pi*2*nodes[...,1]/size[1])',
       '4*np.pi*2/size[2]*np.cos(np.pi*2*nodes[...,2]/size[2])',
       '8*np.pi*2/size[0]*np.cos(np.pi*2*nodes[...,0]/size[0])']),

    (['0.0',
      'np.cos(np.pi*2*nodes[...,0]/size[0])+5*np.cos(np.pi*2*nodes[...,2]/size[2])',
      '0.0'],
     ['5*np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
      '0.0',
      '-np.sin(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]'])
                     ]

    @pytest.mark.parametrize('field_def,curl_def',curl_test_data)
    def test_curl(self,field_def,curl_def):
        size = np.random.random(3)+1.0
        cells = np.random.randint(8,32,(3))

        nodes = grid_filters.coordinates0_point(cells,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),cells) for f in field_def],axis=-1)
        field = field.reshape(tuple(cells) + ((3,3) if len(field_def)==9 else (3,)))
        curl = np.stack([np.broadcast_to(eval(c,globals(),my_locals),cells) for c in curl_def], axis=-1)
        curl = curl.reshape(tuple(cells) + ((3,3) if len(curl_def)==9 else (3,)))

        assert np.allclose(curl,grid_filters.curl(size,field))


    div_test_data =[
    (['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0',
      '0.0'                                 , '0.0', '0.0',
      '0.0'                                 , '0.0', '0.0'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]','0.0', '0.0']),

    (['0.0', '0.0',                                  '0.0',
      '0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0',
      '0.0', '0.0',                                  '0.0'],
     ['0.0', '-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0']),

    (['1.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '2*np.cos(np.pi*2*nodes[...,2]/size[2])' ],
     ['0.0', '0.0', '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']
     ),

    ([ '23.0', '0.0',   'np.sin(np.pi*2*nodes[...,2]/size[2])',
       '0.0',  '100.0', 'np.sin(np.pi*2*nodes[...,2]/size[2])',
       '0.0',  '0.0',   'np.sin(np.pi*2*nodes[...,2]/size[2])'],
      ['np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',\
       'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]', \
       'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),

    (['400.0',                                '0.0',                                  '0.0',
      'np.sin(np.pi*2*nodes[...,0]/size[0])', 'np.sin(np.pi*2*nodes[...,1]/size[1])', 'np.sin(np.pi*2*nodes[...,2]/size[2])',
      '0.0',                                  '10.0',                                 '6.0'],
     ['0.0','np.sum(np.cos(np.pi*2*nodes/size)*np.pi*2/size,axis=-1)', '0.0' ]),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]',]),

    (['0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0' ],
     ['-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]'])
     ]

    @pytest.mark.parametrize('field_def,div_def',div_test_data)

    def test_div(self,field_def,div_def):
        size = np.random.random(3)+1.0
        cells = np.random.randint(8,32,(3))

        nodes = grid_filters.coordinates0_point(cells,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),cells) for f in field_def],axis=-1)
        field = field.reshape(tuple(cells) + ((3,3) if len(field_def)==9 else (3,)))
        div = np.stack([np.broadcast_to(eval(c,globals(),my_locals),cells) for c in div_def], axis=-1)
        if len(div_def)==3:
            div = div.reshape(tuple(cells) + ((3,)))
        else:
            div=div.reshape(tuple(cells))

        assert np.allclose(div,grid_filters.divergence(size,field))
