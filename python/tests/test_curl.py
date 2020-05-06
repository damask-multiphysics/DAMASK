import pytest
import numpy as np
from damask import grid_filters


class TestGridFilters:


    @pytest.mark.parametrize('field_def,curl_def',
                             [(['np.sin(np.pi*2*nodes[...,2]/size[2])', '0.0', '0.0',
                                '0.0',                                  '0.0', '0.0',
                                '0.0',                                  '0.0', '0.0' ],
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
                              (['4*np.sin(np.pi*2*nodes[...,2]/size[2])', '8*np.sin(np.pi*2*nodes[...,0]/size[0])', '16*np.sin(np.pi*2*nodes[...,1]/size[1])'],
                               ['16*np.pi*2/size[1]*np.cos(np.pi*2*nodes[...,1]/size[1])', \
                                '4*np.pi*2/size[2]*np.cos(np.pi*2*nodes[...,2]/size[2])', \
                                '8*np.pi*2/size[0]*np.cos(np.pi*2*nodes[...,0]/size[0])']
                              ),
                              (['0.0', 'np.cos(np.pi*2*nodes[...,0]/size[0])+5*np.cos(np.pi*2*nodes[...,2]/size[2])', '0.0'],
                               ['5*np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',\
                                '0.0',\
                                '-np.sin(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]']
                              )
                             ])
    def test_curl(self,field_def,curl_def):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))

        nodes = grid_filters.cell_coord0(grid,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),grid) for f in field_def],axis=-1)
        field = field.reshape(tuple(grid) + ((3,3) if len(field_def)==9 else (3,)))
        curl = np.stack([np.broadcast_to(eval(c,globals(),my_locals),grid) for c in curl_def], axis=-1)
        curl = curl.reshape(tuple(grid) + ((3,3) if len(curl_def)==9 else (3,)))

        assert np.allclose(curl,grid_filters.curl(size,field))
        
        