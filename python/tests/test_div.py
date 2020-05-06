import pytest
import numpy as np
from damask import grid_filters


class TestGridFilters:


    div_test_data =[(['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0',
                      '0.0'                                 , '0.0', '0.0',
                      '0.0'                                 , '0.0', '0.0'],
                     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]' ,'0.0', '0.0']
                     ),
                    (['0.0',                                   '0.0', '0.0',
                      '0.0',  'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0',
                      '0.0',                                   '0.0', '0.0'],
                     ['0.0', '-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0']
                     ),
                    (['1.0',                          '0.0',                            '0.0',
                      '0.0',                          '0.0',                            '0.0',           
                      '0.0',                          '0.0',                            '2*np.cos(np.pi*2*nodes[...,2]/size[2])'            ],
                     ['0.0',                          '0.0',                            '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']
                     ),
                    ([ '23.0',                         '0.0',                            'np.sin(np.pi*2*nodes[...,2]/size[2])',    
                       '0.0',                          '100.0',                          'np.sin(np.pi*2*nodes[...,2]/size[2])',              
                       '0.0',                          '0.0',                            'np.sin(np.pi*2*nodes[...,2]/size[2])'],
                      ['np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]','np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',    'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']
                     ),
                    (['400.0',                        '0.0',                            '0.0',
                      'np.sin(np.pi*2*nodes[...,0]/size[0])',            'np.sin(np.pi*2*nodes[...,1]/size[1])',              'np.sin(np.pi*2*nodes[...,2]/size[2])', 
                      '0.0',          '10.0',                           '6.0'                              ],
                     ['0.0','np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]+np.cos(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]+np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]','0.0' ]
                     ),
                    (['np.sin(np.pi*2*nodes[...,0]/size[0])',            '0.0',                            '0.0'],
                     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]',]
                     ),
                    (['0.0',                          'np.cos(np.pi*2*nodes[...,1]/size[1])',              '0.0' ],
                     ['-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]']
                     )
                    
                    
                     ]


    @pytest.mark.parametrize('field_def,div_def',div_test_data)                        
                             
    def test_div(self,field_def,div_def):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))

        nodes = grid_filters.cell_coord0(grid,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),grid) for f in field_def],axis=-1)
        field = field.reshape(tuple(grid) + ((3,3) if len(field_def)==9 else (3,)))
        div = np.stack([np.broadcast_to(eval(c,globals(),my_locals),grid) for c in div_def], axis=-1)
        if len(div_def)==3:
            div = div.reshape(tuple(grid) + ((3,)))
        else:
            div=div.reshape(tuple(grid))

        assert np.allclose(div,grid_filters.divergence(size,field))
        