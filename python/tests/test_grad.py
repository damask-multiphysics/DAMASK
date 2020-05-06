import pytest
import numpy as np
from damask import grid_filters


class TestGridFilters:

    grad_test_data = [(['np.sin(np.pi*2*nodes[...,0]/size[0])',                '0.0',                            '0.0'],
                       ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]','0.0',                            '0.0',
                        '0.0',                                                 '0.0',                            '0.0',
                        '0.0',                                                 '0.0',                            '0.0']),              
                      (['0.0',           'np.cos(np.pi*2*nodes[...,1]/size[1])',                   '0.0' ],
                       ['0.0',           '0.0',                                                    '0.0',
                        '0.0',           '-np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])',  '0.0',
                        '0.0',           '0.0',                                                    '0.0' ]),
                      (['1.0',           '0.0',                '2.0*np.cos(np.pi*2*nodes[...,2]/size[2])'],
                       ['0.0',           '0.0',                '0.0',
                        '0.0',           '0.0',                '0.0',
                        '0.0',           '0.0',                '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']),                
                      (['np.cos(np.pi*2*nodes[...,2]/size[2])', '3.0',  'np.sin(np.pi*2*nodes[...,2]/size[2])'],
                       ['0.0',                                  '0.0',  '-np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
                        '0.0',                                  '0.0',  '0.0',
                        '0.0',                                  '0.0',  'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),
                      (['np.sin(np.pi*2*nodes[...,0]/size[0])','np.sin(np.pi*2*nodes[...,1]/size[1])',\
                                                              'np.sin(np.pi*2*nodes[...,2]/size[2])'],
                       ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0',                            '0.0',  
                        '0.0',          'np.cos(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]',   '0.0',
                        '0.0',          '0.0',                          'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),
                      (['np.sin(np.pi*2*nodes[...,0]/size[0])'                                                               ],
                       ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0',                            '0.0'      ]),
                      (['8.0'                                                                                                ],
                       ['0.0',                          '0.0',                            '0.0'                              ])]

    @pytest.mark.parametrize('field_def,grad_def',grad_test_data)

    def test_grad(self,field_def,grad_def):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))
 
        nodes = grid_filters.cell_coord0(grid,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),grid) for f in field_def],axis=-1)
        field = field.reshape(tuple(grid) + ((3,) if len(field_def)==3 else (1,)))
        grad = np.stack([np.broadcast_to(eval(c,globals(),my_locals),grid) for c in grad_def], axis=-1)
        grad = grad.reshape(tuple(grid) + ((3,3) if len(grad_def)==9 else (3,)))

        assert np.allclose(grad,grid_filters.gradient(size,field))
        