import pytest
import numpy as np

from damask import grid_filters

class TestGridFilters:
   
    def test_cell_coord0(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         coord = grid_filters.cell_coord0(grid,size)
         assert np.allclose(coord[0,0,0],size/grid*.5) and coord.shape == tuple(grid[::-1]) + (3,)

    def test_node_coord0(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         coord = grid_filters.node_coord0(grid,size)
         assert np.allclose(coord[-1,-1,-1],size) and coord.shape == tuple(grid[::-1]+1) + (3,)

    def test_coord0(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         c = grid_filters.cell_coord0(grid+1,size+size/grid)
         n = grid_filters.node_coord0(grid,size) + size/grid*.5
         assert np.allclose(c,n)

    @pytest.mark.parametrize('mode',[('cell'),('node')])
    def test_displacement_avg_vanishes(self,mode):
         """Ensure that random fluctuations in F do not result in average displacement."""
         size = np.random.random(3)                                                                 # noqa
         grid = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(grid)+(3,3))
         F   += np.eye(3) - np.average(F,axis=(0,1,2))
         assert np.allclose(eval('grid_filters.{}_displacement_avg(size,F)'.format(mode)),0.0)

    @pytest.mark.parametrize('mode',[('cell'),('node')])
    def test_displacement_fluct_vanishes(self,mode):
         """Ensure that constant F does not result in fluctuating displacement."""
         size = np.random.random(3)                                                                 # noqa
         grid = np.random.randint(8,32,(3))
         F    = np.broadcast_to(np.random.random((3,3)), tuple(grid)+(3,3))                         # noqa
         assert np.allclose(eval('grid_filters.{}_displacement_fluct(size,F)'.format(mode)),0.0)
