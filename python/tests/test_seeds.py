import pytest
import numpy as np
from scipy.spatial import cKDTree

from damask import seeds
from damask import grid_filters
from damask import Grid

class TestSeeds:

    @pytest.mark.parametrize('cells',[None,np.ones(3,dtype='i')*10])
    def test_from_random(self,cells):
        N_seeds = np.random.randint(30,300)
        size = np.ones(3) + np.random.random(3)
        coords = seeds.from_random(size,N_seeds,cells)
        assert (0<=coords).all() and (coords<size).all()

    @pytest.mark.parametrize('periodic',[True,False])
    def test_from_Poisson_disc(self,periodic):
        N_seeds = np.random.randint(30,300)
        N_candidates = N_seeds//15
        distance = np.random.random()
        size = np.ones(3)*distance*N_seeds
        coords = seeds.from_Poisson_disc(size,N_seeds,N_candidates,distance,periodic=periodic)
        min_dists, _ = cKDTree(coords,boxsize=size).query(coords, 2) if periodic else \
                       cKDTree(coords).query(coords, 2)
        assert (0<= coords).all() and (coords<size).all() and np.min(min_dists[:,1])>=distance

    @pytest.mark.parametrize('periodic',[True,False])
    def test_from_Poisson_disc_invalid(self,periodic):
        N_seeds = np.random.randint(30,300)
        N_candidates = N_seeds//15
        distance = np.random.random()
        size = np.ones(3)*distance
        with pytest.raises(ValueError):
            seeds.from_Poisson_disc(size,N_seeds,N_candidates,distance,periodic=periodic)

    def test_from_grid_reconstruct(self):
        cells = np.random.randint(10,20,3)
        N_seeds = np.random.randint(30,300)
        size = np.ones(3) + np.random.random(3)
        coords = seeds.from_random(size,N_seeds,cells)
        grid_1 = Grid.from_Voronoi_tessellation(cells,size,coords)
        coords,material = seeds.from_grid(grid_1)
        grid_2 = Grid.from_Voronoi_tessellation(cells,size,coords,material)
        assert (grid_2.material==grid_1.material).all()

    @pytest.mark.parametrize('periodic',[True,False])
    @pytest.mark.parametrize('average',[True,False])
    def test_from_grid_grid(self,periodic,average):
        cells = np.random.randint(10,20,3)
        size  = np.ones(3) + np.random.random(3)
        coords = grid_filters.coordinates0_point(cells,size).reshape(-1,3)
        np.random.shuffle(coords)
        grid_1 = Grid.from_Voronoi_tessellation(cells,size,coords)
        coords,material = seeds.from_grid(grid_1,average=average,periodic=periodic)
        grid_2 = Grid.from_Voronoi_tessellation(cells,size,coords,material)
        assert (grid_2.material==grid_1.material).all()

    @pytest.mark.parametrize('periodic',[True,False])
    @pytest.mark.parametrize('average',[True,False])
    @pytest.mark.parametrize('invert',[True,False])
    def test_from_grid_selection(self,periodic,average,invert):
        cells = np.random.randint(10,20,3)
        N_seeds = np.random.randint(30,300)
        size = np.ones(3) + np.random.random(3)
        coords = seeds.from_random(size,N_seeds,cells)
        grid = Grid.from_Voronoi_tessellation(cells,size,coords)
        selection=np.random.randint(N_seeds)+1
        coords,material = seeds.from_grid(grid,average=average,periodic=periodic,invert=invert,selection=[selection])
        assert selection not in material if invert else (selection==material).all()
