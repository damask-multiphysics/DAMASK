import pytest
import numpy as np
from scipy.spatial import cKDTree

from damask import seeds
from damask import grid_filters
from damask import GeomGrid


@pytest.mark.parametrize('cells',[None,np.ones(3,dtype='i')*10])
def test_from_random(np_rng,cells):
    N_seeds = np_rng.integers(30,300)
    size = np.ones(3) + np_rng.random(3)
    coords = seeds.from_random(size,N_seeds,cells,rng_seed=np_rng)
    assert (0<=coords).all() and (coords<size).all()

@pytest.mark.parametrize('periodic',[True,False])
def test_from_Poisson_disc(np_rng,periodic):
    N_seeds = np_rng.integers(30,300)
    N_candidates = N_seeds//15
    distance = np_rng.random()
    size = np.ones(3)*distance*N_seeds
    coords = seeds.from_Poisson_disc(size,N_seeds,N_candidates,distance,periodic=periodic)
    min_dists, _ = cKDTree(coords,boxsize=size).query(coords, 2) if periodic else \
                    cKDTree(coords).query(coords, 2)
    assert (0<= coords).all() and (coords<size).all() and np.min(min_dists[:,1])>=distance

@pytest.mark.parametrize('periodic',[True,False])
def test_from_Poisson_disc_invalid(np_rng,periodic):
    N_seeds = np_rng.integers(30,300)
    N_candidates = N_seeds//15
    distance = np_rng.random()
    size = np.ones(3)*distance
    with pytest.raises(ValueError):
        seeds.from_Poisson_disc(size,N_seeds,N_candidates,distance,periodic=periodic)

def test_from_grid_reconstruct(np_rng):
    cells = np_rng.integers(10,20,3)
    N_seeds = np_rng.integers(30,300)
    size = np.ones(3) + np_rng.random(3)
    coords = seeds.from_random(size,N_seeds,cells,rng_seed=np_rng)
    grid_1 = GeomGrid.from_Voronoi_tessellation(cells,size,coords)
    coords,material = seeds.from_grid(grid_1)
    grid_2 = GeomGrid.from_Voronoi_tessellation(cells,size,coords,material)
    assert (grid_2.material==grid_1.material).all()

@pytest.mark.parametrize('periodic',[True,False])
@pytest.mark.parametrize('average',[True,False])
def test_from_grid_grid(np_rng,periodic,average):
    cells = np_rng.integers(10,20,3)
    size  = np.ones(3) + np_rng.random(3)
    coords = grid_filters.coordinates0_point(cells,size).reshape(-1,3)
    np_rng.shuffle(coords)
    grid_1 = GeomGrid.from_Voronoi_tessellation(cells,size,coords)
    coords,material = seeds.from_grid(grid_1,average=average,periodic=periodic)
    grid_2 = GeomGrid.from_Voronoi_tessellation(cells,size,coords,material)
    assert (grid_2.material==grid_1.material).all()

@pytest.mark.parametrize('periodic',[True,False])
@pytest.mark.parametrize('average',[True,False])
@pytest.mark.parametrize('invert',[True,False])
def test_from_grid_selection(np_rng,periodic,average,invert):
    cells = np_rng.integers(10,20,3)
    N_seeds = np_rng.integers(30,300)
    size = np.ones(3) + np_rng.random(3)
    coords = seeds.from_random(size,N_seeds,cells,rng_seed=np_rng)
    grid = GeomGrid.from_Voronoi_tessellation(cells,size,coords)
    selection=np_rng.integers(N_seeds)+1
    coords,material = seeds.from_grid(grid,average=average,periodic=periodic,invert_selection=invert,selection=[selection])
    assert selection not in material if invert else (selection==material).all()
