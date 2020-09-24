import pytest
import numpy as np
from scipy.spatial import cKDTree

from damask import seeds

class TestSeeds:

    @pytest.mark.parametrize('grid',[None,np.ones(3,dtype='i')*10])
    def test_from_random(self,grid):
        N_seeds = np.random.randint(30,300)
        size = np.ones(3) + np.random.random(3)
        coords = seeds.from_random(size,N_seeds,grid)
        assert (0<=coords).all() and (coords<size).all()

    def test_from_Poisson_disc(self):
        N_seeds = np.random.randint(30,300)
        N_candidates = N_seeds//15
        distance = np.random.random()
        size = np.ones(3)*distance*N_seeds
        coords = seeds.from_Poisson_disc(size,N_seeds,N_candidates,distance)
        min_dists, _ = cKDTree(coords).query(coords, 2)
        assert (0<= coords).all() and (coords<size).all() and np.min(min_dists[:,1])>=distance
