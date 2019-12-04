import shutil
import os

import pytest
import numpy as np

from damask import DADF5
from damask import mechanics

@pytest.fixture
def default(tmp_path,reference_dir):
    """Small DADF5 file in temp location for modification."""
    fname = '12grains6x7x8_tensionY.hdf5'
    shutil.copy(os.path.join(reference_dir,fname),tmp_path)
    f = DADF5(os.path.join(tmp_path,fname))
    f.set_by_time(20.0,20.0)
    return f

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'DADF5')


class TestDADF5:

    def test_time_increments(self,default):
        shape = default.read_dataset(default.get_dataset_location('F'),0).shape
        default.set_by_time(0.0,20.0)
        for i in default.iter_visible('increments'):
            assert shape == default.read_dataset(default.get_dataset_location('F'),0).shape

    def test_add_deviator(self,default):
        default.add_deviator('P')
        loc = {'P'  :default.get_dataset_location('P'),
               's_P':default.get_dataset_location('s_P')}
        in_memory = mechanics.deviatoric_part(default.read_dataset(loc['P'],0))
        in_file   = default.read_dataset(loc['s_P'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_Cauchy(self,default):
        default.add_Cauchy('P','F')
        loc = {'F':    default.get_dataset_location('F'),
               'P':    default.get_dataset_location('P'),
               'sigma':default.get_dataset_location('sigma')}
        in_memory = mechanics.Cauchy(default.read_dataset(loc['F'],0),
                                     default.read_dataset(loc['P'],0))
        in_file   = default.read_dataset(loc['sigma'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_determinant(self,default):
        default.add_determinant('P')
        loc = {'P':     default.get_dataset_location('P'),
               'det(P)':default.get_dataset_location('det(P)')}
        in_memory = np.linalg.det(default.read_dataset(loc['P'],0))
        in_file   = default.read_dataset(loc['det(P)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_norm(self,default):
        default.add_norm('F',1)
        loc = {'F':    default.get_dataset_location('F'),
               '|F|_1':default.get_dataset_location('|F|_1')}
        in_memory = np.linalg.norm(default.read_dataset(loc['F'],0),ord=1,axis=(1,2),keepdims=True)
        in_file   = default.read_dataset(loc['|F|_1'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_absolute(self,default):
        default.add_absolute('Fe')
        loc = {'Fe':   default.get_dataset_location('Fe'),
               '|Fe|': default.get_dataset_location('|Fe|')}
        in_memory = np.abs(default.read_dataset(loc['Fe'],0))
        in_file   = default.read_dataset(loc['|Fe|'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_determinant(self,default):
        default.add_determinant('P')
        loc = {'P':      default.get_dataset_location('P'),
               'det(P)': default.get_dataset_location('det(P)')}
        in_memory = np.linalg.det(default.read_dataset(loc['P'],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc['det(P)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_spherical(self,default):
        default.add_spherical('P')
        loc = {'P':   default.get_dataset_location('P'),
               'p_P': default.get_dataset_location('p_P')}
        in_memory = mechanics.spherical_part(default.read_dataset(loc['P'],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc['p_P'],0)
        assert np.allclose(in_memory,in_file)
