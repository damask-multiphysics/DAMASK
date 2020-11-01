import os

import pytest
import numpy as np

from damask import ConfigMaterial
from damask import Table

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return reference_dir_base/'ConfigMaterial'


class TestConfigMaterial:

    @pytest.mark.parametrize('fname',[None,'test.yaml'])
    def test_load_save(self,reference_dir,tmp_path,fname):
        reference = ConfigMaterial.load(reference_dir/'material.yaml')
        os.chdir(tmp_path)
        if fname is None:
            reference.save()
            new = ConfigMaterial.load('material.yaml')
        else:
            reference.save(fname)
            new = ConfigMaterial.load(fname)
        assert reference == new

    def test_valid_complete(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        assert material_config.is_valid and material_config.is_complete

    def test_invalid_lattice(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        material_config['phase']['Aluminum']['lattice']='fxc'
        assert not material_config.is_valid

    def test_invalid_orientation(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        material_config['material'][0]['constituents'][0]['O']=[0,0,0,0]
        assert not material_config.is_valid

    def test_invalid_fraction(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        material_config['material'][0]['constituents'][0]['fraction']=.9
        assert not material_config.is_valid

    @pytest.mark.parametrize('item',['homogenization','phase','material'])
    def test_incomplete_missing(self,reference_dir,item):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        del material_config[item]
        assert not material_config.is_complete

    @pytest.mark.parametrize('item',['O','phase'])
    def test_incomplete_material_constituent(self,reference_dir,item):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        del material_config['material'][0]['constituents'][0][item]
        assert not material_config.is_complete

    def test_incomplete_material_homogenization(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        del material_config['material'][0]['homogenization']
        assert not material_config.is_complete

    def test_incomplete_phase_lattice(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        del material_config['phase']['Aluminum']['lattice']
        assert not material_config.is_complete

    def test_incomplete_wrong_phase(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        new = material_config.material_rename_phase({'Steel':'FeNbC'})
        assert not new.is_complete

    def test_incomplete_wrong_homogenization(self,reference_dir):
        material_config = ConfigMaterial.load(reference_dir/'material.yaml')
        new = material_config.material_rename_homogenization({'Taylor':'isostrain'})
        assert not new.is_complete

    def test_from_table(self):
        N = np.random.randint(3,10)
        a = np.vstack((np.hstack((np.arange(N),np.arange(N)[::-1])),np.ones(N*2),np.zeros(N*2),np.ones(N*2))).T
        t = Table(a,{'varying':2,'constant':2})
        c = ConfigMaterial.from_table(t,constituents={'a':'varying','b':'1_constant'},c='2_constant')
        assert len(c['material']) == N
        for i,m in enumerate(c['material']):
            c = m['constituents'][0]
            assert m['c'] == 1 and c['b'] == 0 and c['a'] == [i,1]

    def test__constituents(self):
        c = ConfigMaterial._constituents(c=1,v=[2,3])
        assert c[0][0]['c'] == c[1][0]['c']    == 1
        assert c[0][0]['v'] == c[1][0]['v'] -1 ==2
