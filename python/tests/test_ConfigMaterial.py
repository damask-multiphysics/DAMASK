import os

import pytest
import numpy as np

from damask import ConfigMaterial
from damask import Table

@pytest.fixture
def ref_path(ref_path_base):
    """Directory containing reference results."""
    return ref_path_base/'ConfigMaterial'


class TestConfigMaterial:

    @pytest.mark.parametrize('fname',[None,'test.yaml'])
    def test_load_save(self,ref_path,tmp_path,fname):
        reference = ConfigMaterial.load(ref_path/'material.yaml')
        os.chdir(tmp_path)
        if fname is None:
            reference.save()
            new = ConfigMaterial.load('material.yaml')
        else:
            reference.save(fname)
            new = ConfigMaterial.load(fname)
        assert reference == new

    def test_valid_complete(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        assert material_config.is_valid and material_config.is_complete

    def test_invalid_lattice(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        material_config['phase']['Aluminum']['lattice']='fxc'
        assert not material_config.is_valid

    def test_invalid_orientation(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        material_config['material'][0]['constituents'][0]['O']=[0,0,0,0]
        assert not material_config.is_valid

    def test_invalid_fraction(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        material_config['material'][0]['constituents'][0]['v']=.9
        assert not material_config.is_valid

    @pytest.mark.parametrize('item',['homogenization','phase','material'])
    def test_incomplete_missing(self,ref_path,item):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config[item]
        assert not material_config.is_complete

    @pytest.mark.parametrize('item',['O','phase'])
    def test_incomplete_material_constituent(self,ref_path,item):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config['material'][0]['constituents'][0][item]
        assert not material_config.is_complete

    def test_incomplete_material_homogenization(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config['material'][0]['homogenization']
        assert not material_config.is_complete

    def test_incomplete_homogenization_N_constituents(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        for h in material_config['homogenization'].keys():
            del material_config['homogenization'][h]['N_constituents']
        assert not material_config.is_complete

    def test_incomplete_phase_lattice(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config['phase']['Aluminum']['lattice']
        assert not material_config.is_complete

    def test_incomplete_wrong_phase(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        new = material_config.material_rename_phase({'Steel':'FeNbC'})
        assert not new.is_complete

    def test_incomplete_wrong_homogenization(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
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
            assert m['c'] == 1 and c['b'] == 0 and (c['a'] == [i,1]).all()

    def test_constituents(self):
        c = ConfigMaterial._constituents(c=1,v=[2,3])
        assert c[0][0]['c'] == c[1][0]['c']    == 1
        assert c[0][0]['v'] == c[1][0]['v'] -1 ==2

    @pytest.mark.parametrize('constituents',[{'W':1,'X':[2,3]},{'Y':4},{'Z':[5,6]}])
    @pytest.mark.parametrize('a',[[7.,8.],9.])
    @pytest.mark.parametrize('b',['bd',['efg','hi']])
    def test_material_add(self,tmp_path,constituents,a,b):
        len_c = len(ConfigMaterial()._constituents(1,**constituents))
        len_a = len(a) if isinstance(a,list) else 1
        len_b = len(b) if isinstance(b,list) else 1
        m = ConfigMaterial().material_add(constituents,a=a,b=b)
        m.save()
        assert len(m['material']) == np.max([len_a,len_b,len_c])

    @pytest.mark.parametrize('constituents',[{'W':1,'X':np.array([2,3])},{'Y':4},{'Z':np.array([5,6])}])
    @pytest.mark.parametrize('a',[np.array([7,8]),9])
    def test_material_add_np(self,tmp_path,constituents,a):
        len_c = len(ConfigMaterial()._constituents(1,**constituents))
        len_a = len(a) if isinstance(a,np.ndarray) else 1
        m = ConfigMaterial().material_add(constituents,ld=a)
        m.save()
        assert len(m['material']) == np.max([len_a,len_c])

    @pytest.mark.parametrize('constituents',[{'X':np.array([2,3,4,5])},{'Y':4}])
    @pytest.mark.parametrize('a',[np.array([1,2,3]),[4,5,6]])
    @pytest.mark.parametrize('b',[np.array([6.,7.]),[8.,9.]])
    def test_material_add_invalid(self,constituents,a,b):
        with pytest.raises(ValueError):
            ConfigMaterial().material_add(constituents,a=a,u=b)
