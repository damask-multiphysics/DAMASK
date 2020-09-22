import os

import pytest

from damask import Material

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return reference_dir_base/'Material'


class TestMaterial:
    
    @pytest.mark.parametrize('fname',[None,'test.yaml'])
    def test_load_save(self,reference_dir,tmp_path,fname):
        reference = Material.load(reference_dir/'material.yaml')
        os.chdir(tmp_path)
        if fname is None:
            reference.save()
            new = Material.load('material.yaml')
        else:
            reference.save(fname)
            new = Material.load(fname)
        assert reference == new

    def test_valid_complete(self,reference_dir):
        material_config = Material.load(reference_dir/'material.yaml')
        assert material_config.is_valid and material_config.is_complete

    def test_invalid_lattice(self,reference_dir):
        material_config = Material.load(reference_dir/'material.yaml')
        material_config['phase']['Aluminum']['lattice']='fxc'
        assert not material_config.is_valid
    
    def test_invalid_orientation(self,reference_dir):
        material_config = Material.load(reference_dir/'material.yaml')
        material_config['microstructure'][0]['constituents'][0]['orientation']=[0,0,0,0]
        assert not material_config.is_valid
    
    def test_invalid_fraction(self,reference_dir):
        material_config = Material.load(reference_dir/'material.yaml')
        material_config['microstructure'][0]['constituents'][0]['fraction']=.9
        assert not material_config.is_valid
    

    @pytest.mark.parametrize('item',['homogenization','phase','microstructure'])
    def test_incomplete_missing(self,reference_dir,item):
        material_config = Material.load(reference_dir/'material.yaml')
        del material_config[item]
        assert not material_config.is_complete
    
    def test_incomplete_wrong_phase(self,reference_dir):
        material_config = Material.load(reference_dir/'material.yaml')
        new = material_config.microstructure_rename_phase({'Steel':'FeNbC'})
        assert not new.is_complete
    
    def test_incomplete_wrong_homogenization(self,reference_dir):
        material_config = Material.load(reference_dir/'material.yaml')
        new = material_config.microstructure_rename_homogenization({'Taylor':'isostrain'})
        assert not new.is_complete
