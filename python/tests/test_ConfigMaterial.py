import sys
import os
import pytest
import numpy as np

from damask import ConfigMaterial
from damask import Table
from damask import Rotation
from damask import Grid

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'ConfigMaterial'


class TestConfigMaterial:

    def test_init_empty(self):
        c = ConfigMaterial()
        assert len(c) == 3
        assert c['homogenization'] == {}
        assert c['phase'] == {}
        assert c['material'] == []

    def test_init_d(self):
        c = ConfigMaterial(config={'phase':4})
        assert len(c) == 1
        assert c['phase'] == 4

    @pytest.mark.parametrize('kwargs',[{'homogenization':{'SX':{}}},
                                       {'phase':{'Aluminum':{}}},
                                       {'material':[{'A':1},{'B':2}]}])
    def test_init_some(self,kwargs):
        c = ConfigMaterial(**kwargs)
        assert len(c) == 3
        for k,v in kwargs.items():
            if k in kwargs: assert v == kwargs[k]

    @pytest.mark.parametrize('fname',[None,'test.yaml'])
    def test_load_save(self,res_path,tmp_path,fname):
        reference = ConfigMaterial.load(res_path/'material.yaml')
        os.chdir(tmp_path)
        if fname is None:
            reference.save()
            new = ConfigMaterial.load('material.yaml')
        else:
            reference.save(fname)
            new = ConfigMaterial.load(fname)
        assert reference == new

    def test_valid_complete(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        assert material_config.is_valid and material_config.is_complete

    def test_invalid_lattice(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        material_config['phase']['Aluminum']['lattice']='fxc'
        assert not material_config.is_valid

    def test_invalid_orientation(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        material_config['material'][0]['constituents'][0]['O']=[0,0,0,0]
        assert not material_config.is_valid

    @pytest.mark.xfail(sys.platform == 'win32', reason='utf8 "not equal" might cause trouble')
    def test_invalid_fraction(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        material_config['material'][0]['constituents'][0]['v']=.9
        assert not material_config.is_valid

    @pytest.mark.parametrize('item',['homogenization','phase','material'])
    def test_incomplete_missing(self,res_path,item):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        del material_config[item]
        assert not material_config.is_complete

    @pytest.mark.parametrize('item',['O','phase'])
    def test_incomplete_material_constituent(self,res_path,item):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        del material_config['material'][0]['constituents'][0][item]
        assert not material_config.is_complete

    def test_incomplete_material_homogenization(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        del material_config['material'][0]['homogenization']
        assert not material_config.is_complete

    def test_incomplete_wrong_phase(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        new = material_config.material_rename_phase({'Steel':'FeNbC'})
        assert not new.is_complete

    def test_incomplete_wrong_homogenization(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        new = material_config.material_rename_homogenization({'Taylor':'isostrain'})
        assert not new.is_complete

    def test_empty_phase(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        material_config['phase'] = None
        assert not material_config.is_complete

    def test_empty_homogenization(self,res_path):
        material_config = ConfigMaterial.load(res_path/'material.yaml')
        material_config['homogenization'] = None
        assert not material_config.is_complete

    def test_from_table(self):
        N = np.random.randint(3,10)
        a = np.vstack((np.hstack((np.arange(N),np.arange(N)[::-1])),
                       np.zeros(N*2),np.ones(N*2),np.zeros(N*2),np.zeros(N*2),
                       np.ones(N*2),
                      )).T
        t = Table({'varying':1,'constant':4,'ones':1},a)
        c = ConfigMaterial.from_table(t,**{'phase':'varying','O':'constant','homogenization':'ones'})
        assert len(c['material']) == N
        for i,m in enumerate(c['material']):
            assert m['homogenization'] == 1 and (m['constituents'][0]['O'] == [0,1,0,0]).all()

    def test_updated_dicts(self,res_path):
        m1 = ConfigMaterial().material_add(phase=['Aluminum'],O=[1.0,0.0,0.0,0.0],homogenization='SX')
        m2 = ConfigMaterial.load(res_path/'material.yaml')
        for k in m2['phase']:
            m2 = m2.material_add(phase=[k],O=[1.0,0.0,0.0,0.0],homogenization='SX')
            assert not m2['phase'].get(k) is None
        assert m1['phase'].get('Aluminum') is None
        assert m1['homogenization'].get('SX') is None

    def test_from_table_with_constant(self):
        N = np.random.randint(3,10)
        a = np.vstack((np.hstack((np.arange(N),np.arange(N)[::-1])),
                       np.zeros(N*2),np.ones(N*2),np.zeros(N*2),np.zeros(N*2),
                       np.ones(N*2),
                      )).T
        t = Table({'varying':1,'constant':4,'ones':1},a)
        c = ConfigMaterial.from_table(t,**{'phase':'varying','O':'constant','homogenization':1})
        assert len(c['material']) == N
        for i,m in enumerate(c['material']):
            assert m['homogenization'] == 1 and (m['constituents'][0]['O'] == [0,1,0,0]).all()

    @pytest.mark.parametrize('N,n,kw',[
                                        (1,1,{'phase':'Gold',
                                              'O':[1,0,0,0],
                                              'V_e':np.eye(3),
                                              'homogenization':'SX'}),
                                        (3,1,{'phase':'Gold',
                                              'O':Rotation.from_random(3),
                                              'V_e':np.broadcast_to(np.eye(3),(3,3,3)),
                                              'homogenization':'SX'}),
                                        (2,3,{'phase':np.broadcast_to(['a','b','c'],(2,3)),
                                              'O':Rotation.from_random((2,3)),
                                              'V_e':np.broadcast_to(np.eye(3),(2,3,3,3)),
                                              'homogenization':['SX','PX']}),
                                        ])
    def test_material_add(self,kw,N,n):
        m = ConfigMaterial().material_add(**kw)
        assert len(m['material']) == N
        assert len(m['material'][0]['constituents']) == n

    @pytest.mark.parametrize('shape',[(),(4,),(5,2)])
    @pytest.mark.parametrize('kw',[{'V_e':np.random.rand(3,3)},
                                   {'O':np.random.rand(4)},
                                   {'v':np.array(2)}])
    def test_material_add_invalid(self,kw,shape):
        kw = {arg:np.broadcast_to(val,shape+val.shape) for arg,val in kw.items()}
        with pytest.raises(ValueError):
            ConfigMaterial().material_add(**kw)

    @pytest.mark.parametrize('v',[2,np.ones(3)*2,np.ones((2,2))])
    def test_material_add_invalid_v(self,v):
        with pytest.raises(ValueError):
            ConfigMaterial().material_add(v=v)

    @pytest.mark.parametrize('cell_ensemble_data',[None,'CellEnsembleData'])
    def test_load_DREAM3D(self,res_path,cell_ensemble_data):
        grain_c = ConfigMaterial.load_DREAM3D(res_path/'2phase_irregularGrid.dream3d','Grain Data',
                  cell_ensemble_data = cell_ensemble_data)
        point_c = ConfigMaterial.load_DREAM3D(res_path/'2phase_irregularGrid.dream3d',
                  cell_ensemble_data = cell_ensemble_data)

        assert point_c.is_valid and grain_c.is_valid and \
               len(point_c['material'])+1 == len(grain_c['material'])

        grain_m = Grid.load_DREAM3D(res_path/'2phase_irregularGrid.dream3d','FeatureIds').material.flatten()
        point_m = Grid.load_DREAM3D(res_path/'2phase_irregularGrid.dream3d').material.flatten()

        for i in np.unique(point_m):
            j = int(grain_m[(point_m==i).nonzero()[0][0]])
            assert np.allclose(point_c['material'][i]['constituents'][0]['O'],
                               grain_c['material'][j]['constituents'][0]['O'])
            assert point_c['material'][i]['constituents'][0]['phase'] == \
                   grain_c['material'][j]['constituents'][0]['phase']


    def test_load_DREAM3D_reference(self,tmp_path,res_path,update):
        cur = ConfigMaterial.load_DREAM3D(res_path/'measured.dream3d')
        ref = ConfigMaterial.load(res_path/'measured.material.yaml')
        if update:
            cur.save(res_path/'measured.material.yaml')
        for i,m in enumerate(ref['material']):
            assert Rotation(m['constituents'][0]['O']).isclose(Rotation(cur['material'][i]['constituents'][0]['O']))
        assert cur.is_valid and cur['phase'] == ref['phase'] and cur['homogenization'] == ref['homogenization']
