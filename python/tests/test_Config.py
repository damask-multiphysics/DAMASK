import pytest
import numpy as np

from damask import Config
from damask import Rotation
from damask import Orientation

class TestConfig:

    @pytest.mark.parametrize('flow_style',[None,True,False])
    def test_load_save_str(self,tmp_path,flow_style):
        config = Config()
        config['A'] = 1
        config['B'] = [2,3]
        config.save(tmp_path/'config.yaml',default_flow_style=flow_style)
        assert Config.load(tmp_path/'config.yaml') == config

    def test_load_save_file(self,tmp_path):
        config = Config()
        config['A'] = 1
        config['B'] = [2,3]
        with open(tmp_path/'config.yaml','w') as f:
            config.save(f)
        with open(tmp_path/'config.yaml') as f:
            assert Config.load(f) == config

    def test_add_remove(self):
        dummy = {'hello':'world','foo':'bar'}
        config = Config()
        config |= dummy
        assert config == Config() | dummy
        config = config.delete(dummy)
        assert config == Config()
        assert (config |        dummy ).delete(        'hello'            ) == config | {'foo':'bar'}
        assert (config |        dummy ).delete([       'hello',  'foo'   ]) == config
        assert (config | Config(dummy)).delete({       'hello':1,'foo':2 }) == config
        assert (config | Config(dummy)).delete(Config({'hello':1        })) == config | {'foo':'bar'}


    def test_repr(self,tmp_path):
        config = Config()
        config['A'] = 1
        config['B'] = [2,3]
        with open(tmp_path/'config.yaml','w') as f:
            f.write(config.__repr__())
        assert Config.load(tmp_path/'config.yaml') == config

    def test_numpy(self,tmp_path):
        assert Config({'A':np.ones(3,'i')}).__repr__() == Config({'A':[1,1,1]}).__repr__()

    def test_abstract_is_valid(self):
        with pytest.raises(NotImplementedError):
            Config().is_valid

    def test_abstract_is_complete(self):
        with pytest.raises(NotImplementedError):
            Config().is_complete

    @pytest.mark.parametrize('data',[Rotation.from_random(),Orientation.from_random(lattice='cI')])
    def test_rotation_orientation(self,data):
        assert str(Config(a=data)) == str(Config(a=data.as_quaternion()))
