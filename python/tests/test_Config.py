import pytest
import numpy as np

from damask import Config

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
        assert Config().is_valid is None

    def test_abstract_is_complete(self):
        assert Config().is_complete is None
