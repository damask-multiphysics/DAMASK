import pytest
import numpy as np

from damask import YAML
from damask import Rotation
from damask import Orientation

def test_init_keyword():
    assert YAML(p=4)['p'] == 4

@pytest.mark.parametrize('config',[{'p':1},'{p: 1}'])
def test_init_config(config):
    assert YAML(config)['p'] == 1

@pytest.mark.parametrize('config',[{'p':1},'{p: 1}'])
def test_init_both(config):
    assert YAML(config,p=2)['p'] == 2

@pytest.mark.parametrize('flow_style',[None,True,False])
def test_load_save_path(tmp_path,flow_style):
    config = YAML()
    config['A'] = 1
    config['B'] = [2,3]
    config.save(tmp_path/'config.yaml',default_flow_style=flow_style)
    assert YAML.load(tmp_path/'config.yaml') == config

def test_load_save_file(tmp_path):
    config = YAML()
    config['A'] = 1
    config['B'] = [2,3]
    with open(tmp_path/'config.yaml','w') as f:
        config.save(f)
    with open(tmp_path/'config.yaml') as f:
        assert YAML.load(f) == config

def test_add_remove():
    dummy = {'hello':'world','foo':'bar'}
    config = YAML()
    config |= dummy
    assert config == YAML() | dummy
    config = config.delete(dummy)
    assert config == YAML()
    assert (config |      dummy ).delete(        'hello'            ) == config | {'foo':'bar'}
    assert (config |      dummy ).delete([       'hello',  'foo'   ]) == config
    assert (config | YAML(dummy)).delete({       'hello':1,'foo':2 }) == config
    assert (config | YAML(dummy)).delete(YAML({'hello':1        }))   == config | {'foo':'bar'}

def test_repr(tmp_path):
    config = YAML()
    config['A'] = 1
    config['B'] = [2,3]
    with open(tmp_path/'config.yaml','w') as f:
        f.write(config.__repr__())
    assert YAML.load(tmp_path/'config.yaml') == config

def test_utf8_minimal_escaping(tmp_path):
    config = YAML()
    config['key'] = 'Ælżβéth 🌸 Михаи́л ᚠᛚᚨ🧙‍♂️'
    repr = config.__repr__()
    assert r'\x' not in repr
    assert r'\u' not in repr

def test_numpy(tmp_path):
    assert YAML({'A':np.ones(3,'i'), 'B':np.ones(1)[0]}).__repr__() == \
           YAML({'A':[1,1,1],        'B':1.0}).__repr__()

def test_abstract_is_valid():
    with pytest.raises(NotImplementedError):
        YAML().is_valid

def test_abstract_is_complete():
    with pytest.raises(NotImplementedError):
        YAML().is_complete

@pytest.mark.parametrize('data',[Rotation.from_random(),Orientation.from_random(lattice='cI')])
def test_rotation_orientation(data):
    assert str(YAML(a=data)) == str(YAML(a=data.as_quaternion()))

def test_initialize():
    yml = """
    a:
        - 1
        - 2
    """
    assert YAML(yml) == YAML('{"a":[1,2]}') == YAML(a=[1,2]) == YAML(dict(a=[1,2]))
