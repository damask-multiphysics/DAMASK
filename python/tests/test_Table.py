import pytest
import numpy as np

from damask import Table

@pytest.fixture
def default():
    """Simple Table."""
    x = np.ones((5,13))
    return Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['test data','contains only ones'])


class TestTable:
    
    def test_get_tensor(self,default):
        d = default.get_array('F')
        assert np.allclose(d,1.0) and d.shape[1:] == (3,3) 

    def test_get_vector(self,default):
        d = default.get_array('v')
        assert np.allclose(d,1.0) and d.shape[1:] == (3,)
  
    def test_write_read_str(self,default,tmpdir):
        default.to_ASCII(str(tmpdir.join('default.txt')))
        new = Table.from_ASCII(str(tmpdir.join('default.txt')))
        assert all(default.data==new.data)

    def test_write_read_file(self,default,tmpdir):
        with open(tmpdir.join('default.txt'),'w') as f:
            default.to_ASCII(f)
        with open(tmpdir.join('default.txt')) as f:
            new = Table.from_ASCII(f)
        assert all(default.data==new.data)
    
    def test_set_array(self,default):
        default.set_array('F',np.zeros((5,3,3)),'set to zero')
        d=default.get_array('F')
        assert np.allclose(d,0.0) and d.shape[1:] == (3,3)

    def test_get_labels(self,default):
        assert default.get_labels() == ['F','v','s']
        
    def test_add_array(self,default):
        d = np.random.random((5,9))
        default.add_array('nine',d,'random data')
        assert np.allclose(d,default.get_array('nine'))


    def test_invalid_initialization(self,default):
        x = default.get_array('v')
        with pytest.raises(IndexError):
            Table(x,{'F':(3,3)})

    def test_invalid_set(self,default):
        x = default.get_array('v')
        with pytest.raises(ValueError):
            default.set_array('F',x,'does not work')
