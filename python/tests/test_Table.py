import os

import pytest
import numpy as np

from damask import Table


@pytest.fixture
def default():
    """Simple Table."""
    x = np.ones((5,13),dtype=float)
    return Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['test data','contains only ones'])

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Table')

class TestTable:
    
    def test_get_scalar(self,default):
        d = default.get('s')
        assert np.allclose(d,1.0) and d.shape[1:] == (1,)

    def test_get_vector(self,default):
        d = default.get('v')
        assert np.allclose(d,1.0) and d.shape[1:] == (3,)

    def test_get_tensor(self,default):
        d = default.get('F')
        assert np.allclose(d,1.0) and d.shape[1:] == (3,3) 

    def test_get_component(self,default):
        d = default.get('5_F')
        assert np.allclose(d,1.0) and d.shape[1:] == (1,)
  
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

    @pytest.mark.parametrize('fname',['datatype-mix.txt','whitespace-mix.txt'])
    def test_read_strange(self,reference_dir,fname):
        with open(os.path.join(reference_dir,fname)) as f:
            Table.from_ASCII(f)
    
    def test_set(self,default):
        default.set('F',np.zeros((5,3,3)),'set to zero')
        d=default.get('F')
        assert np.allclose(d,0.0) and d.shape[1:] == (3,3)

    def test_labels(self,default):
        assert default.labels == ['F','v','s']
        
    def test_add(self,default):
        d = np.random.random((5,9))
        default.add('nine',d,'random data')
        assert np.allclose(d,default.get('nine'))

    def test_rename_equivalent(self):
        x = np.random.random((5,13))
        t = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        s = t.get('s')
        t.rename('s','u')
        u = t.get('u')
        assert np.all(s == u)

    def test_rename_gone(self,default):
        default.rename('v','V')
        with pytest.raises(KeyError):
            default.get('v')

    def test_delete(self,default):
        default.delete('v')
        with pytest.raises(KeyError):
            default.get('v')


    def test_invalid_initialization(self):
        x = np.random.random((5,10))
        with pytest.raises(ValueError):
            Table(x,{'F':(3,3)})

    def test_invalid_set(self,default):
        x = default.get('v')
        with pytest.raises(ValueError):
            default.set('F',x,'does not work')

    def test_invalid_get(self,default):
        with pytest.raises(KeyError):
            default.get('n')

    def test_sort_scalar(self):
        x = np.random.random((5,13))
        t = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        unsort = t.get('s')
        t.sort_by('s')
        sort   = t.get('s')
        assert np.all(np.sort(unsort,0)==sort)

    def test_sort_component(self):
        x = np.random.random((5,12))
        t = Table(x,{'F':(3,3),'v':(3,)},['random test data'])
        unsort = t.get('4_F')
        t.sort_by('4_F')
        sort = t.get('4_F')
        assert np.all(np.sort(unsort,0)==sort)

    def test_sort_revert(self):
        x = np.random.random((5,12))
        t = Table(x,{'F':(3,3),'v':(3,)},['random test data'])
        t.sort_by('4_F',ascending=False)
        sort = t.get('4_F')
        assert np.all(np.sort(sort,0)==sort[::-1,:])

    def test_sort(self):
        t = Table(np.array([[0,1,],[2,1,]]),
                  {'v':(2,)},
                  ['test data'])
        t.add('s',np.array(['b','a']))
        t.sort_by('s')
        assert np.all(t.get('1_v') == np.array([2,0]).reshape((2,1)))
