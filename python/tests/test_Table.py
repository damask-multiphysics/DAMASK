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
    return reference_dir_base/'Table'

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

    @pytest.mark.parametrize('mode',['str','path'])
    def test_write_read(self,default,tmpdir,mode):
        default.to_ASCII(tmpdir/'default.txt')
        if   mode == 'path':
            new = Table.from_ASCII(tmpdir/'default.txt')
        elif mode == 'str':
            new = Table.from_ASCII(str(tmpdir/'default.txt'))
        assert all(default.data==new.data) and default.shapes == new.shapes

    def test_write_read_file(self,default,tmpdir):
        with open(tmpdir.join('default.txt'),'w') as f:
            default.to_ASCII(f)
        with open(tmpdir.join('default.txt')) as f:
            new = Table.from_ASCII(f)
        assert all(default.data==new.data) and default.shapes == new.shapes

    def test_write_read_new_style(self,default,tmpdir):
        with open(tmpdir.join('new_style.txt'),'w') as f:
            default.to_ASCII(f,new_style=True)
        with open(tmpdir.join('new_style.txt')) as f:
            new = Table.from_ASCII(f)
        assert all(default.data==new.data) and default.shapes == new.shapes


    @pytest.mark.parametrize('mode',['str','path'])
    def test_read_ang(self,reference_dir,mode):
        if   mode == 'path':
            new = Table.from_ang(reference_dir/'simple.ang')
        elif mode == 'str':
            new = Table.from_ang(str(reference_dir/'simple.ang'))
        assert new.data.shape == (4,10) and \
               new.labels == ['eu', 'pos', 'IQ', 'CI', 'ID', 'intensity', 'fit']

    def test_read_ang_file(self,reference_dir):
        f = open(reference_dir/'simple.ang')
        new = Table.from_ang(f)
        assert new.data.shape == (4,10) and \
               new.labels == ['eu', 'pos', 'IQ', 'CI', 'ID', 'intensity', 'fit']

    @pytest.mark.parametrize('fname',['datatype-mix.txt','whitespace-mix.txt'])
    def test_read_strange(self,reference_dir,fname):
        with open(reference_dir/fname) as f:
            Table.from_ASCII(f)

    def test_set(self,default):
        default.set('F',np.zeros((5,3,3)),'set to zero')
        d=default.get('F')
        assert np.allclose(d,0.0) and d.shape[1:] == (3,3)

    def test_set_component(self,default):
        default.set('1_F',np.zeros((5)),'set to zero')
        d=default.get('F')
        assert np.allclose(d[...,0,0],0.0) and d.shape[1:] == (3,3)

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
        assert 'v' not in default.shapes and 'v' not in default.data.columns
        with pytest.raises(KeyError):
            default.get('v')

    def test_delete(self,default):
        default.delete('v')
        assert 'v' not in default.shapes and 'v' not in default.data.columns
        with pytest.raises(KeyError):
            default.get('v')

    def test_join(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        y = np.random.random((5,3))
        b = Table(y,{'u':(3,)},['random test data'])
        a.join(b)
        assert np.array_equal(a.get('u'), b.get('u'))

    def test_join_invalid(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        with pytest.raises(KeyError):
            a.join(a)

    def test_append(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        a.append(a)
        assert np.array_equal(a.data[:5].to_numpy(),a.data[5:].to_numpy())

    def test_append_invalid(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        b = Table(x,{'F':(3,3),'u':(3,),'s':(1,)},['random test data'])
        with pytest.raises(KeyError):
            a.append(b)

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
        assert np.all(t.get('1_v') == np.array([2,0]).reshape(2,1))
