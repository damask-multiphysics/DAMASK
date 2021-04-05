import pytest
import numpy as np

from damask import Table


@pytest.fixture
def default():
    """Simple Table."""
    x = np.ones((5,13),dtype=float)
    return Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['test data','contains five rows of only ones'])

@pytest.fixture
def ref_path(ref_path_base):
    """Directory containing reference results."""
    return ref_path_base/'Table'

class TestTable:

    def test_repr(self,default):
        print(default)

    @pytest.mark.parametrize('N',[10,40])
    def test_len(self,N):
        assert len(Table(np.random.rand(N,3),{'X':3})) == N

    def test_get_scalar(self,default):
        d = default.get('s')
        assert np.allclose(d,1.0) and d.shape[1:] == (1,)

    def test_get_vector(self,default):
        d = default.get('v')
        assert np.allclose(d,1.0) and d.shape[1:] == (3,)

    def test_get_tensor(self,default):
        d = default.get('F')
        assert np.allclose(d,1.0) and d.shape[1:] == (3,3)

    def test_set(self,default):
        d = default.set('F',np.zeros((5,3,3)),'set to zero').get('F')
        assert np.allclose(d,0.0) and d.shape[1:] == (3,3)

    def test_set_component(self,default):
        d = default.set('F[0,0]',np.zeros((5)),'set to zero').get('F')
        assert np.allclose(d[...,0,0],0.0) and d.shape[1:] == (3,3)

    def test_labels(self,default):
        assert default.labels == ['F','v','s']

    def test_add(self,default):
        d = np.random.random((5,9))
        assert np.allclose(d,default.add('nine',d,'random data').get('nine'))

    def test_isclose(self,default):
        assert default.isclose(default).all()

    def test_allclose(self,default):
        assert default.allclose(default)

    @pytest.mark.parametrize('N',[1,3,4])
    def test_slice(self,default,N):
        assert len(default[:N]) == 1+N
        assert len(default[:N,['F','s']]) == 1+N
        assert default[N:].get('F').shape == (len(default)-N,3,3)
        assert (default[:N,['v','s']].data == default['v','s'][:N].data).all().all()

    @pytest.mark.parametrize('mode',['str','path'])
    def test_write_read(self,default,tmp_path,mode):
        default.save(tmp_path/'default.txt')
        if   mode == 'path':
            new = Table.load(tmp_path/'default.txt')
        elif mode == 'str':
            new = Table.load(str(tmp_path/'default.txt'))
        assert all(default.data==new.data) and default.shapes == new.shapes

    def test_write_read_file(self,default,tmp_path):
        with open(tmp_path/'default.txt','w') as f:
            default.save(f)
        with open(tmp_path/'default.txt') as f:
            new = Table.load(f)
        assert all(default.data==new.data) and default.shapes == new.shapes

    def test_write_invalid_format(self,default,tmp_path):
        with pytest.raises(TypeError):
            default.save(tmp_path/'shouldnotbethere.txt',format='invalid')

    @pytest.mark.parametrize('mode',['str','path'])
    def test_read_ang(self,ref_path,mode):
        if   mode == 'path':
            new = Table.load_ang(ref_path/'simple.ang')
        elif mode == 'str':
            new = Table.load_ang(str(ref_path/'simple.ang'))
        assert new.data.shape == (4,10) and \
               new.labels == ['eu', 'pos', 'IQ', 'CI', 'ID', 'intensity', 'fit']

    def test_read_ang_file(self,ref_path):
        f = open(ref_path/'simple.ang')
        new = Table.load_ang(f)
        assert new.data.shape == (4,10) and \
               new.labels == ['eu', 'pos', 'IQ', 'CI', 'ID', 'intensity', 'fit']

    @pytest.mark.parametrize('fname',['datatype-mix.txt','whitespace-mix.txt'])
    def test_read_strange(self,ref_path,fname):
        with open(ref_path/fname) as f:
            Table.load(f)

    def test_rename_equivalent(self):
        x = np.random.random((5,13))
        t = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        s = t.get('s')
        u = t.rename('s','u').get('u')
        assert np.all(s == u)

    def test_rename_gone(self,default):
        gone = default.rename('v','V')
        assert 'v' not in gone.shapes and 'v' not in gone.data.columns
        with pytest.raises(KeyError):
            gone.get('v')

    def test_delete(self,default):
        delete = default.delete('v')
        assert 'v' not in delete.shapes and 'v' not in delete.data.columns
        with pytest.raises(KeyError):
            delete.get('v')

    def test_join(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        y = np.random.random((5,3))
        b = Table(y,{'u':(3,)},['random test data'])
        c = a.join(b)
        assert np.array_equal(c.get('u'), b.get('u'))

    def test_join_invalid(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        with pytest.raises(KeyError):
            a.join(a)

    def test_append(self):
        x = np.random.random((5,13))
        a = Table(x,{'F':(3,3),'v':(3,),'s':(1,)},['random test data'])
        b = a.append(a)
        assert np.array_equal(b.data[:5].to_numpy(),b.data[5:].to_numpy())

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
        sort   = t.sort_by('s').get('s')
        assert np.all(np.sort(unsort,0)==sort)

    def test_sort_component(self):
        x = np.random.random((5,12))
        t = Table(x,{'F':(3,3),'v':(3,)},['random test data'])
        unsort = t.get('F')[:,1,0]
        sort   = t.sort_by('F[1,0]').get('F')[:,1,0]
        assert np.all(np.sort(unsort,0)==sort)

    def test_sort_revert(self):
        x = np.random.random((5,12))
        t = Table(x,{'F':(3,3),'v':(3,)},['random test data'])
        sort = t.sort_by('F[1,0]',ascending=False).get('F')[:,1,0]
        assert np.all(np.sort(sort,0)==sort[::-1])

    def test_sort(self):
        t = Table(np.array([[0,1,],[2,1,]]),
                  {'v':(2,)},
                  ['test data'])\
            .add('s',np.array(['b','a']))\
            .sort_by('s')
        assert np.all(t.get('v')[:,0] == np.array([2,0]))
