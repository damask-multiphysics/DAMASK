import pytest
import numpy as np

from damask import Table


@pytest.fixture
def default():
    """Simple Table."""
    x = np.ones((5,13),dtype=float)
    return Table({'F':(3,3),'v':(3,),'s':(1,)},
                 x,
                 ['test data','contains five rows of only ones'])

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Table'


def test_repr(default):
    print(default)

@pytest.mark.parametrize('N',[10,40])
def test_len(np_rng,N):
    assert len(Table({'X':3},np_rng.random((N,3)))) == N

def test_get_scalar(default):
    d = default.get('s')
    assert np.allclose(d,1.0) and d.shape[1:] == (1,)

def test_get_vector(default):
    d = default.get('v')
    assert np.allclose(d,1.0) and d.shape[1:] == (3,)

def test_get_tensor(default):
    d = default.get('F')
    assert np.allclose(d,1.0) and d.shape[1:] == (3,3)

def test_empty_init():
    N = 3
    D = dict(
             scal=np.arange(10),
             vctr=np.arange(10*N).reshape((10,N)),
            )
    t = Table()
    for label,data in D.items():
        t = t.set(label,data)
    assert np.allclose(t.get('scal').flatten()*3,t.get('vctr')[:,0])

def test_set_tensor(default):
    d = default.set('F',np.zeros((5,3,3)),'set to zero').get('F')
    assert np.allclose(d,0.0) and d.shape[1:] == (3,3)

def test_set_scalar(np_rng,default):
    d = np_rng.random((5,9))
    assert np.allclose(d,
                        default
                        .set('nine',d,'random data')
                        .get('nine'))

def test_set_overwrite(np_rng,default):
    d = np_rng.random((5,9))
    assert np.allclose(d,
                       default
                              .set('nine',np.ones((5,9)),'zeros')
                              .set('nine',d,'random data')
                              .get('nine'))

def test_set_component(default):
    d = default.set('F[0,0]',np.zeros((5)),'set to zero').get('F')
    assert np.allclose(d[...,0,0],0.0) and d.shape[1:] == (3,3)

def test_labels(default):
    assert default.labels == ['F','v','s']

def test_isclose(default):
    assert default.isclose(default).all()

def test_allclose(default):
    assert default.allclose(default)

@pytest.mark.parametrize('N',[1,3,4])
def test_slice(np_rng,default,N):
    mask = np_rng.choice([True,False],len(default))
    arr = np_rng.choice(len(default),N)
    assert len(default[N]) == 1
    assert len(default[1:N]) == N
    assert len(default[arr]) == N
    assert len(default[:N]) == 1+N
    assert len(default[:N,['F','s']]) == 1+N
    assert len(default[mask,['F','s']]) == np.count_nonzero(mask)
    assert default[mask,['F','s']] == default[mask][['F','s']] == default[['F','s']][mask]
    assert default[np.logical_not(mask),['F','s']] != default[mask][['F','s']]
    assert default[N:].get('F').shape == (len(default)-N,3,3)
    assert default[:N,['v','s']].data.equals(default['v','s'][:N].data)

@pytest.mark.parametrize('mode',['str','path','file'])
def test_write_read_mode(default,tmp_path,mode):
    path = tmp_path/'default.txt'
    if mode == 'file':
        default.save(open(path,'w'))
        new = Table.load(open(path))
    else:
        default.save(str(path) if mode == 'str' else path)
        new = Table.load(str(path) if mode == 'str' else path)
    assert all(default.data == new.data) and default.shapes == new.shapes

def test_write_read_file(default,tmp_path):
    with open(tmp_path/'default.txt','w') as f:
        default.save(f)
    with open(tmp_path/'default.txt') as f:
        new = Table.load(f)
    assert all(default.data == new.data) and default.shapes == new.shapes

def test_write_invalid_format(default,tmp_path):
    with pytest.raises(TypeError):
        default.save(tmp_path/'shouldnotbethere.txt',format='invalid')

@pytest.mark.parametrize('mode',['str','path','file'])
def test_read_ang_mode(res_path,mode):
    where = res_path/'simple.ang'
    fname = {'path': where,
             'str':  str(where),
             'file': open(where)}[mode]
    new = Table.load_ang(fname)
    assert new.data.shape == (4,10) and \
           new.labels == ['eu', 'pos', 'IQ', 'CI', 'ID', 'intensity', 'fit']

def test_read_ang_shapes(res_path):
    new = Table.load_ang(str(res_path/'simple.ang'),shapes={})
    assert new.data.shape == (4,10) and new.labels == ['unknown']

def test_save_ang(res_path,tmp_path):
    orig = Table.load_ang(res_path/'simple.ang')
    orig.save(tmp_path/'simple.ang',with_labels=False)
    saved = Table.load_ang(tmp_path/'simple.ang')
    assert saved == orig

def test_ang_no_header(res_path,tmp_path):
    orig = Table.load_ang(res_path/'simple.ang')
    with open(res_path/'simple.ang') as f_in, open(tmp_path/'no_header.ang','w') as f_out:
        for l in f_in:
            if not l.startswith('#'): f_out.write(l)
    no_header = Table.load_ang(tmp_path/'no_header.ang')
    assert no_header == orig


@pytest.mark.parametrize('fname',['datatype-mix.txt','whitespace-mix.txt'])
def test_read_strange(res_path,fname):
    with open(res_path/fname) as f:
        Table.load(f)

def test_rename_equivalent(np_rng):
    x = np_rng.random((5,13))
    t = Table({'F':(3,3),'v':(3,),'s':(1,)},x,['random test data'])
    s = t.get('s')
    u = t.rename('s','u').get('u')
    assert np.all(s == u)

def test_rename_gone(default):
    gone = default.rename('v','V')
    assert 'v' not in gone.shapes and 'v' not in gone.data.columns
    with pytest.raises(KeyError):
        gone.get('v')

def test_delete(default):
    delete = default.delete('v')
    assert 'v' not in delete.shapes and 'v' not in delete.data.columns
    with pytest.raises(KeyError):
        delete.get('v')

def test_join(np_rng):
    x = np_rng.random((5,13))
    a = Table({'F':(3,3),'v':(3,),'s':(1,)},x,['random test data'])
    y = np_rng.random((5,3))
    b = Table({'u':(3,)},y,['random test data'])
    c = a.join(b)
    assert np.array_equal(c.get('u'), b.get('u'))

def test_join_invalid(np_rng):
    x = np_rng.random((5,13))
    a = Table({'F':(3,3),'v':(3,),'s':(1,)},x,['random test data'])
    with pytest.raises(KeyError):
        a.join(a)

def test_append(np_rng):
    x = np_rng.random((5,13))
    a = Table({'F':(3,3),'v':(3,),'s':(1,)},x,['random test data'])
    b = a.append(a)
    assert np.array_equal(b.data[:5].to_numpy(),b.data[5:].to_numpy())

def test_append_invalid(np_rng):
    x = np_rng.random((5,13))
    a = Table({'F':(3,3),'v':(3,),'s':(1,)},x,['random test data'])
    b = Table({'F':(3,3),'u':(3,),'s':(1,)},x,['random test data'])
    with pytest.raises(KeyError):
        a.append(b)

def test_invalid_initialization(np_rng):
    x = np_rng.random((5,10))
    with pytest.raises(ValueError):
        Table({'F':(3,3)},x)

def test_invalid_set(default):
    x = default.get('v')
    with pytest.raises(ValueError):
        default.set('F',x,'does not work')

def test_invalid_get(default):
    with pytest.raises(KeyError):
        default.get('n')

def test_sort_scalar(np_rng):
    x = np_rng.random((5,13))
    t = Table({'F':(3,3),'v':(3,),'s':(1,)},x,['random test data'])
    unsort = t.get('s')
    sort   = t.sort_by('s').get('s')
    assert np.all(np.sort(unsort,0)==sort)

def test_sort_component(np_rng):
    x = np_rng.random((5,12))
    t = Table({'F':(3,3),'v':(3,)},x,['random test data'])
    unsort = t.get('F')[:,1,0]
    sort   = t.sort_by('F[1,0]').get('F')[:,1,0]
    assert np.all(np.sort(unsort,0)==sort)

def test_sort_revert(np_rng):
    x = np_rng.random((5,12))
    t = Table({'F':(3,3),'v':(3,)},x,['random test data'])
    sort = t.sort_by('F[1,0]',ascending=False).get('F')[:,1,0]
    assert np.all(np.sort(sort,0)==sort[::-1])

def test_sort():
    t = Table({'v':(2,)},
              np.array([[0,1,],[2,1,]]),
              ['test data']
             )\
             .set('s',np.array(['b','a']))\
             .sort_by('s')
    assert np.all(t.get('v')[:,0] == np.array([2,0]))

def test_unique(np_rng):
    x = np_rng.choice([True,False],size=(np_rng.integers(12,18),2))
    t = Table({'a':(1,),'b':(1,)},x,['random choices'])
    u,idx,inv,cnt = t.unique(
        return_index=True,
        return_inverse=True,
        return_counts=True,
        )
    assert len(inv.shape) == 1, f'return_inverse must be 1D but got {inv.shape=}'
    assert u.allclose(t[idx])
    assert t.allclose(u[inv])
    assert (cnt[np.lexsort(u.data.to_numpy()[:, ::-1].T)] == np.unique(x,return_counts=True,axis=0)[1]).all()
