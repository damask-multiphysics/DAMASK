import sys
import random
import pydoc

import pytest
import numpy as np
from scipy import stats
import h5py

from damask import util


@pytest.mark.xfail(sys.platform == 'win32', reason='echo is not a Windows command')
def test_run_direct():
    out,err = util.run('echo test')
    assert out=='test\n' and err==''

@pytest.mark.xfail(sys.platform == 'win32', reason='echo is not a Windows command')
def test_run_env():
    out,err = util.run('sh -c "echo $test_for_execute"',env={'test_for_execute':'test'})
    assert out=='test\n' and err==''

@pytest.mark.xfail(sys.platform == 'win32', reason='false is not a Windows command')
def test_run_runtime_error():
    with pytest.raises(RuntimeError):
        util.run('false')

@pytest.mark.parametrize('input,glue,quote,output',
                        [
                         (None,'',False,'None'),
                         ([None,None],'\n',False,'None\nNone'),
                         ([-0.5,0.5],'=',False,'-0.5=0.5'),
                         ([1,2,3],'_',False,'1_2_3'),
                         ([1,2,3],'/',True,'"1"/"2"/"3"'),
                        ])
def test_srepr(input,glue,quote,output):
    assert output == util.srepr(input,glue,quote)


@pytest.mark.parametrize('N',[5,6,7,8,9,10,11,12,13,14,15,16,20,30,40,50])
@pytest.mark.parametrize('input,output',
                        [
                         ([0,-2],[0,-1]),
                         ([-0.5,0.5],[-1,1]),
                         ([1./2.,1./3.],[3,2]),
                         ([2./3.,1./2.,1./3.],[4,3,2]),
                         ([0.666666666666,-0.33333333333,-0.33333],[2,-1,-1]),
                         ([1./3., 1./4., 1./22],[536870912, 402653184,  73209669]),
                        ])
def test_scale2coprime(input,output,N):
    res = util.scale_to_coprime(input,N)
    assert np.allclose(res/np.max(np.abs(res)),output/np.max(np.abs(output)),atol=1e-2,rtol=0)


@pytest.mark.parametrize('rv',[stats.rayleigh(),stats.weibull_min(1.2),stats.halfnorm(),stats.pareto(2.62)])
def test_hybridIA_distribution(rv):
    bins = np.linspace(0,10,100000)
    centers = (bins[1:]+bins[:-1])/2
    N_samples = bins.shape[0]-1000
    dist = rv.pdf(centers)
    selected = util.hybrid_IA(dist,N_samples)
    dist_sampled = np.histogram(centers[selected],bins)[0]/N_samples*np.sum(dist)
    assert np.sqrt(((dist - dist_sampled) ** 2).mean()) < .025 and selected.shape[0]==N_samples

def test_hybridIA_constant(np_rng):
    N_bins = np_rng.integers(20,400)
    m = np_rng.integers(1,20)
    N_samples = m * N_bins
    dist = np.ones(N_bins)*np_rng.random()
    assert np.all(np.sort(util.hybrid_IA(dist,N_samples))==np.arange(N_samples).astype(int)//m)

def test_hybridIA_linear(np_rng):
    N_points = np_rng.integers(10,200)
    m = np_rng.integers(1,20)
    dist = np.arange(N_points)
    N_samples = m * np.sum(dist)
    assert np.all(np.bincount(util.hybrid_IA(dist*np_rng.random(),N_samples)) == dist*m)


@pytest.mark.parametrize('point,direction,normalize,keepdims,answer',
                         [
                          ([1,0,0],'z',False,True, [1,0,0]),
                          ([1,0,0],'z',True, False,[1,0]),
                          ([0,1,1],'z',False,True, [0,0.5,0]),
                          ([0,1,1],'y',True, False,[0.41421356,0]),
                          ([1,1,0],'x',False,False,[0.5,0]),
                          ([1,1,1],'y',True, True, [0.3660254, 0,0.3660254]),
                         ])
def test_project_equal_angle(point,direction,normalize,keepdims,answer):
    assert np.allclose(util.project_equal_angle(np.array(point),direction=direction,
                                                normalize=normalize,keepdims=keepdims),answer)

@pytest.mark.parametrize('point,direction,normalize,keepdims,answer',
                         [
                          ([1,0,0],'z',False,True, [1,0,0]),
                          ([1,0,0],'z',True, False,[1,0]),
                          ([0,1,1],'z',False,True, [0,0.70710678,0]),
                          ([0,1,1],'y',True, False,[0.5411961,0]),
                          ([1,1,0],'x',False,False,[0.70710678,0]),
                          ([1,1,1],'y',True, True, [0.45970084,0,0.45970084]),
                         ])
def test_project_equal_area(point,direction,normalize,keepdims,answer):
    assert np.allclose(util.project_equal_area(np.array(point),direction=direction,
                                                normalize=normalize,keepdims=keepdims),answer)

@pytest.mark.parametrize('fro,to,mode,answer',
                         [
                          ((),(1,),'left',(1,)),
                          ((1,),(7,),'right',(1,)),
                          ((1,2),(1,1,2,2),'right',(1,1,2,1)),
                          ((1,2),(1,1,2,2),'left',(1,1,1,2)),
                          ((1,2,3),(1,1,2,3,4),'right',(1,1,2,3,1)),
                          ((10,2),(10,3,2,2,),'right',(10,1,2,1)),
                          ((10,2),(10,3,2,2,),'left',(10,1,1,2)),
                          ((2,2,3),(2,2,2,3,4),'left',(1,2,2,3,1)),
                          ((2,2,3),(2,2,2,3,4),'right',(2,2,1,3,1)),
                         ])
def test_shapeshifter(fro,to,mode,answer):
    assert util.shapeshifter(fro,to,mode) == answer

@pytest.mark.parametrize('fro,to,mode',
                         [
                          ((10,3,4),(10,3,2,2),'left'),
                          ((2,3),(10,3,2,2),'right'),
                         ])
def test_invalid_shapeshifter(fro,to,mode):
    with pytest.raises(ValueError):
        util.shapeshifter(fro,to,mode)

@pytest.mark.parametrize('a,b,ones,answer',
                         [
                          ((),(1,),True,(1,)),
                          ((1,),(),False,(1,)),
                          ((1,1),(7,),False,(1,7)),
                          ((1,),(7,),False,(7,)),
                          ((1,),(7,),True,(1,7)),
                          ((2,),(2,2),False,(2,2)),
                          ((1,3),(2,3),False,(2,3)),
                          ((1,1,2),(2,2),False,(1,2,2)),
                          ((1,1,2),(2,2),True,(1,1,2,2)),
                          ((1,2,3),(2,3,4),False,(1,2,3,4)),
                          ((1,2,3),(1,2,3),False,(1,2,3)),
                          ((2,3,1,1),(2,3),False,(2,3,2,3)),
                          ((2,3,1,1),(2,3),True,(2,3,1,1,2,3)),
                         ])
def test_shapeblender(a,b,ones,answer):
    assert util.shapeblender(a,b,ones) == answer

@pytest.mark.parametrize('style',[util.emph,util.deemph,util.warn,util.strikeout])
def test_decorate(style):
    assert 'DAMASK' in style('DAMASK')

@pytest.mark.parametrize('complete',[True,False])
@pytest.mark.parametrize('fhandle',[True,False])
def test_D3D_base_group(np_rng,tmp_path,complete,fhandle):
    random.seed(int(np_rng.integers(np.iinfo(int).max)))
    base_group = ''.join(random.choices('DAMASK', k=10))
    with h5py.File(tmp_path/'base_group.dream3d','w') as f:
        f.create_group('/'.join((base_group,'_SIMPL_GEOMETRY')))
        if complete:
            f['/'.join((base_group,'_SIMPL_GEOMETRY'))].create_dataset('SPACING',data=np.ones(3))

    fname = tmp_path/'base_group.dream3d'
    if fhandle: fname = h5py.File(fname)
    if complete:
        assert base_group == util.DREAM3D_base_group(fname)
    else:
        with pytest.raises(ValueError):
            util.DREAM3D_base_group(fname)

@pytest.mark.parametrize('complete',[True,False])
@pytest.mark.parametrize('fhandle',[True,False])
def test_D3D_cell_data_group(np_rng,tmp_path,complete,fhandle):
    random.seed(int(np_rng.integers(np.iinfo(int).max)))
    base_group = ''.join(random.choices('DAMASK', k=10))
    cell_data_group = ''.join(random.choices('KULeuven', k=10))
    cells = np_rng.integers(1,50,3)
    with h5py.File(tmp_path/'cell_data_group.dream3d','w') as f:
        f.create_group('/'.join((base_group,'_SIMPL_GEOMETRY')))
        f['/'.join((base_group,'_SIMPL_GEOMETRY'))].create_dataset('SPACING',data=np.ones(3))
        f['/'.join((base_group,'_SIMPL_GEOMETRY'))].create_dataset('DIMENSIONS',data=cells[::-1])
        f[base_group].create_group(cell_data_group)
        if complete:
            f['/'.join((base_group,cell_data_group))].create_dataset('data',shape=np.append(cells,1))

    fname = tmp_path/'cell_data_group.dream3d'
    if fhandle: fname = h5py.File(fname)
    if complete:
        assert cell_data_group == util.DREAM3D_cell_data_group(fname)
    else:
        with pytest.raises(ValueError):
            util.DREAM3D_cell_data_group(fname)


@pytest.mark.parametrize('full,reduced',[({},                           {}),
                                         ({'A':{}},                     {}),
                                         ({'A':{'B':{}}},               {}),
                                         ({'A':{'B':'C'}},)*2,
                                         ({'A':{'B':{},'C':'D'}},       {'A':{'C':'D'}})])
def test_prune(full,reduced):
    assert util.dict_prune(full) == reduced


@pytest.mark.parametrize('full,reduced',[({},                           {}),
                                         ({'A':{}},                     {}),
                                         ({'A':'F'},                    'F'),
                                         ({'A':{'B':{}}},               {}),
                                         ({'A':{'B':'C'}},              'C'),
                                         ({'A':1,'B':2},)*2,
                                         ({'A':{'B':'C','D':'E'}},      {'B':'C','D':'E'}),
                                         ({'B':'C','D':'E'},)*2,
                                         ({'A':{'B':{},'C':'D'}},       {'B':{},'C':'D'})])
def test_flatten(full,reduced):
    assert util.dict_flatten(full) == reduced


def test_double_Bravais_to_Miller():
    with pytest.raises(KeyError):
        util.Bravais_to_Miller(uvtw=np.ones(4),hkil=np.ones(4))

def test_double_Miller_to_Bravais():
    with pytest.raises(KeyError):
        util.Miller_to_Bravais(uvw=np.ones(4),hkl=np.ones(4))

@pytest.mark.parametrize('key_value',[{'uvtw':[1.,0.,-1.,1.1]},
                                      {'hkil':[1.,0.,-1.,1.1]}])
def test_float_Bravais_to_Miller(key_value):
    with pytest.raises(ValueError):
        util.Bravais_to_Miller(**key_value)

@pytest.mark.parametrize('key_value',[{'uvw':[1.,0.,-1.1]},
                                      {'hkl':[1.,0.,-9.4]}])
def test_float_Miller_to_Bravais(key_value):
    with pytest.raises(ValueError):
        util.Miller_to_Bravais(**key_value)


@pytest.mark.parametrize('key_value',[{'uvtw':[1,0,0,0]},
                                      {'hkil':[1,0,0,0]}])
def test_invalid_MillerBravais(key_value):
    with pytest.raises(ValueError):
        util.Bravais_to_Miller(**key_value)

@pytest.mark.parametrize('vector',np.array([
                                            [1,0,0],
                                            [1,1,0],
                                            [1,1,1],
                                            [1,0,-2],
                                           ]))
@pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
def test_Miller_Bravais_Miller(vector,kw_Miller,kw_Bravais):
    assert np.all(vector == util.Bravais_to_Miller(**{kw_Bravais:util.Miller_to_Bravais(**{kw_Miller:vector})}))

@pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
def test_Miller_Bravais_Miller_random(np_rng,kw_Miller,kw_Bravais):
    vector = np_rng.integers(-25,26,size=(5,6,3))
    vector //= np.gcd.reduce(vector,axis=-1,keepdims=True)
    assert np.all(vector == util.Bravais_to_Miller(**{kw_Bravais:util.Miller_to_Bravais(**{kw_Miller:vector})}))


@pytest.mark.parametrize('vector',np.array([
                                            [1,0,-1,2],
                                            [1,-1,0,3],
                                            [1,1,-2,-3],
                                            [0,0,0,1],
                                           ]))
@pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
def test_Bravais_Miller_Bravais(vector,kw_Miller,kw_Bravais):
    assert np.all(vector == util.Miller_to_Bravais(**{kw_Miller:util.Bravais_to_Miller(**{kw_Bravais:vector})}))

@pytest.mark.parametrize('dim',(None,1,4))
def test_standardize_MillerBravais(np_rng,dim):
    shape = tuple(np_rng.integers(1,6,dim))+(3,) if dim is not None else (3,)
    idx_red = np_rng.integers(-10,11,shape)
    idx_full = np.block([idx_red[...,:2], -np.sum(idx_red[...,:2],axis=-1,keepdims=True), idx_red[...,2:]])
    idx_missing = idx_full.astype(object)
    idx_missing[...,2][idx_full[...,3]>0] = ...
    assert np.equal(idx_full,util._standardize_MillerBravais(idx_red)).all() and \
            np.equal(idx_full,util._standardize_MillerBravais(idx_missing)).all()
