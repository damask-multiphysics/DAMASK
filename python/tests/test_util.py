import random
import os

import pytest
import numpy as np
from scipy import stats
import h5py

from damask import util


class TestUtil:

    def test_execute_direct(self):
        out,err = util.execute('echo test')
        assert out=='test\n' and err==''

    def test_execute_env(self):
        out,err = util.execute('sh -c "echo $test_for_execute"',env={'test_for_execute':'test'})
        assert out=='test\n' and err==''

    def test_execute_invalid(self):
        with pytest.raises(RuntimeError):
            util.execute('/bin/false')

    @pytest.mark.parametrize('input,output',
                            [
                            ([0,-2],[0,-1]),
                            ([-0.5,0.5],[-1,1]),
                            ([1./2.,1./3.],[3,2]),
                            ([2./3.,1./2.,1./3.],[4,3,2]),
                            ])

    def test_scale2coprime(self,input,output):
        assert np.allclose(util.scale_to_coprime(np.array(input)),
                                                 np.array(output).astype(int))

    def test_lackofprecision(self):
        with pytest.raises(ValueError):
            util.scale_to_coprime(np.array([1/333.333,1,1]))


    @pytest.mark.parametrize('rv',[stats.rayleigh(),stats.weibull_min(1.2),stats.halfnorm(),stats.pareto(2.62)])
    def test_hybridIA(self,rv):
        bins = np.linspace(0,10,100000)
        centers = (bins[1:]+bins[:-1])/2
        N_samples = bins.shape[0]-1000
        dist = rv.pdf(centers)
        selected = util.hybrid_IA(dist,N_samples)
        dist_sampled = np.histogram(centers[selected],bins)[0]/N_samples*np.sum(dist)
        assert np.sqrt(((dist - dist_sampled) ** 2).mean()) < .025 and selected.shape[0]==N_samples

    @pytest.mark.parametrize('point,direction,normalize,keepdims,answer',
                             [
                              ([1,0,0],'z',False,True, [1,0,0]),
                              ([1,0,0],'z',True, False,[1,0]),
                              ([0,1,1],'z',False,True, [0,0.5,0]),
                              ([0,1,1],'y',True, False,[0.41421356,0]),
                              ([1,1,0],'x',False,False,[0.5,0]),
                              ([1,1,1],'y',True, True, [0.3660254, 0,0.3660254]),
                             ])
    def test_project_stereographic(self,point,direction,normalize,keepdims,answer):
        assert np.allclose(util.project_stereographic(np.array(point),direction=direction,
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
    def test_shapeshifter(self,fro,to,mode,answer):
        assert util.shapeshifter(fro,to,mode) == answer

    @pytest.mark.parametrize('fro,to,mode',
                             [
                              ((10,3,4),(10,3,2,2),'left'),
                              ((2,3),(10,3,2,2),'right'),
                             ])
    def test_invalid_shapeshifter(self,fro,to,mode):
        with pytest.raises(ValueError):
            util.shapeshifter(fro,to,mode)

    @pytest.mark.parametrize('a,b,answer',
                             [
                              ((),(1,),(1,)),
                              ((1,),(),(1,)),
                              ((1,),(7,),(1,7)),
                              ((2,),(2,2),(2,2)),
                              ((1,2),(2,2),(1,2,2)),
                              ((1,2,3),(2,3,4),(1,2,3,4)),
                              ((1,2,3),(1,2,3),(1,2,3)),
                             ])
    def test_shapeblender(self,a,b,answer):
        assert util.shapeblender(a,b) == answer

    @pytest.mark.parametrize('style',[util.emph,util.deemph,util.warn,util.strikeout])
    def test_decorate(self,style):
        assert 'DAMASK' in style('DAMASK')

    @pytest.mark.parametrize('complete',[True,False])
    def test_D3D_base_group(self,tmp_path,complete):
        base_group = ''.join(random.choices('DAMASK', k=10))
        with h5py.File(tmp_path/'base_group.dream3d','w') as f:
            f.create_group(os.path.join(base_group,'_SIMPL_GEOMETRY'))
            if complete:
                f[os.path.join(base_group,'_SIMPL_GEOMETRY')].create_dataset('SPACING',data=np.ones(3))

        if complete:
            assert base_group == util.DREAM3D_base_group(tmp_path/'base_group.dream3d')
        else:
            with pytest.raises(ValueError):
                util.DREAM3D_base_group(tmp_path/'base_group.dream3d')

    @pytest.mark.parametrize('complete',[True,False])
    def test_D3D_cell_data_group(self,tmp_path,complete):
        base_group = ''.join(random.choices('DAMASK', k=10))
        cell_data_group = ''.join(random.choices('KULeuven', k=10))
        cells = np.random.randint(1,50,3)
        with h5py.File(tmp_path/'cell_data_group.dream3d','w') as f:
            f.create_group(os.path.join(base_group,'_SIMPL_GEOMETRY'))
            f[os.path.join(base_group,'_SIMPL_GEOMETRY')].create_dataset('SPACING',data=np.ones(3))
            f[os.path.join(base_group,'_SIMPL_GEOMETRY')].create_dataset('DIMENSIONS',data=cells[::-1])
            f[base_group].create_group(cell_data_group)
            if complete:
                f[os.path.join(base_group,cell_data_group)].create_dataset('data',shape=np.append(cells,1))

        if complete:
            assert cell_data_group == util.DREAM3D_cell_data_group(tmp_path/'cell_data_group.dream3d')
        else:
            with pytest.raises(ValueError):
                util.DREAM3D_cell_data_group(tmp_path/'cell_data_group.dream3d')


    @pytest.mark.parametrize('full,reduced',[({},                           {}),
                                             ({'A':{}},                     {}),
                                             ({'A':{'B':{}}},               {}),
                                             ({'A':{'B':'C'}},)*2,
                                             ({'A':{'B':{},'C':'D'}},       {'A':{'C':'D'}})])
    def test_prune(self,full,reduced):
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
    def test_flatten(self,full,reduced):
        assert util.dict_flatten(full) == reduced
