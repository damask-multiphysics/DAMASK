import pytest
import numpy as np
from scipy import stats

from damask import util


class TestUtil:

    def test_execute_direct(self):
        out,err = util.execute('echo test')
        assert out=='test\n' and err==''

    def test_execute_env(self):
        out,err = util.execute('sh -c "echo $test_for_execute"',env={'test_for_execute':'test'})
        assert out=='test\n' and err==''

    def test_croak(self):
        util.croak('Burp!')

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

    @pytest.mark.parametrize('point,normalize,answer',
                             [
                              ([1,0,0],False,[1,0,0]),
                              ([1,0,0],True, [1,0,0]),
                              ([0,1,1],False,[0,0.5,0]),
                              ([0,1,1],True, [0,0.41421356,0]),
                              ([1,1,1],False,[0.5,0.5,0]),
                              ([1,1,1],True, [0.3660254, 0.3660254, 0]),
                             ])
    def test_project_stereographic(self,point,normalize,answer):
        assert np.allclose(util.project_stereographic(np.array(point),normalize=normalize),answer)

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
