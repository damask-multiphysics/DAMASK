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
