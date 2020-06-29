import pytest
import numpy as np
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
                            ([2,0],[1,0]),
                            ([0.5,0.5],[1,1]),
                            ([1./2.,1./3.],[3,2]),
                            ([2./3.,1./2.,1./3.],[4,3,2]),
                            ])

    def test_scale2coprime(self,input,output):
        assert np.allclose(util.scale_to_coprime(np.array(input)),
                                                 np.array(output).astype(int))

    def test_lackofprecision(self):
        with pytest.raises(ValueError):
            util.scale_to_coprime(np.array([1/3333,1,1]))
