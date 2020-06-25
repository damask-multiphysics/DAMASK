from damask import util

class TestUtil:

    def test_execute_direct(self):
        out,err = util.execute('echo test')
        assert out=='test\n' and err==''
    
    def test_execute_env(self):
        out,err = util.execute('sh -c "echo $test_for_execute"',env={'test_for_execute':'test'})
        assert out=='test\n' and err==''
