import sys
import random
import pydoc

import pytest
import numpy as np
from scipy import stats
import h5py

from damask import util

class TestUtil:

    @pytest.mark.xfail(sys.platform == 'win32', reason='echo is not a Windows command')
    def test_run_direct(self):
        out,err = util.run('echo test')
        assert out=='test\n' and err==''

    @pytest.mark.xfail(sys.platform == 'win32', reason='echo is not a Windows command')
    def test_run_env(self):
        out,err = util.run('sh -c "echo $test_for_execute"',env={'test_for_execute':'test'})
        assert out=='test\n' and err==''

    @pytest.mark.xfail(sys.platform == 'win32', reason='false is not a Windows command')
    def test_run_runtime_error(self):
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
    def test_srepr(self,input,glue,quote,output):
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
    def test_scale2coprime(self,input,output,N):
        res = util.scale_to_coprime(input,N)
        assert np.allclose(res/np.max(np.abs(res)),output/np.max(np.abs(output)),atol=1e-2,rtol=0)


    @pytest.mark.parametrize('rv',[stats.rayleigh(),stats.weibull_min(1.2),stats.halfnorm(),stats.pareto(2.62)])
    def test_hybridIA_distribution(self,rv):
        bins = np.linspace(0,10,100000)
        centers = (bins[1:]+bins[:-1])/2
        N_samples = bins.shape[0]-1000
        dist = rv.pdf(centers)
        selected = util.hybrid_IA(dist,N_samples)
        dist_sampled = np.histogram(centers[selected],bins)[0]/N_samples*np.sum(dist)
        assert np.sqrt(((dist - dist_sampled) ** 2).mean()) < .025 and selected.shape[0]==N_samples

    def test_hybridIA_constant(self):
       N_bins = np.random.randint(20,400)
       m = np.random.randint(1,20)
       N_samples = m * N_bins
       dist = np.ones(N_bins)*np.random.rand()
       assert np.all(np.sort(util.hybrid_IA(dist,N_samples))==np.arange(N_samples).astype(int)//m)

    def test_hybridIA_linear(self):
       N_points = np.random.randint(10,200)
       m = np.random.randint(1,20)
       dist = np.arange(N_points)
       N_samples = m * np.sum(dist)
       assert np.all(np.bincount(util.hybrid_IA(dist*np.random.rand(),N_samples)) == dist*m)


    @pytest.mark.parametrize('point,direction,normalize,keepdims,answer',
                             [
                              ([1,0,0],'z',False,True, [1,0,0]),
                              ([1,0,0],'z',True, False,[1,0]),
                              ([0,1,1],'z',False,True, [0,0.5,0]),
                              ([0,1,1],'y',True, False,[0.41421356,0]),
                              ([1,1,0],'x',False,False,[0.5,0]),
                              ([1,1,1],'y',True, True, [0.3660254, 0,0.3660254]),
                             ])
    def test_project_equal_angle(self,point,direction,normalize,keepdims,answer):
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
    def test_project_equal_area(self,point,direction,normalize,keepdims,answer):
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
    def test_shapeblender(self,a,b,ones,answer):
        assert util.shapeblender(a,b,ones) == answer

    @pytest.mark.parametrize('style',[util.emph,util.deemph,util.warn,util.strikeout])
    def test_decorate(self,style):
        assert 'DAMASK' in style('DAMASK')

    @pytest.mark.parametrize('complete',[True,False])
    @pytest.mark.parametrize('fhandle',[True,False])
    def test_D3D_base_group(self,tmp_path,complete,fhandle):
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
    def test_D3D_cell_data_group(self,tmp_path,complete,fhandle):
        base_group = ''.join(random.choices('DAMASK', k=10))
        cell_data_group = ''.join(random.choices('KULeuven', k=10))
        cells = np.random.randint(1,50,3)
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


    def test_double_Bravais_to_Miller(self):
        with pytest.raises(KeyError):
            util.Bravais_to_Miller(uvtw=np.ones(4),hkil=np.ones(4))

    def test_double_Miller_to_Bravais(self):
        with pytest.raises(KeyError):
            util.Miller_to_Bravais(uvw=np.ones(4),hkl=np.ones(4))


    @pytest.mark.parametrize('vector',np.array([
                                                [1,0,0],
                                                [1,1,0],
                                                [1,1,1],
                                                [1,0,-2],
                                               ]))
    @pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
    def test_Miller_Bravais_Miller(self,vector,kw_Miller,kw_Bravais):
        assert np.all(vector == util.Bravais_to_Miller(**{kw_Bravais:util.Miller_to_Bravais(**{kw_Miller:vector})}))

    @pytest.mark.parametrize('vector',np.array([
                                                [1,0,-1,2],
                                                [1,-1,0,3],
                                                [1,1,-2,-3],
                                                [0,0,0,1],
                                               ]))
    @pytest.mark.parametrize('kw_Miller,kw_Bravais',[('uvw','uvtw'),('hkl','hkil')])
    def test_Bravais_Miller_Bravais(self,vector,kw_Miller,kw_Bravais):
        assert np.all(vector == util.Miller_to_Bravais(**{kw_Miller:util.Bravais_to_Miller(**{kw_Bravais:vector})}))

    @pytest.mark.parametrize('adopted_parameters',[
            pytest.param("""
            p2 : str, optional
                p2 description 1
                p2 description 2
            """,
            id = 'standard'),
            pytest.param("""

                        p2 : str, optional
                            p2 description 1
                            p2 description 2

            """,
            id = 'indented'),
            pytest.param("""
p2 : str, optional
    p2 description 1
    p2 description 2
            """,
            id = 'no_indent')])
    def test_extend_docstring_parameters_string(self,adopted_parameters):
        test_docstring = """
            Function description.

            Parameters
            ----------
            p0 : numpy.ndarray, shape (...,4)
                p0 description 1
                p0 description 2
            p1 : int, optional
                p1 description

            Remaining description\n"""
        res = util._docstringer(test_docstring,adopted_parameters)
        assert res ==\
        """
            Function description.

            Parameters
            ----------
            p0 : numpy.ndarray, shape (...,4)
                p0 description 1
                p0 description 2
            p1 : int, optional
                p1 description
            p2 : str, optional
                p2 description 1
                p2 description 2

            Remaining description\n"""

    def test_extend_docstring_parameters_function(self):
        test_docstring = """
            Function description.

            Parameters
            ----------
            p0 : numpy.ndarray, shape (...,4)
                p0 description 1

        """

        def testfunction_1():
            """
            Function description.

            Parameters
            ----------
            p1 : int, optional
                p1 description

            Notes
            -----
            Function Notes 1
            Function Notes 2

            References
            ----------
            Reference 1
            <reference link>
            Reference 2
            <reference link>

            """
            pass

        test_docstring = util._docstringer(test_docstring,adopted_references = testfunction_1)
        assert test_docstring == \
        """
            Function description.

            Parameters
            ----------
            p0 : numpy.ndarray, shape (...,4)
                p0 description 1

            References
            ----------
            Reference 1
            <reference link>
            Reference 2
            <reference link>\n"""
        test_docstring = util._docstringer(test_docstring,adopted_notes = testfunction_1)
        assert test_docstring == \
        """
            Function description.

            Parameters
            ----------
            p0 : numpy.ndarray, shape (...,4)
                p0 description 1

            Notes
            -----
            Function Notes 1
            Function Notes 2

            References
            ----------
            Reference 1
            <reference link>
            Reference 2
            <reference link>\n"""

        def testfunction_2():
            """
            Function description.

            References
            ----------
            Reference 3
            <reference link>

            """
        test_docstring = util._docstringer(test_docstring,adopted_references = testfunction_2)
        assert test_docstring == \
        """
            Function description.

            Parameters
            ----------
            p0 : numpy.ndarray, shape (...,4)
                p0 description 1

            Notes
            -----
            Function Notes 1
            Function Notes 2

            References
            ----------
            Reference 1
            <reference link>
            Reference 2
            <reference link>
            Reference 3
            <reference link>\n"""

    def return_bound_method():
        class TestClassDecorated:
            def decorated_func_bound(self) -> 'TestClassDecorated':
                pass
        return TestClassDecorated.decorated_func_bound

    def return_simple_function():
        class TestClassDecorated:
            pass

        def decorated_func() -> TestClassDecorated:
            pass
        return decorated_func

    @pytest.mark.parametrize('adopted_return',[
            pytest.param(return_simple_function(),
            id = 'decorated_func'),
            pytest.param(return_bound_method(),
            id = 'decorated_func_bound'),
            pytest.param('test_util.TestClassDecorated',
            id = 'decorated_func_bound')])
    def test_replace_docstring_return(self,adopted_return):
        class TestClassOriginal:
            pass

        def original_func() -> TestClassOriginal:
            pass

        original_func.__doc__ = """
            Function description.

            Returns
            -------
            Return value : test_util.TestClassOriginal

            Remaining description
            """


        assert util._docstringer(original_func,adopted_return=adopted_return) == """
            Function description.

            Returns
            -------
            Return value : test_util.TestClassDecorated

            Remaining description\n"""


    @pytest.mark.parametrize('adopted_func_doc',[
            pytest.param("""
            Function description.

            Parameters
            ----------
            b : float 4
                b description
            c : float 5
                c description differing
            d : float 6
                d description

            Remaining description\n
            """,
            id = 'append'),
            pytest.param("""
            Function description.

            Parameters
            ----------
            d : float 6
                d description
            a : float 7
                a description\n""",
            id = 'insert')])
    def test_extend_docstring_overlapping_section_content(self,adopted_func_doc):

        original_func_doc = """
            Function description.

            Parameters
            ----------
            a : float 1
                a description
            b : float 2
                b description
            c : float 3
                c description

            Remaining description\n"""

        expected = """
            Function description.

            Parameters
            ----------
            a : float 1
                a description
            b : float 2
                b description
            c : float 3
                c description
            d : float 6
                d description

            Remaining description\n"""

        assert util._docstringer(original_func_doc,adopted_parameters=adopted_func_doc) == expected

    def test_passon_result(self):
        def testfunction_inner(a=None,b=None):
            return a+b

        @util.pass_on('inner_result',testfunction_inner)
        def testfunction_outer(**kwargs):
            return kwargs['inner_result']+";"+kwargs['c']+kwargs['d']
        assert testfunction_outer(a='1',b='2',c='3',d='4',e='5') == '12;34'

    def test_passon_signature(self):
        def testfunction_inner(a='1',b='2'):
            return a+b

        def testfunction_extra(e='5',f='6'):
            return e+f

        @util.pass_on('inner_result', testfunction_inner, wrapped=testfunction_extra)
        def testfunction_outer(**kwargs):
            return kwargs['inner_result']+";"+kwargs['c']+kwargs['d']
        assert [(param.name, param.default) for param in testfunction_outer.__signature__.parameters.values()] == \
               [('a', '1'), ('b', '2'), ('e', '5'), ('f', '6')]

    def test_passon_help(self):
        def testfunction_inner(a=None,b=None):
            return a+b

        def testfunction_extra(*,c=None,d=None):
            return c+d

        @util.pass_on('inner_result', testfunction_inner, wrapped=testfunction_extra)
        def testfunction_outer(**kwargs) -> int:
            return kwargs['inner_result']+kwargs['c']+kwargs['d']

        assert pydoc.render_doc(testfunction_outer, renderer=pydoc.plaintext).split("\n")[-2] ==\
              'testfunction_outer(*, a=None, b=None, c=None, d=None) -> int'
