import bz2
import pickle
import time
import shutil
import os
import sys
import hashlib
from datetime import datetime

import pytest
import numpy as np

from damask import Result
from damask import Rotation
from damask import Orientation
from damask import tensor
from damask import mechanics
from damask import grid_filters

@pytest.fixture
def default(tmp_path,ref_path):
    """Small Result file in temp location for modification."""
    fname = '12grains6x7x8_tensionY.hdf5'
    shutil.copy(ref_path/fname,tmp_path)
    f = Result(tmp_path/fname)
    return f.view('times',20.0)

@pytest.fixture
def single_phase(tmp_path,ref_path):
    """Single phase Result file in temp location for modification."""
    fname = '6grains6x7x8_single_phase_tensionY.hdf5'
    shutil.copy(ref_path/fname,tmp_path)
    return Result(tmp_path/fname)

@pytest.fixture
def ref_path(ref_path_base):
    """Directory containing reference results."""
    return ref_path_base/'Result'

def dict_equal(d1, d2):
    for k in d1:
        if (k not in d2):
            return False
        else:
            if type(d1[k]) is dict:
                return dict_equal(d1[k],d2[k])
            else:
                if not np.allclose(d1[k],d2[k]):
                    return False
    return True

class TestResult:

    def test_self_report(self,default):
        print(default)


    def test_view_all(self,default):
        a = default.view('increments',True).get('F')

        assert dict_equal(a,default.view('increments','*').get('F'))
        assert dict_equal(a,default.view('increments',default.increments_in_range(0,np.iinfo(int).max)).get('F'))

        assert dict_equal(a,default.view('times',True).get('F'))
        assert dict_equal(a,default.view('times','*').get('F'))
        assert dict_equal(a,default.view('times',default.times_in_range(0.0,np.inf)).get('F'))

    @pytest.mark.parametrize('what',['increments','times','phases','fields'])                       # ToDo: discuss homogenizations
    def test_view_none(self,default,what):
        n0 = default.view(what,False)
        n1 = default.view(what,[])

        label = 'increments' if what == 'times' else what

        assert n0.get('F') is n1.get('F') is None and \
               len(n0.visible[label]) == len(n1.visible[label]) == 0

    @pytest.mark.parametrize('what',['increments','times','phases','fields'])                       # ToDo: discuss homogenizations
    def test_view_more(self,default,what):
        empty = default.view(what,False)

        a = empty.view_more(what,'*').get('F')
        b = empty.view_more(what,True).get('F')

        assert dict_equal(a,b)

    @pytest.mark.parametrize('what',['increments','times','phases','fields'])                       # ToDo: discuss homogenizations
    def test_view_less(self,default,what):
        full = default.view(what,True)

        n0 = full.view_less(what,'*')
        n1 = full.view_less(what,True)

        label = 'increments' if what == 'times' else what

        assert n0.get('F') is n1.get('F') is None and \
               len(n0.visible[label]) == len(n1.visible[label]) == 0

    def test_view_invalid(self,default):
        with pytest.raises(AttributeError):
            default.view('invalid',True)

    def test_add_absolute(self,default):
        default.add_absolute('F_e')
        in_memory = np.abs(default.place('F_e'))
        in_file   = default.place('|F_e|')
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('mode',['direct','function'])
    def test_add_calculation(self,default,tmp_path,mode):

        if mode == 'direct':
            default.add_calculation('2.0*np.abs(#F#)-1.0','x','-','my notes')
        else:
            with open(tmp_path/'f.py','w') as f:
                f.write("import numpy as np\ndef my_func(field):\n  return 2.0*np.abs(field)-1.0\n")
            sys.path.insert(0,str(tmp_path))
            import f
            default.enable_user_function(f.my_func)
            default.add_calculation('my_func(#F#)','x','-','my notes')

        in_memory = 2.0*np.abs(default.place('F'))-1.0
        in_file   = default.place('x')
        assert np.allclose(in_memory,in_file)

    def test_add_stress_Cauchy(self,default):
        default.add_stress_Cauchy('P','F')
        in_memory = mechanics.stress_Cauchy(default.place('P'), default.place('F'))
        in_file   = default.place('sigma')
        assert np.allclose(in_memory,in_file)

    def test_add_determinant(self,default):
        default.add_determinant('P')
        in_memory = np.linalg.det(default.place('P'))
        in_file   = default.place('det(P)')
        assert np.allclose(in_memory,in_file)

    def test_add_deviator(self,default):
        default.add_deviator('P')
        in_memory = tensor.deviatoric(default.place('P'))
        in_file   = default.place('s_P')
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('eigenvalue,function',[('max',np.amax),('min',np.amin)])
    def test_add_eigenvalue(self,default,eigenvalue,function):
        default.add_stress_Cauchy('P','F')
        default.add_eigenvalue('sigma',eigenvalue)
        in_memory = function(tensor.eigenvalues(default.place('sigma')),axis=1)
        in_file   = default.place(f'lambda_{eigenvalue}(sigma)')
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('eigenvalue,idx',[('max',2),('mid',1),('min',0)])
    def test_add_eigenvector(self,default,eigenvalue,idx):
        default.add_stress_Cauchy('P','F')
        default.add_eigenvector('sigma',eigenvalue)
        in_memory = tensor.eigenvectors(default.place('sigma'))[:,idx]
        in_file   = default.place(f'v_{eigenvalue}(sigma)')
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('d',[[1,0,0],[0,1,0],[0,0,1]])
    def test_add_IPF_color(self,default,d):
        default.add_IPF_color(d,'O')
        qu = default.place('O')
        crystal_structure = qu.dtype.metadata['lattice']
        c = Orientation(rotation=qu,lattice=crystal_structure)
        in_memory = np.uint8(c.IPF_color(np.array(d))*255)
        in_file = default.place('IPFcolor_({} {} {})'.format(*d))
        assert np.allclose(in_memory,in_file)

    def test_add_maximum_shear(self,default):
        default.add_stress_Cauchy('P','F')
        default.add_maximum_shear('sigma')
        in_memory = mechanics.maximum_shear(default.place('sigma'))
        in_file   = default.place('max_shear(sigma)')
        assert np.allclose(in_memory,in_file)

    def test_add_Mises_strain(self,default):
        t = ['V','U'][np.random.randint(0,2)]
        m = np.random.random()*2.0 - 1.0
        default.add_strain('F',t,m)
        label = f'epsilon_{t}^{m}(F)'
        default.add_equivalent_Mises(label)
        in_memory = mechanics.equivalent_strain_Mises(default.place(label))
        in_file   = default.place(label+'_vM')
        assert np.allclose(in_memory,in_file)

    def test_add_Mises_stress(self,default):
        default.add_stress_Cauchy('P','F')
        default.add_equivalent_Mises('sigma')
        in_memory = mechanics.equivalent_stress_Mises(default.place('sigma'))
        in_file   = default.place('sigma_vM')
        assert np.allclose(in_memory,in_file)

    def test_add_Mises_invalid(self,default):
        default.add_stress_Cauchy('P','F')
        default.add_calculation('#sigma#','sigma_y',unit='y')
        default.add_equivalent_Mises('sigma_y')
        assert default.get('sigma_y_vM') is None

    def test_add_Mises_stress_strain(self,default):
        default.add_stress_Cauchy('P','F')
        default.add_calculation('#sigma#','sigma_y',unit='y')
        default.add_calculation('#sigma#','sigma_x',unit='x')
        default.add_equivalent_Mises('sigma_y',kind='strain')
        default.add_equivalent_Mises('sigma_x',kind='stress')
        assert not np.allclose(default.place('sigma_y_vM'),default.place('sigma_x_vM'))

    @pytest.mark.parametrize('ord',[1,2])
    @pytest.mark.parametrize('dataset,axis',[('F',(1,2)),('xi_sl',(1,))])
    def test_add_norm(self,default,ord,dataset,axis):
        default.add_norm(dataset,ord)
        in_memory = np.linalg.norm(default.place(dataset),ord=ord,axis=axis,keepdims=True)
        in_file   = default.place(f'|{dataset}|_{ord}')
        assert np.allclose(in_memory,in_file)

    def test_add_stress_second_Piola_Kirchhoff(self,default):
        default.add_stress_second_Piola_Kirchhoff('P','F')
        in_memory = mechanics.stress_second_Piola_Kirchhoff(default.place('P'),default.place('F'))
        in_file   = default.place('S')
        assert np.allclose(in_memory,in_file)

    @pytest.mark.skip(reason='requires rework of lattice.f90')
    @pytest.mark.parametrize('polar',[True,False])
    def test_add_pole(self,default,polar):
        pole = np.array([1.,0.,0.])
        default.add_pole('O',pole,polar)
        rot = Rotation(default.place('O'))
        rotated_pole = rot * np.broadcast_to(pole,rot.shape+(3,))
        xy = rotated_pole[:,0:2]/(1.+abs(pole[2]))
        in_memory = xy if not polar else \
                    np.block([np.sqrt(xy[:,0:1]*xy[:,0:1]+xy[:,1:2]*xy[:,1:2]),np.arctan2(xy[:,1:2],xy[:,0:1])])
        in_file = default.place('p^{}_[1 0 0)'.format(u'rφ' if polar else 'xy'))
        assert np.allclose(in_memory,in_file)

    def test_add_rotation(self,default):
        default.add_rotation('F')
        in_memory = mechanics.rotation(default.place('F')).as_matrix()
        in_file   = default.place('R(F)')
        assert np.allclose(in_memory,in_file)

    def test_add_spherical(self,default):
        default.add_spherical('P')
        in_memory = tensor.spherical(default.place('P'),False)
        in_file   = default.place('p_P')
        assert np.allclose(in_memory,in_file)

    def test_add_strain(self,default):
        t = ['V','U'][np.random.randint(0,2)]
        m = np.random.random()*2.0 - 1.0
        default.add_strain('F',t,m)
        label = f'epsilon_{t}^{m}(F)'
        in_memory = mechanics.strain(default.place('F'),t,m)
        in_file   = default.place(label)
        assert np.allclose(in_memory,in_file)

    def test_add_stretch_right(self,default):
        default.add_stretch_tensor('F','U')
        in_memory = mechanics.stretch_right(default.place('F'))
        in_file   = default.place('U(F)')
        assert np.allclose(in_memory,in_file)

    def test_add_stretch_left(self,default):
        default.add_stretch_tensor('F','V')
        in_memory = mechanics.stretch_left(default.place('F'))
        in_file   = default.place('V(F)')
        assert np.allclose(in_memory,in_file)

    def test_add_invalid(self,default):
        with pytest.raises(TypeError):
            default.add_calculation('#invalid#*2')


    @pytest.mark.parametrize('shape',['vector','tensor'])
    def test_add_curl(self,default,shape):
        if shape == 'vector': default.add_calculation('#F#[:,:,0]','x','1','just a vector')
        if shape == 'tensor': default.add_calculation('#F#[:,:,:]','x','1','just a tensor')
        x = default.place('x')
        default.add_curl('x')
        in_file   = default.place('curl(x)')
        in_memory = grid_filters.curl(default.size,x.reshape(tuple(default.cells)+x.shape[1:])).reshape(in_file.shape)
        assert (in_file==in_memory).all()

    @pytest.mark.parametrize('shape',['vector','tensor'])
    def test_add_divergence(self,default,shape):
        if shape == 'vector': default.add_calculation('#F#[:,:,0]','x','1','just a vector')
        if shape == 'tensor': default.add_calculation('#F#[:,:,:]','x','1','just a tensor')
        x = default.place('x')
        default.add_divergence('x')
        in_file   = default.place('divergence(x)')
        in_memory = grid_filters.divergence(default.size,x.reshape(tuple(default.cells)+x.shape[1:])).reshape(in_file.shape)
        assert (in_file==in_memory).all()

    @pytest.mark.parametrize('shape',['scalar','pseudo_scalar','vector'])
    def test_add_gradient(self,default,shape):
        if shape == 'pseudo_scalar': default.add_calculation('#F#[:,0,0:1]','x','1','a pseudo scalar')
        if shape == 'scalar': default.add_calculation('#F#[:,0,0]','x','1','just a scalar')
        if shape == 'vector': default.add_calculation('#F#[:,:,1]','x','1','just a vector')
        x = default.place('x').reshape((np.product(default.cells),-1))
        default.add_gradient('x')
        in_file   = default.place('gradient(x)')
        in_memory = grid_filters.gradient(default.size,x.reshape(tuple(default.cells)+x.shape[1:])).reshape(in_file.shape)
        assert (in_file==in_memory).all()

    @pytest.mark.parametrize('overwrite',['off','on'])
    def test_add_overwrite(self,default,overwrite):
        last = default.view('increments',-1)

        last.add_stress_Cauchy()

        created_first = last.place('sigma').dtype.metadata['created']
        created_first = datetime.strptime(created_first,'%Y-%m-%d %H:%M:%S%z')

        if overwrite == 'on':
            last = last.modification_enable()
        else:
            last = last.modification_disable()

        time.sleep(2.)
        try:
            last.add_calculation('#sigma#*0.0+311.','sigma','not the Cauchy stress')
        except ValueError:
            pass

        created_second = last.place('sigma').dtype.metadata['created']
        created_second = datetime.strptime(created_second,'%Y-%m-%d %H:%M:%S%z')

        if overwrite == 'on':
            assert created_first < created_second and np.allclose(last.place('sigma'),311.)
        else:
            assert created_first == created_second and not np.allclose(last.place('sigma'),311.)

    @pytest.mark.parametrize('allowed',['off','on'])
    def test_rename(self,default,allowed):
        if allowed == 'on':
            F = default.place('F')
            default = default.modification_enable()
            default.rename('F','new_name')
            assert np.all(F == default.place('new_name'))
            default = default.modification_disable()

        with pytest.raises(PermissionError):
            default.rename('P','another_new_name')

    @pytest.mark.parametrize('allowed',['off','on'])
    def test_remove(self,default,allowed):
        if allowed == 'on':
            unsafe = default.modification_enable()
            unsafe.remove('F')
            assert unsafe.get('F') is None
        else:
            with pytest.raises(PermissionError):
                default.remove('F')

    @pytest.mark.parametrize('mode',['cell','node'])
    def test_coordinates(self,default,mode):
         if   mode == 'cell':
             a = grid_filters.coordinates0_point(default.cells,default.size,default.origin)
             b = default.coordinates0_point.reshape(tuple(default.cells)+(3,),order='F')
         elif mode == 'node':
             a = grid_filters.coordinates0_node(default.cells,default.size,default.origin)
             b = default.coordinates0_node.reshape(tuple(default.cells+1)+(3,),order='F')
         assert np.allclose(a,b)

    # need to wait for writing in parallel, output order might change if select more then one
    @pytest.mark.parametrize('output',['F','*',['P']],ids=range(3))
    @pytest.mark.parametrize('fname',['12grains6x7x8_tensionY.hdf5'],ids=range(1))
    @pytest.mark.parametrize('inc',[4,0],ids=range(2))
    def test_vtk(self,request,tmp_path,ref_path,update,patch_execution_stamp,patch_datetime_now,output,fname,inc):
        result = Result(ref_path/fname).view('increments',inc)
        os.chdir(tmp_path)
        result.save_VTK(output)
        fname = fname.split('.')[0]+f'_inc{(inc if type(inc) == int else inc[0]):0>2}.vtr'
        last = ''
        for i in range(10):
            if os.path.isfile(tmp_path/fname):
                with open(fname) as f:
                    cur = hashlib.md5(f.read().encode()).hexdigest()
                    if cur == last:
                        break
                    else:
                        last = cur
            time.sleep(.5)
        if update:
            with open((ref_path/'save_VTK'/request.node.name).with_suffix('.md5'),'w') as f:
                f.write(cur)
        with open((ref_path/'save_VTK'/request.node.name).with_suffix('.md5')) as f:
            assert cur == f.read()

    @pytest.mark.parametrize('mode',['point','cell'])
    @pytest.mark.parametrize('output',[False,True])
    def test_vtk_marc(self,tmp_path,ref_path,mode,output):
        os.chdir(tmp_path)
        result = Result(ref_path/'check_compile_job1.hdf5')
        result.save_VTK(output,mode)

    def test_marc_coordinates(self,ref_path):
        result = Result(ref_path/'check_compile_job1.hdf5').view('increments',-1)
        c_n = result.coordinates0_node + result.get('u_n')
        c_p = result.coordinates0_point + result.get('u_p')
        assert len(c_n) > len(c_p)

    @pytest.mark.parametrize('mode',['point','cell'])
    def test_vtk_mode(self,tmp_path,single_phase,mode):
        os.chdir(tmp_path)
        single_phase.save_VTK(mode=mode)

    def test_XDMF(self,tmp_path,single_phase,update,ref_path):
        for shape in [('scalar',()),('vector',(3,)),('tensor',(3,3)),('matrix',(12,))]:
            for dtype in ['f4','f8','i1','i2','i4','i8','u1','u2','u4','u8']:
                 single_phase.add_calculation(f"np.ones(np.shape(#F#)[0:1]+{shape[1]},'{dtype}')",f'{shape[0]}_{dtype}')
        fname = os.path.splitext(os.path.basename(single_phase.fname))[0]+'.xdmf'
        os.chdir(tmp_path)
        single_phase.save_XDMF()
        if update:
            shutil.copy(tmp_path/fname,ref_path/fname)
        assert sorted(open(tmp_path/fname).read()) == sorted(open(ref_path/fname).read())           # XML is not ordered

    def test_XDMF_invalid(self,default):
        with pytest.raises(TypeError):
            default.save_XDMF()

    @pytest.mark.parametrize('view,output,flatten,prune',
            [({},['F','P','F','L_p','F_e','F_p'],True,True),
             ({'increments':3},'F',True,True),
             ({'increments':[1,8,3,4,5,6,7]},['F','P'],True,True),
             ({'phases':['A','B']},['F','P'],True,True),
             ({'phases':['A','C'],'homogenizations':False},['F','P','O'],True,True),
             ({'phases':False,'homogenizations':False},['F','P','O'],True,True),
             ({'phases':False},['Delta_V'],True,True),
             ({},['u_p','u_n'],False,False)],
            ids=list(range(8)))
    def test_get(self,update,request,ref_path,view,output,flatten,prune):
        result = Result(ref_path/'4grains2x4x3_compressionY.hdf5')
        for key,value in view.items():
            result = result.view(key,value)

        fname = request.node.name
        cur = result.get(output,flatten,prune)
        if update:
            with bz2.BZ2File((ref_path/'get'/fname).with_suffix('.pbz2'),'w') as f:
                pickle.dump(cur,f)

        with bz2.BZ2File((ref_path/'get'/fname).with_suffix('.pbz2')) as f:
            ref = pickle.load(f)
            assert cur is None if ref is None else dict_equal(cur,ref)

    @pytest.mark.parametrize('view,output,flatten,constituents,prune',
            [({},['F','P','F','L_p','F_e','F_p'],True,True,None),
             ({'increments':3},'F',True,True,[0,1,2,3,4,5,6,7]),
             ({'increments':[1,8,3,4,5,6,7]},['F','P'],True,True,1),
             ({'phases':['A','B']},['F','P'],True,True,[1,2]),
             ({'phases':['A','C'],'homogenizations':False},['F','P','O'],True,True,[0,7]),
             ({'phases':False,'homogenizations':False},['F','P','O'],True,True,[1,2,3,4]),
             ({'phases':False},['Delta_V'],True,True,[1,2,4]),
             ({},['u_p','u_n'],False,False,None)],
            ids=list(range(8)))
    def test_place(self,update,request,ref_path,view,output,flatten,prune,constituents):
        result = Result(ref_path/'4grains2x4x3_compressionY.hdf5')
        for key,value in view.items():
            result = result.view(key,value)

        fname = request.node.name
        cur = result.place(output,flatten,prune,constituents)
        if update:
            with bz2.BZ2File((ref_path/'place'/fname).with_suffix('.pbz2'),'w') as f:
                pickle.dump(cur,f)

        with bz2.BZ2File((ref_path/'place'/fname).with_suffix('.pbz2')) as f:
            ref = pickle.load(f)
            assert cur is None if ref is None else dict_equal(cur,ref)
