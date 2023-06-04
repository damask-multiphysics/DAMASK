import bz2
import pickle
import time
import shutil
import os
import sys
import hashlib
import fnmatch
import random
from datetime import datetime

import pytest
from vtkmodules.vtkIOXML import vtkXMLImageDataReader
from vtkmodules.vtkCommonCore import vtkVersion
try:
    from vtkmodules.vtkIOXdmf2 import vtkXdmfReader
except ImportError:
    vtkXdmfReader=None                                                                              # noqa type: ignore
import h5py
import numpy as np

from damask import Result
from damask import Orientation
from damask import VTK
from damask import tensor
from damask import mechanics
from damask import grid_filters


@pytest.fixture
def default(tmp_path,res_path):
    """Small Result file in temp location for modification."""
    fname = '12grains6x7x8_tensionY.hdf5'
    shutil.copy(res_path/fname,tmp_path)
    return Result(tmp_path/fname).view(times=20.0)

@pytest.fixture
def single_phase(tmp_path,res_path):
    """Single phase Result file in temp location for modification."""
    fname = '6grains6x7x8_single_phase_tensionY.hdf5'
    shutil.copy(res_path/fname,tmp_path)
    return Result(tmp_path/fname)

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Result'

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
        a = default.view(increments=True).get('F')

        assert dict_equal(a,default.view(increments='*').get('F'))
        assert dict_equal(a,default.view(increments=default.increments_in_range(0,np.iinfo(int).max)).get('F'))

        assert dict_equal(a,default.view(times=True).get('F'))
        assert dict_equal(a,default.view(times='*').get('F'))
        assert dict_equal(a,default.view(times=default.times_in_range(0.0,np.inf)).get('F'))

    @pytest.mark.parametrize('what',['increments','times','phases','fields'])                       # ToDo: discuss homogenizations
    def test_view_none(self,default,what):
        n0 = default.view(**{what:False})
        n1 = default.view(**{what:[]})

        label = 'increments' if what == 'times' else what

        assert n0.get('F') is n1.get('F') is None and \
               len(n0.visible[label]) == len(n1.visible[label]) == 0

    @pytest.mark.parametrize('what',['increments','times','phases','fields'])                       # ToDo: discuss homogenizations
    def test_view_more(self,default,what):
        empty = default.view(**{what:False})

        a = empty.view_more(**{what:'*'}).get('F')
        b = empty.view_more(**{what:True}).get('F')

        assert dict_equal(a,b)

    @pytest.mark.parametrize('what',['increments','times','phases','fields'])                       # ToDo: discuss homogenizations
    def test_view_less(self,default,what):
        full = default.view(**{what:True})

        n0 = full.view_less(**{what:'*'})
        n1 = full.view_less(**{what:True})

        label = 'increments' if what == 'times' else what

        assert n0.get('F') is n1.get('F') is None and \
               len(n0.visible[label]) == len(n1.visible[label]) == 0

    def test_view_invalid_incstimes(self,default):
        with pytest.raises(ValueError):
            default.view(increments=0,times=0)

    @pytest.mark.parametrize('inc',[0,10])
    @pytest.mark.parametrize('sign',[+1,-1])
    def test_view_approxtimes(self,default,inc,sign):
        eps = sign*1e-3
        assert [default.increments[inc]] == default.view(times=default.times[inc]+eps).visible['increments']

    def test_add_invalid(self,default):
        default.add_absolute('xxxx')

    def test_add_absolute(self,default):
        default.add_absolute('F_e')
        in_memory = np.abs(default.place('F_e'))
        in_file   = default.place('|F_e|')
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('mode',
        ['direct',pytest.param('function',marks=pytest.mark.xfail(sys.platform in ['darwin','win32'], reason='n/a'))])
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

    def test_add_calculation_invalid(self,default):
        default.add_calculation('np.linalg.norm(#F#,axis=0)','wrong_dim')
        assert default.get('wrong_dim') is None

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

    @pytest.mark.parametrize('options',[{'uvw':[1,0,0],'with_symmetry':False},
                                        {'uvw':[1,1,0],'with_symmetry':True},
                                        {'hkl':[0,1,1],'with_symmetry':True},
                                        {'hkl':[1,1,1],'with_symmetry':False},
                                        ])
    def test_add_pole(self,default,options):
        default.add_pole(**options)
        rot = default.place('O')
        in_memory = Orientation(rot,lattice=rot.dtype.metadata['lattice']).to_pole(**options)
        brackets = [['[[]','[]]'],'()','⟨⟩','{}'][('hkl' in options)*1+(options['with_symmetry'])*2]   # escape fnmatch
        label = 'p^{}{} {} {}{}'.format(brackets[0],
                                        *(list(options.values())[0]),
                                        brackets[-1])
        in_file = default.place(label)
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

    def test_add_invalid_dataset(self,default):
        with pytest.raises(TypeError):
            default.add_calculation('#invalid#*2')

    def test_add_generic_grid_invalid(self,res_path):
        result = Result(res_path/'4grains2x4x3_compressionY.hdf5')
        with pytest.raises(NotImplementedError):
            result.add_curl('F')


    @pytest.mark.parametrize('shape',['vector','tensor'])
    def test_add_curl(self,default,shape):
        if shape == 'vector': default.add_calculation('#F#[:,:,0]','x','1','just a vector')
        if shape == 'tensor': default.add_calculation('#F#[:,:,:]','x','1','just a tensor')
        x = default.place('x')
        default.add_curl('x')
        in_file   = default.place('curl(x)')
        in_memory = grid_filters.curl(default.size,x.reshape(tuple(default.cells)+x.shape[1:])).reshape(in_file.shape)
        assert (in_file == in_memory).all()

    @pytest.mark.parametrize('shape',['vector','tensor'])
    def test_add_divergence(self,default,shape):
        if shape == 'vector': default.add_calculation('#F#[:,:,0]','x','1','just a vector')
        if shape == 'tensor': default.add_calculation('#F#[:,:,:]','x','1','just a tensor')
        x = default.place('x')
        default.add_divergence('x')
        in_file   = default.place('divergence(x)')
        in_memory = grid_filters.divergence(default.size,x.reshape(tuple(default.cells)+x.shape[1:])).reshape(in_file.shape)
        assert (in_file == in_memory).all()

    @pytest.mark.parametrize('shape',['scalar','pseudo_scalar','vector'])
    def test_add_gradient(self,default,shape):
        if shape == 'pseudo_scalar': default.add_calculation('#F#[:,0,0:1]','x','1','a pseudo scalar')
        if shape == 'scalar': default.add_calculation('#F#[:,0,0]','x','1','just a scalar')
        if shape == 'vector': default.add_calculation('#F#[:,:,1]','x','1','just a vector')
        x = default.place('x').reshape((np.product(default.cells),-1))
        default.add_gradient('x')
        in_file   = default.place('gradient(x)')
        in_memory = grid_filters.gradient(default.size,x.reshape(tuple(default.cells)+x.shape[1:])).reshape(in_file.shape)
        assert (in_file == in_memory).all()

    @pytest.mark.parametrize('overwrite',['off','on'])
    def test_add_overwrite(self,default,overwrite):
        last = default.view(increments=-1)

        last.add_stress_Cauchy()

        created_first = last.place('sigma').dtype.metadata['created']
        created_first = datetime.strptime(created_first,'%Y-%m-%d %H:%M:%S%z')

        last = last.view(protected=overwrite != 'on')

        time.sleep(2)
        try:
            last.add_calculation('#sigma#*0.0+311.','sigma','not the Cauchy stress')
        except ValueError:
            pass

        created_second = last.place('sigma').dtype.metadata['created']
        created_second = datetime.strptime(created_second,'%Y-%m-%d %H:%M:%S%z')

        if overwrite == 'on':
            assert created_first  < created_second and     np.allclose(last.place('sigma'),311.)
        else:
            assert created_first == created_second and not np.allclose(last.place('sigma'),311.)

    @pytest.mark.parametrize('allowed',['off','on'])
    def test_rename(self,default,allowed):
        if allowed == 'on':
            F = default.place('F')
            default = default.view(protected=False)
            default.rename('F','new_name')
            assert np.all(F == default.place('new_name'))
            default = default.view(protected=True)

        with pytest.raises(PermissionError):
            default.rename('P','another_new_name')

    @pytest.mark.parametrize('allowed',['off','on'])
    def test_remove(self,default,allowed):
        if allowed == 'on':
            unsafe = default.view(protected=False)
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

    @pytest.mark.parametrize('output',['F','*',['P'],['P','F']],ids=range(4))
    @pytest.mark.parametrize('fname',['12grains6x7x8_tensionY.hdf5',
                                      '4grains2x4x3_compressionY.hdf5',
                                      '6grains6x7x8_single_phase_tensionY.hdf5'],ids=range(3))
    @pytest.mark.parametrize('inc',[4,0],ids=range(2))
    @pytest.mark.xfail(vtkVersion.GetVTKMajorVersion()<9, reason='missing "Direction" attribute')
    def test_export_vtk(self,request,tmp_path,res_path,update,patch_execution_stamp,patch_datetime_now,output,fname,inc):
        result = Result(res_path/fname).view(increments=inc)
        result.export_VTK(output,target_dir=tmp_path,parallel=False)
        fname = fname.split('.')[0]+f'_inc{(inc if type(inc) == int else inc[0]):0>2}.vti'
        v = VTK.load(tmp_path/fname)
        v.comments = ['n/a']
        v.save(tmp_path/fname,parallel=False)
        with open(tmp_path/fname,'rb') as f:
            cur = hashlib.md5(f.read()).hexdigest()
        if update:
            with open((res_path/'export_VTK'/request.node.name).with_suffix('.md5'),'w') as f:
                f.write(cur+'\n')
        with open((res_path/'export_VTK'/request.node.name).with_suffix('.md5')) as f:
            assert cur == f.read().strip('\n')

    @pytest.mark.parametrize('mode',['point','cell'])
    @pytest.mark.parametrize('output',[False,True])
    def test_export_vtk_marc(self,tmp_path,res_path,mode,output):
        os.chdir(tmp_path)
        result = Result(res_path/'check_compile_job1.hdf5')
        result.export_VTK(output,mode)

    def test_marc_coordinates(self,res_path):
        result = Result(res_path/'check_compile_job1.hdf5').view(increments=-1)
        c_n = result.coordinates0_node + result.get('u_n')
        c_p = result.coordinates0_point + result.get('u_p')
        assert len(c_n) > len(c_p)

    @pytest.mark.parametrize('mode',['point','cell'])
    def test_vtk_mode(self,tmp_path,single_phase,mode):
        os.chdir(tmp_path)
        single_phase.export_VTK(mode=mode)

    def test_vtk_invalid_mode(self,single_phase):
        with pytest.raises(ValueError):
            single_phase.export_VTK(mode='invalid')

    def test_vtk_custom_path(self,tmp_path,single_phase):
        export_dir = tmp_path/'export_dir'
        single_phase.export_VTK(mode='point',target_dir=export_dir,parallel=False)
        assert set(os.listdir(export_dir)) == set([f'{single_phase.fname.stem}_inc{i:02}.vtp' for i in range(0,40+1,4)])

    def test_XDMF_datatypes(self,tmp_path,single_phase,update,res_path):
        for what,shape in {'scalar':(),'vector':(3,),'tensor':(3,3),'matrix':(12,)}.items():
            for dtype in ['f4','f8','i1','i2','i4','i8','u1','u2','u4','u8']:
                single_phase.add_calculation(f"np.ones(np.shape(#F#)[0:1]+{shape},'{dtype}')",f'{what}_{dtype}')
        xdmf_path = tmp_path/single_phase.fname.with_suffix('.xdmf').name
        single_phase.export_XDMF(target_dir=tmp_path)
        if update:
            shutil.copy(xdmf_path,res_path/xdmf_path.name)
        assert sorted(open(xdmf_path).read()) == sorted(open(res_path/xdmf_path.name).read())

    @pytest.mark.skipif(not hasattr(vtkXdmfReader,'GetOutput'),reason='https://discourse.vtk.org/t/2450')
    def test_XDMF_shape(self,tmp_path,single_phase):
        single_phase.export_XDMF(target_dir=single_phase.fname.parent)
        fname = single_phase.fname.with_suffix('.xdmf')
        reader_xdmf = vtkXdmfReader()
        reader_xdmf.SetFileName(fname)
        reader_xdmf.Update()
        dim_xdmf = reader_xdmf.GetOutput().GetDimensions()
        bounds_xdmf = reader_xdmf.GetOutput().GetBounds()

        single_phase.view(increments=0).export_VTK(target_dir=single_phase.fname.parent,parallel=False)
        fname = single_phase.fname.with_name(single_phase.fname.stem+'_inc00.vti')
        reader_vti = vtkXMLImageDataReader()
        reader_vti.SetFileName(fname)
        reader_vti.Update()
        dim_vti = reader_vti.GetOutput().GetDimensions()
        bounds_vti = reader_vti.GetOutput().GetBounds()
        assert dim_vti == dim_xdmf and bounds_vti == bounds_xdmf

    def test_XDMF_invalid(self,default):
        with pytest.raises(TypeError):
            default.export_XDMF()

    def test_XDMF_custom_path(self,single_phase,tmp_path):
        os.chdir(tmp_path)
        single_phase.export_XDMF()
        assert single_phase.fname.with_suffix('.xdmf').name in os.listdir(tmp_path)
        export_dir = tmp_path/'export_dir'
        single_phase.export_XDMF(target_dir=export_dir)
        assert single_phase.fname.with_suffix('.xdmf').name in os.listdir(export_dir)

    @pytest.mark.skipif(not hasattr(vtkXdmfReader,'GetOutput'),reason='https://discourse.vtk.org/t/2450')
    def test_XDMF_relabs_path(self,single_phase,tmp_path):
        def dims(xdmf):
            reader_xdmf = vtkXdmfReader()
            reader_xdmf.SetFileName(xdmf)
            reader_xdmf.Update()
            return reader_xdmf.GetOutput().GetDimensions()

        single_phase.export_XDMF(target_dir=tmp_path)
        xdmfname = single_phase.fname.with_suffix('.xdmf').name
        ref_dims = dims(tmp_path/xdmfname)

        for (d,info) in {
             'A': dict(absolute_path=True,
                       mv='..',
                       ),
             'B': dict(absolute_path=False,
                       mv='../A',
                       ),
            }.items():
            sub = tmp_path/d; sub.mkdir(exist_ok=True)
            single_phase.export_XDMF(target_dir=sub,absolute_path=info['absolute_path'])
            os.replace(sub/xdmfname,sub/info['mv']/xdmfname)
            assert ref_dims == dims(sub/info['mv']/xdmfname)

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
    def test_get(self,update,request,res_path,view,output,flatten,prune):
        result = Result(res_path/'4grains2x4x3_compressionY.hdf5')
        for key,value in view.items():
            result = result.view(**{key:value})

        fname = request.node.name
        cur = result.get(output,flatten,prune)
        if update:
            with bz2.BZ2File((res_path/'get'/fname).with_suffix('.pbz2'),'w') as f:
                pickle.dump(cur,f)

        with bz2.BZ2File((res_path/'get'/fname).with_suffix('.pbz2')) as f:
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
    def test_place(self,update,request,res_path,view,output,flatten,prune,constituents):
        result = Result(res_path/'4grains2x4x3_compressionY.hdf5')
        for key,value in view.items():
            result = result.view(**{key:value})

        fname = request.node.name
        cur = result.place(output,flatten,prune,constituents)
        if update:
            with bz2.BZ2File((res_path/'place'/fname).with_suffix('.pbz2'),'w') as f:
                pickle.dump(cur,f)

        with bz2.BZ2File((res_path/'place'/fname).with_suffix('.pbz2')) as f:
            ref = pickle.load(f)
            assert cur is None if ref is None else dict_equal(cur,ref)

    def test_simulation_setup_files(self,default):
        assert set(default.simulation_setup_files) == set(['12grains6x7x8.vti',
                                                            'material.yaml',
                                                            'tensionY.yaml',
                                                            'previous/12grains6x7x8.vti',
                                                            'previous/material.yaml',
                                                            'previous/tensionY.yaml'])

    def test_export_simulation_setup_files(self,tmp_path,default):
        sub = 'deep/down'
        default.export_simulation_setup(target_dir=tmp_path/sub,overwrite=True)
        for f in default.simulation_setup_files:
            assert (tmp_path/sub/f).exists()

    def test_export_simulation_setup_overwrite(self,tmp_path,default):
        os.chdir(tmp_path)
        default.export_simulation_setup('material.yaml',overwrite=True)
        with pytest.raises(PermissionError):
            default.export_simulation_setup('material.yaml',overwrite=False)

    @pytest.mark.parametrize('output',['12grains6x7x8.vti',
                                       'tensionY.yaml',
                                      ])
    def test_export_simulation_setup_content(self,res_path,tmp_path,default,output):
        default.export_simulation_setup(output,target_dir=tmp_path,overwrite=True)
        assert open(tmp_path/output).read() == open(res_path/output).read()

    @pytest.mark.parametrize('fname',['4grains2x4x3_compressionY.hdf5',
                                      '6grains6x7x8_single_phase_tensionY.hdf5'])
    @pytest.mark.parametrize('output',['material.yaml','*'])
    def test_export_simulation_setup_consistency(self,res_path,tmp_path,fname,output):
        r = Result(res_path/fname)
        r.export_simulation_setup(output,target_dir=tmp_path)
        with h5py.File(res_path/fname,'r') as f_hdf5:
            for file in fnmatch.filter(f_hdf5['setup'].keys(),output):
                with open(tmp_path/file) as f:
                    assert f_hdf5[f'setup/{file}'][()][0].decode() == f.read()

    def test_export_simulation_setup_custom_path(self,res_path,tmp_path):
        subdir = 'export_dir'
        absdir = tmp_path/subdir
        absdir.mkdir(exist_ok=True)

        r = Result(res_path/'4grains2x4x3_compressionY.hdf5')
        for t,cwd in zip([absdir,subdir,None],[tmp_path,tmp_path,absdir]):
            os.chdir(cwd)
            r.export_simulation_setup('material.yaml',target_dir=t)
            assert 'material.yaml' in os.listdir(absdir); (absdir/'material.yaml').unlink()

    @pytest.mark.parametrize('fname',['4grains2x4x3_compressionY.hdf5',
                                      '6grains6x7x8_single_phase_tensionY.hdf5',
                                      '12grains6x7x8_tensionY.hdf5',
                                      'check_compile_job1.hdf5',])
    def test_export_DADF5(self,res_path,tmp_path,fname):
        r = Result(res_path/fname)
        r = r.view(phases = random.sample(r.phases,1))
        r = r.view(increments = random.sample(r.increments,np.random.randint(1,len(r.increments))))
        r.export_DADF5(tmp_path/fname)
        r_exp = Result(tmp_path/fname)
        assert str(r.get()) == str(r_exp.get())
        assert str(r.place()) == str(r_exp.place())

    @pytest.mark.parametrize('fname',['4grains2x4x3_compressionY.hdf5',
                                      '6grains6x7x8_single_phase_tensionY.hdf5'])
    def test_export_DADF5_name_clash(self,res_path,tmp_path,fname):
        r = Result(res_path/fname)
        with pytest.raises(PermissionError):
            r.export_DADF5(r.fname)

    @pytest.mark.parametrize('fname',['4grains2x4x3_compressionY.hdf5',
                                      '6grains6x7x8_single_phase_tensionY.hdf5',
                                      '12grains6x7x8_tensionY.hdf5'])
    def test_export_DADF5_regrid(self,res_path,tmp_path,fname):
        r = Result(res_path/fname)
        m = grid_filters.regrid(r.size,np.broadcast_to(np.eye(3),tuple(r.cells)+(3,3)),r.cells*2)
        r.export_DADF5(tmp_path/'regridded.hdf5',mapping=m)
        assert np.all(Result(tmp_path/'regridded.hdf5').cells == r.cells*2)
