import time
import shutil
import os
import sys
from datetime import datetime

import pytest
import numpy as np
import h5py

from damask import Result
from damask import Rotation
from damask import Orientation
from damask import mechanics
from damask import grid_filters

@pytest.fixture
def default(tmp_path,reference_dir):
    """Small Result file in temp location for modification."""
    fname = '12grains6x7x8_tensionY.hdf5'
    shutil.copy(reference_dir/fname,tmp_path)
    f = Result(tmp_path/fname)
    f.pick('times',20.0)
    return f

@pytest.fixture
def single_phase(tmp_path,reference_dir):
    """Single phase Result file in temp location for modification."""
    fname = '6grains6x7x8_single_phase_tensionY.hdf5'
    shutil.copy(reference_dir/fname,tmp_path)
    return Result(tmp_path/fname)

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return reference_dir_base/'Result'


class TestResult:

    def test_self_report(self,default):
        print(default)


    def test_pick_all(self,default):
        default.pick('increments',True)
        a = default.get_dataset_location('F')
        default.pick('increments','*')
        b = default.get_dataset_location('F')
        default.pick('increments',default.incs_in_range(0,np.iinfo(int).max))
        c = default.get_dataset_location('F')

        default.pick('times',True)
        d = default.get_dataset_location('F')
        default.pick('times','*')
        e = default.get_dataset_location('F')
        default.pick('times',default.times_in_range(0.0,np.inf))
        f = default.get_dataset_location('F')
        assert a == b == c == d == e ==f

    @pytest.mark.parametrize('what',['increments','times','constituents'])                          # ToDo: discuss materialpoints
    def test_pick_none(self,default,what):
        default.pick(what,False)
        a = default.get_dataset_location('F')
        default.pick(what,[])
        b = default.get_dataset_location('F')

        assert a == b == []

    @pytest.mark.parametrize('what',['increments','times','constituents'])                          # ToDo: discuss materialpoints
    def test_pick_more(self,default,what):
        default.pick(what,False)
        default.pick_more(what,'*')
        a = default.get_dataset_location('F')

        default.pick(what,True)
        b = default.get_dataset_location('F')

        assert a == b

    @pytest.mark.parametrize('what',['increments','times','constituents'])                          # ToDo: discuss materialpoints
    def test_pick_less(self,default,what):
        default.pick(what,True)
        default.pick_less(what,'*')
        a = default.get_dataset_location('F')

        default.pick(what,False)
        b = default.get_dataset_location('F')

        assert a == b == []

    def test_pick_invalid(self,default):
        with pytest.raises(AttributeError):
            default.pick('invalid',True)

    def test_add_absolute(self,default):
        default.add_absolute('Fe')
        loc = {'Fe':   default.get_dataset_location('Fe'),
               '|Fe|': default.get_dataset_location('|Fe|')}
        in_memory = np.abs(default.read_dataset(loc['Fe'],0))
        in_file   = default.read_dataset(loc['|Fe|'],0)
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('mode',['direct','function'])
    def test_add_calculation(self,default,tmp_path,mode):

        if mode == 'direct':
            default.add_calculation('x','2.0*np.abs(#F#)-1.0','-','my notes')
        else:
            with open(tmp_path/'f.py','w') as f:
                f.write("import numpy as np\ndef my_func(field):\n  return 2.0*np.abs(field)-1.0\n")
            sys.path.insert(0,str(tmp_path))
            import f
            default.enable_user_function(f.my_func)
            default.add_calculation('x','my_func(#F#)','-','my notes')

        loc = {'F':    default.get_dataset_location('F'),
               'x':    default.get_dataset_location('x')}
        in_memory = 2.0*np.abs(default.read_dataset(loc['F'],0))-1.0
        in_file   = default.read_dataset(loc['x'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_Cauchy(self,default):
        default.add_Cauchy('P','F')
        loc = {'F':    default.get_dataset_location('F'),
               'P':    default.get_dataset_location('P'),
               'sigma':default.get_dataset_location('sigma')}
        in_memory = mechanics.Cauchy(default.read_dataset(loc['P'],0),
                                     default.read_dataset(loc['F'],0))
        in_file   = default.read_dataset(loc['sigma'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_determinant(self,default):
        default.add_determinant('P')
        loc = {'P':     default.get_dataset_location('P'),
               'det(P)':default.get_dataset_location('det(P)')}
        in_memory = np.linalg.det(default.read_dataset(loc['P'],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc['det(P)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_deviator(self,default):
        default.add_deviator('P')
        loc = {'P'  :default.get_dataset_location('P'),
               's_P':default.get_dataset_location('s_P')}
        in_memory = mechanics.deviatoric_part(default.read_dataset(loc['P'],0))
        in_file   = default.read_dataset(loc['s_P'],0)
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('eigenvalue,function',[('max',np.amax),('min',np.amin)])
    def test_add_eigenvalue(self,default,eigenvalue,function):
        default.add_Cauchy('P','F')
        default.add_eigenvalue('sigma',eigenvalue)
        loc = {'sigma' :default.get_dataset_location('sigma'),
               'lambda':default.get_dataset_location(f'lambda_{eigenvalue}(sigma)')}
        in_memory = function(mechanics.eigenvalues(default.read_dataset(loc['sigma'],0)),axis=1,keepdims=True)
        in_file   = default.read_dataset(loc['lambda'],0)
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('eigenvalue,idx',[('max',2),('mid',1),('min',0)])
    def test_add_eigenvector(self,default,eigenvalue,idx):
        default.add_Cauchy('P','F')
        default.add_eigenvector('sigma',eigenvalue)
        loc = {'sigma'   :default.get_dataset_location('sigma'),
               'v(sigma)':default.get_dataset_location(f'v_{eigenvalue}(sigma)')}
        in_memory = mechanics.eigenvectors(default.read_dataset(loc['sigma'],0))[:,idx]
        in_file   = default.read_dataset(loc['v(sigma)'],0)
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('d',[[1,0,0],[0,1,0],[0,0,1]])
    def test_add_IPF_color(self,default,d):
        default.add_IPF_color('orientation',d)
        loc = {'orientation': default.get_dataset_location('orientation'),
               'color':       default.get_dataset_location('IPFcolor_[{} {} {}]'.format(*d))}
        qu = default.read_dataset(loc['orientation']).view(np.double).reshape(-1,4)
        crystal_structure = default.get_crystal_structure()
        in_memory = np.empty((qu.shape[0],3),np.uint8)
        for i,q in enumerate(qu):
            o = Orientation(q,crystal_structure).reduced
            in_memory[i] = np.uint8(o.IPF_color(np.array(d))*255)
        in_file = default.read_dataset(loc['color'])
        assert np.allclose(in_memory,in_file)

    def test_add_maximum_shear(self,default):
        default.add_Cauchy('P','F')
        default.add_maximum_shear('sigma')
        loc = {'sigma'           :default.get_dataset_location('sigma'),
               'max_shear(sigma)':default.get_dataset_location('max_shear(sigma)')}
        in_memory = mechanics.maximum_shear(default.read_dataset(loc['sigma'],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc['max_shear(sigma)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_Mises_strain(self,default):
        t = ['V','U'][np.random.randint(0,2)]
        m = np.random.random()*2.0 - 1.0
        default.add_strain_tensor('F',t,m)
        label = f'epsilon_{t}^{m}(F)'
        default.add_Mises(label)
        loc = {label      :default.get_dataset_location(label),
               label+'_vM':default.get_dataset_location(label+'_vM')}
        in_memory = mechanics.Mises_strain(default.read_dataset(loc[label],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc[label+'_vM'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_Mises_stress(self,default):
        default.add_Cauchy('P','F')
        default.add_Mises('sigma')
        loc = {'sigma'   :default.get_dataset_location('sigma'),
               'sigma_vM':default.get_dataset_location('sigma_vM')}
        in_memory = mechanics.Mises_stress(default.read_dataset(loc['sigma'],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc['sigma_vM'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_norm(self,default):
        default.add_norm('F',1)
        loc = {'F':    default.get_dataset_location('F'),
               '|F|_1':default.get_dataset_location('|F|_1')}
        in_memory = np.linalg.norm(default.read_dataset(loc['F'],0),ord=1,axis=(1,2),keepdims=True)
        in_file   = default.read_dataset(loc['|F|_1'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_PK2(self,default):
        default.add_PK2('P','F')
        loc = {'F':default.get_dataset_location('F'),
               'P':default.get_dataset_location('P'),
               'S':default.get_dataset_location('S')}
        in_memory = mechanics.PK2(default.read_dataset(loc['P'],0),
                                  default.read_dataset(loc['F'],0))
        in_file   = default.read_dataset(loc['S'],0)
        assert np.allclose(in_memory,in_file)

    @pytest.mark.parametrize('polar',[True,False])
    def test_add_pole(self,default,polar):
        pole = np.array([1.,0.,0.])
        default.add_pole('orientation',pole,polar)
        loc = {'orientation': default.get_dataset_location('orientation'),
               'pole':        default.get_dataset_location('p^{}_[1 0 0)'.format(u'rφ' if polar else 'xy'))}
        rot = Rotation(default.read_dataset(loc['orientation']).view(np.double))
        rotated_pole = rot * np.broadcast_to(pole,rot.shape+(3,))
        xy = rotated_pole[:,0:2]/(1.+abs(pole[2]))
        in_memory = xy if not polar else \
                    np.block([np.sqrt(xy[:,0:1]*xy[:,0:1]+xy[:,1:2]*xy[:,1:2]),np.arctan2(xy[:,1:2],xy[:,0:1])])
        in_file = default.read_dataset(loc['pole'])
        assert np.allclose(in_memory,in_file)

    def test_add_rotational_part(self,default):
        default.add_rotational_part('F')
        loc = {'F':    default.get_dataset_location('F'),
               'R(F)': default.get_dataset_location('R(F)')}
        in_memory = mechanics.rotational_part(default.read_dataset(loc['F'],0))
        in_file   = default.read_dataset(loc['R(F)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_spherical(self,default):
        default.add_spherical('P')
        loc = {'P':   default.get_dataset_location('P'),
               'p_P': default.get_dataset_location('p_P')}
        in_memory = mechanics.spherical_part(default.read_dataset(loc['P'],0)).reshape(-1,1)
        in_file   = default.read_dataset(loc['p_P'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_strain(self,default):
        t = ['V','U'][np.random.randint(0,2)]
        m = np.random.random()*2.0 - 1.0
        default.add_strain_tensor('F',t,m)
        label = f'epsilon_{t}^{m}(F)'
        loc = {'F':   default.get_dataset_location('F'),
               label: default.get_dataset_location(label)}
        in_memory = mechanics.strain_tensor(default.read_dataset(loc['F'],0),t,m)
        in_file   = default.read_dataset(loc[label],0)
        assert np.allclose(in_memory,in_file)

    def test_add_stretch_right(self,default):
        default.add_stretch_tensor('F','U')
        loc = {'F':    default.get_dataset_location('F'),
               'U(F)': default.get_dataset_location('U(F)')}
        in_memory = mechanics.right_stretch(default.read_dataset(loc['F'],0))
        in_file   = default.read_dataset(loc['U(F)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_stretch_left(self,default):
        default.add_stretch_tensor('F','V')
        loc = {'F':    default.get_dataset_location('F'),
               'V(F)': default.get_dataset_location('V(F)')}
        in_memory = mechanics.left_stretch(default.read_dataset(loc['F'],0))
        in_file   = default.read_dataset(loc['V(F)'],0)
        assert np.allclose(in_memory,in_file)

    def test_add_invalid(self,default):
        with pytest.raises(TypeError):
            default.add_calculation('#invalid#*2')

    @pytest.mark.parametrize('overwrite',['off','on'])
    def test_add_overwrite(self,default,overwrite):
        default.pick('times',default.times_in_range(0,np.inf)[-1])

        default.add_Cauchy()
        loc = default.get_dataset_location('sigma')
        with h5py.File(default.fname,'r') as f:
            created_first = f[loc[0]].attrs['Created'].decode()
        created_first = datetime.strptime(created_first,'%Y-%m-%d %H:%M:%S%z')

        if overwrite == 'on':
            default.allow_modification()
        else:
            default.disallow_modification()

        time.sleep(2.)
        default.add_calculation('sigma','#sigma#*0.0+311.','not the Cauchy stress')
        with h5py.File(default.fname,'r') as f:
            created_second = f[loc[0]].attrs['Created'].decode()
        created_second = datetime.strptime(created_second,'%Y-%m-%d %H:%M:%S%z')
        if overwrite == 'on':
            assert created_first < created_second and np.allclose(default.read_dataset(loc),311.)
        else:
            assert created_first == created_second and not np.allclose(default.read_dataset(loc),311.)

    @pytest.mark.parametrize('allowed',['off','on'])
    def test_rename(self,default,allowed):
        if allowed == 'on':
            F = default.read_dataset(default.get_dataset_location('F'))
            default.allow_modification()
            default.rename('F','new_name')
            assert np.all(F == default.read_dataset(default.get_dataset_location('new_name')))
            default.disallow_modification()

        with pytest.raises(PermissionError):
            default.rename('P','another_new_name')

    @pytest.mark.parametrize('mode',['cell','node'])
    def test_coordinates(self,default,mode):
         if   mode == 'cell':
             a = grid_filters.cell_coord0(default.grid,default.size,default.origin)
             b = default.cell_coordinates.reshape(tuple(default.grid)+(3,),order='F')
         elif mode == 'node':
             a = grid_filters.node_coord0(default.grid,default.size,default.origin)
             b = default.node_coordinates.reshape(tuple(default.grid+1)+(3,),order='F')
         assert np.allclose(a,b)

    @pytest.mark.parametrize('output',['F',[],['F','P']])
    def test_vtk(self,tmp_path,default,output):
        os.chdir(tmp_path)
        default.to_vtk(output)

    def test_XDMF(self,tmp_path,single_phase):
        os.chdir(tmp_path)
        single_phase.write_XDMF()
