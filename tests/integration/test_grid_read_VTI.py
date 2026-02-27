# SPDX-License-Identifier: AGPL-3.0-or-later
import random
import string
import tempfile
from pathlib import Path

import pytest
import numpy as np
import h5py
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.util.numpy_support import numpy_to_vtk as np_to_vtk

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_read_VTI'

temp_res_path = None


@pytest.mark.parametrize('compress',[True,False])
@pytest.mark.parametrize('dtype',['i4','i8','f4','f8'])
def test_grid_read_VTI(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,np_rng,
                        dtype,compress):
    global temp_res_path
    load = 'none'
    grid = 'grid'
    material = 'material'

    if temp_res_path is None:
        temp_res_path = Path(tempfile.mkdtemp())
        copy_files(res_path,temp_res_path)
        damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load}',
                        wd=temp_res_path)

    # create non-standard.vti
    geom = damask.GeomGrid.load(temp_res_path/f'{grid}.vti')
    vtk_data = vtkImageData()
    vtk_data.SetDimensions(*(geom.cells+1))
    origin = np_rng.random(3)
    vtk_data.SetOrigin(*origin)
    spacing = (geom.size+np_rng.random(3))/geom.cells
    vtk_data.SetSpacing(*spacing)
    ms = np_to_vtk(geom.material.flatten(order='F').astype(dtype))
    ms.SetName('material')
    vtk_data.GetCellData().AddArray(ms)
    noise = np_to_vtk(geom.material.flatten(order='F').astype(dtype)+np_rng.integers(1,32))
    noise.SetName(''.join(random.choices(string.ascii_letters, k=np_rng.integers(1,32))))
    vtk_data.GetCellData().AddArray(noise)

    v = damask.VTK(vtk_data)

    v.save(tmp_path/f'{grid}_{dtype}_{compress}.vti',compress=compress,parallel=False)
    copy_files(res_path,tmp_path,[f'{material}.yaml',f'{load}.yaml'])
    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}_{dtype}_{compress}.vti -m {material}.yaml '\
                    f'-j {grid}_{dtype}_{compress}_{load}',wd=tmp_path)

    with h5py.File(temp_res_path/f'{grid}_{load}.hdf5','r') as ref, \
         h5py.File(tmp_path/f'{grid}_{dtype}_{compress}_{load}.hdf5','r') as cur:

        for (path,dset) in h5py_dataset_iterator(ref):
            if len(dset.dtype) > 0:                                                                 # cell_to
                for e in dset.dtype.fields:
                    try:
                        assert_allclose(dset[e][()],cur[path][e][()])
                    except TypeError:
                        assert np.all(dset[e]==cur[path][e])
            elif dset.dtype.char == 'S':                                                            # setup
                pass
            elif path == 'geometry/cells':
                assert_allclose(geom.cells,cur[path])
            elif path == 'geometry/origin':
                assert_allclose(origin,cur[path])
            elif path == 'geometry/size':
                assert_allclose(spacing*geom.cells,cur[path])
            else:
                assert_allclose(dset,cur[path])

        with open(tmp_path/f'{grid}_{dtype}_{compress}.vti') as f:
            assert cur[f'setup/{grid}_{dtype}_{compress}.vti'][0].decode() == f.read()

def test_grid_read_VTI_VTK8(res_path,tmp_path,copy_files):
    load = 'none'
    grid = 'vtk8'
    material = 'material'

    copy_files(res_path,tmp_path,[f'{material}.yaml',f'{load}.yaml',f'{grid}.vti'])
    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load}',
                    wd=tmp_path)

    with h5py.File(tmp_path/f'{grid}_{load}.hdf5','r') as result, \
         open(tmp_path/f'{grid}.vti') as vti:
        assert result[f'setup/{grid}.vti'][0].decode() == vti.read()
