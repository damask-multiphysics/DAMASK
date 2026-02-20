# SPDX-License-Identifier: AGPL-3.0-or-later
import shutil

import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'YAML_parse'


@pytest.mark.parametrize('fname',np.arange(7))
def test_YAML_parse(res_path,tmp_path,copy_files,fname):

    load = 'none'
    grid = 'simple'
    material = 'material'

    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.vti'])
    shutil.copy(res_path/f'{fname}.yaml',tmp_path/f'{material}.yaml')

    stdout,stderr = damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml',wd=tmp_path)


def test_CRLF(res_path,tmp_path,copy_files):

    load = 'none'
    grid = 'simple'
    material = 'material'

    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.vti'])
    shutil.copy(res_path/'CRLF.yaml',tmp_path/f'{material}.yaml')

    stdout,stderr = damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml',wd=tmp_path)
