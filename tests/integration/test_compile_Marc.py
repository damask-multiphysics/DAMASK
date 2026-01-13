# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import pytest

from damask import solver

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'compile_Marc'


@pytest.mark.parametrize('optimization',['','l','h'])
def test_compile_marc(damask_root,res_path,tmp_path,copy_files,optimization):
    copy_files(res_path,tmp_path)

    s = solver.Marc(damask_root=damask_root)
    os.chdir(tmp_path)
    s.submit_job('check_compile','job1',True,optimization=optimization)
