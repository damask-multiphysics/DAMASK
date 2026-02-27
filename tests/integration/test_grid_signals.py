# SPDX-License-Identifier: AGPL-3.0-or-later
import os
env = os.environ                                                                                    # without MPI (might come with h5py)
import re
import subprocess
import shlex
import time
import signal as sg

import pytest

import damask

# https://software.intel.com/en-us/mpi-developer-reference-linux-hydra-environment-variables
env['I_MPI_JOB_SIGNAL_PROPAGATION'] = 'enable'

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_signals'


# MPI will not always work https://stackoverflow.com/questions/1145741
# Note: If the test hangs, xfail might not be sufficient. In that case, use a marker
# as implemented in and before 3b95297f0de99d6ecc2847c754b40625dfbae568.
@pytest.mark.parametrize('mpi',[False,pytest.param(True,marks=pytest.mark.xfail)])
@pytest.mark.parametrize('signal',[sg.SIGUSR1,sg.SIGUSR2,sg.SIGINT])
def test_grid_signals(res_path,tmp_path,copy_files,mpi_launcher,mpi,signal):
    grid = 'test'
    load = 'tensionX'
    material = 'material'
    job  = f'{grid}_{load}'

    copy_files(res_path,tmp_path)

    match = re.compile(r' Increment 17/25-1/1 @ Iteration 1≤1≤', re.U)
    cmd = f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}'
    if mpi: cmd = f'{mpi_launcher} -n 2 '+cmd

    start = time.time()
    p = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=tmp_path,env=env)
    while p.poll() is None and time.time()-start<30:
        line = p.stdout.readline().decode('utf-8')
        if re.search(match, line):
            p.send_signal(signal)
            break

    p.communicate()
    assert p.returncode == 0

    res = damask.Result(tmp_path/f'{job}.hdf5')

    if signal == sg.SIGUSR1:
        assert res.increments == ['increment_0', 'increment_17', 'increment_25']
    elif signal == sg.SIGUSR2:
        assert os.path.isfile(tmp_path/f'{job}_restart.hdf5')
    elif signal == sg.SIGINT:
        with open(tmp_path/f'{job}.sta') as f:
            assert f.readlines()[-1].split()[0] == '17'
    else:
        assert False
