# SPDX-License-Identifier: AGPL-3.0-or-later
import os
env = os.environ                                                                                    # without MPI (might come with h5py)
import subprocess
import random
import difflib

import pytest

@pytest.mark.parametrize('solver',['damask_mesh','damask_grid'])
def test_output_redirection(tmp_path,mpi_launcher,solver):

    env['DAMASK_LOGFILE'] = '0'
    ref = subprocess.run(solver,cwd=tmp_path,env=env,capture_output=True).stdout.decode('utf-8').splitlines(True)
    subprocess.run(f'{solver} > log_single',cwd=tmp_path,env=env,shell=True)
    subprocess.run(f'{mpi_launcher} -np 2 {solver} > log_MPI',cwd=tmp_path,env=env,shell=True)

    env['DAMASK_LOGFILE'] = random.choice(['1','True','TRUE','true'])
    subprocess.run(f'{mpi_launcher} -np 2 {solver} > log_MPI_empty',cwd=tmp_path,env=env,shell=True)
    subprocess.run(f'{solver} > log_single_empty',cwd=tmp_path,env=env,shell=True)

    allowed_deviations = ['Date', 'Time', 'MPI worldrank', 'MPI worldsize']
    for fname in ['log_single', 'log_MPI','out.0000','out.0001']:
        with open(tmp_path/fname) as f:
            for line_ref,line_cur in zip(ref,f):
                if fname == 'log_MPI' and line_cur.strip().startswith('\u001b'): continue            # TTY cannot be detected on MPI
                assert line_ref==line_cur or set([line_ref.strip().split(':')[0],
                                                  line_cur.strip().split(':')[0]]).issubset(allowed_deviations)

    for fname in ['log_single_empty', 'log_MPI_empty']:
        assert os.path.getsize(tmp_path/fname) == 0
