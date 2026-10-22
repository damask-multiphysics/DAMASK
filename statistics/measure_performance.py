#!/usr/bin/env python3

import os
import sys
import time
import tempfile
import argparse
import shutil
import re
import platform
import multiprocessing
import subprocess
from pathlib import Path

import pandas as pd
import yaml

datafile = 'performance.txt'

parser = argparse.ArgumentParser(Path(__file__).name,
                  description='Measure performance of grid solver',
                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--Nruns', dest='N_runs', metavar='int', type=int,
                    default=2,
                    help='Number of runs')
parser.add_argument('--Nthreads', dest='N_threads', metavar='int', type=int,
                    default=2,
                    help='Number of threads (OMP_NUM_THREADS)')
parser.add_argument('--tag', dest='tag', metavar='string',
                    default=None,
                    help='Tag (label) for current test')
parser.add_argument('--output_dir', dest='output_dir', metavar='string',
                    default='.',
                    help=f'Directory for results, i.e. where "{datafile}" is stored.')
parser.add_argument('--input_dir', dest='input_dir', metavar='string',
                    help='Directory containing test data (grid, load, material.yaml, numerics.yaml)')

args = parser.parse_args()
input_dir  = Path(args.input_dir ).expanduser().resolve(strict=True)
output_dir = Path(args.output_dir).expanduser().resolve(strict=True)

solvers = ['spectral_basic','spectral_polarization','spectral_Galerkin']
load = 'tensionX'
columns = ['epoch','solver','wall time/s','system time/s','user time/s','memory/kB']

time_cmd = {
            'linux': '/usr/bin/time',
            'darwin': '/opt/homebrew/bin/gtime',                                                    # needs `brew install gnu-time`
           }[sys.platform]

try:
    with open('/proc/cpuinfo') as f:
        cpuinfo = f.read()
    cpu = re.search('model name.*: (.*)',cpuinfo).group(1)
except FileNotFoundError:
    cpu = platform.processor()

damask_out = subprocess.run(['damask_grid','-h'],check=True,capture_output=True).stdout.decode()
tag = re.search('Version: (.*)',damask_out).group(1) if args.tag is None else args.tag
fortran_compiler = re.search(r'Fortran Compiler:[ \t]*(?P<compiler>[^\r\n]+)', damask_out)
if fortran_compiler is None:
    raise RuntimeError('Fortran compiler not found in damask_grid output')

info = {'tag':tag,
        'OMP_NUM_THREADS':args.N_threads,
        'compiler':fortran_compiler.group('compiler').strip(),
        'PETSc': re.search('PETSc version: (.*)',damask_out).group(1),
        'OS': f'{platform.system()} {platform.release()}',
        'processor': cpu if cpu else 'unknown',
        'N_proc': multiprocessing.cpu_count()}


for solver in solvers:
    for run in range(args.N_runs):
        with tempfile.TemporaryDirectory() as wd:
            for config in ['material.yaml', 'numerics.yaml']:
                shutil.copy(input_dir/config,wd)

            with open((input_dir/load).with_suffix('.yaml')) as f:
                config_load = yaml.safe_load(f)
            config_load['solver'] = {'mechanical':solver}
            with open((input_dir/load).with_suffix('.yaml'),'w') as f:
                yaml.dump(config_load,f)

            cmd = [time_cmd, '--format=%e,%S,%U,%M', 'damask_grid',
                   '-g', f'{input_dir/"20grains16x16x16.vti"}',
                   '-l', f'{input_dir/load}.yaml',
                   '-m', f'{input_dir/"material.yaml"}',
                   '-n', f'{input_dir/"numerics.yaml"}']
            env = {**os.environ,**{'OMP_NUM_THREADS':str(args.N_threads)}}
            complete_process = subprocess.run(cmd,check=True,capture_output=True,env=env)
            out = complete_process.stdout.decode()
            err = complete_process.stderr.decode()

        assert int(re.search('OMP_NUM_THREADS: (.*)',out).group(1)) == args.N_threads
        assert solver in out

        t_wall,t_sys,t_user,memory = err.split('\n')[-2].split(',')
        measurement = dict(zip(columns,[time.time(),solver,t_wall,t_sys,t_user,memory]))
        for k,v in measurement.items():
            print(f'{k}: {v}')

        fname = output_dir/datafile
        try:
            d = pd.read_csv(fname,sep='\t',index_col=None)
        except FileNotFoundError:
            d = pd.DataFrame(columns=columns+list(info.keys()))
        pd.concat([d,pd.DataFrame({**measurement,**info},index=[0])],ignore_index=True).to_csv(fname,sep='\t',index=False)
