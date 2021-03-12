#!/usr/bin/env python3

import os
import argparse

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('filenames', nargs='+',
                    help='DADF5 files')
parser.add_argument('-d','--dir', dest='dir',default='postProc',metavar='string',
                    help='name of subdirectory relative to the location of the DADF5 file to hold output')
parser.add_argument('--mat', nargs='+',
                    help='labels for homogenization',dest='mat')
parser.add_argument('--con', nargs='+',
                    help='labels for phase',dest='con')

options = parser.parse_args()

if options.mat is None: options.mat=[]
if options.con is None: options.con=[]

for filename in options.filenames:
    results = damask.Result(filename)

    if not results.structured: continue
    coords = damask.grid_filters.coordinates0_point(results.cells,results.size,results.origin).reshape(-1,3,order='F')

    N_digits = int(np.floor(np.log10(int(results.increments[-1][3:]))))+1
    N_digits = 5 # hack to keep test intact
    for inc in damask.util.show_progress(results.iterate('increments'),len(results.increments)):
        table = damask.Table(np.ones(np.product(results.cells),dtype=int)*int(inc[3:]),{'inc':(1,)})\
                      .add('pos',coords.reshape(-1,3))

        results.view('homogenizations',False)
        results.view('phases',True)
        for label in options.con:
            x = results.get_dataset_location(label)
            if len(x) != 0:
                table = table.add(label,results.read_dataset(x,0,plain=True).reshape(results.cells.prod(),-1))

        results.view('phases',False)
        results.view('homogenizations',True)
        for label in options.mat:
            x = results.get_dataset_location(label)
            if len(x) != 0:
                table = table.add(label,results.read_dataset(x,0,plain=True).reshape(results.cells.prod(),-1))

        dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
        if not os.path.isdir(dirname):
            os.mkdir(dirname,0o755)
        file_out = '{}_inc{}.txt'.format(os.path.splitext(os.path.split(filename)[-1])[0],
                                         inc[3:].zfill(N_digits))
        table.save(os.path.join(dirname,file_out),legacy=True)
