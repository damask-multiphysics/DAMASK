#! /usr/bin/env python

import damask
import numpy as np

rng = np.random.default_rng(20191102)
for s in [2,4,8,16,32]:
    cells = (s*8,s,s)
    size = (1.,1./8.,1./8.)
    N_cells = np.prod(cells)
    g = damask.GeomGrid(np.arange(N_cells).reshape(cells,order='F'),size)
    g.save(f'slab_{s}.vti')
    damask.GeomMesh.from_grid(g).save(f'slab_{s}.msh')

    cfg = damask.ConfigMaterial.load('material_base.yaml')
    for m in cfg['phase'].keys():
        O = damask.Rotation.from_random((N_cells//2,2),rng=rng)
        cfg = cfg.material_add(O=O,homogenization='isostrain',phase=[m]*2)
    cfg.save(f'material_{s}.yaml')
