# DAMASK patching

This folder contains patches that modify the functionality of the current version of DAMASK prior to the corresponding inclusion in the official release.

## Usage

```bash
cd DAMASK_ROOT
patch -p1 installation/patch/nameOfPatch
```

## Available patches

  * **fwbw_derivative** switches the default spatial derivative from continuous to forward/backward difference.  
    This generally reduces spurious oscillations in the result as the spatial accuracy of the derivative is then compatible with the underlying solution grid.
  * **petsc3.7** adapts to API changes introduced between PETSc 3.6.x and 3.7.x for setting PETSc options.  
    Use this patch if your system runs PETSc 3.7.x.
