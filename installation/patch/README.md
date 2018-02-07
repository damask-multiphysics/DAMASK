# DAMASK patching

This folder contains patches that modify the functionality of the current development version of DAMASK ahead of the corresponding adoption in the official release.

## Usage

```bash
cd DAMASK_ROOT
patch -p1 < installation/patch/nameOfPatch
```

## Available patches

  * **fwbw_derivative** switches the default spatial derivative from continuous to forward/backward difference.  
    This generally reduces spurious oscillations in the result as the spatial accuracy of the derivative is then compatible with the underlying solution grid.

  * **PETSc-3.8** adjusts all includes and calls to PETSc to follow the 3.8.x API.
    This allows to use the most recent version of PETSc.

## Create patch
commit your changes

```bash
git format-patch PATH_TO_COMPARE --stdout >
```
