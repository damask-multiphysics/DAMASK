# DAMASK patching

This folder contains patches that modify the functionality of the current development version of DAMASK ahead of the corresponding adoption in the official release.

## Usage

```bash
cd DAMASK_ROOT
patch -p1 < installation/patch/nameOfPatch
```

## Create patch
commit your changes

```bash
git format-patch PATH_TO_COMPARE --stdout >
```
