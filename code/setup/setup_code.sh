#!/usr/bin/env sh

${DAMASK_ROOT}/code/setup/modify_Files.py
${DAMASK_ROOT}/code/setup/compile_SpectralSolver.py $@
${DAMASK_ROOT}/code/setup/symlink_Code.py
