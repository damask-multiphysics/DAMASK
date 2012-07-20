@echo off
chcp 1252
Title Düsseldorf Adavanced MAterials Simulation Kit - DAMASK, MPIE Düsseldorf

echo.
echo Düsseldorf Adavanced MAterials Simulation Kit - DAMASK
echo Max-Planck-Institut für Eisenforschung Düsseldorf
echo http://damask.mpie.de
echo.

::echo %cd%

echo Preparing environment ...
SET DAMASK_ROOT=%~dp0
::SET PYTHONPATH=%DAMASK_ROOT%\lib;%PYTHONPATH%
SET PYTHONPATH=%DAMASK_ROOT%\lib;%PYTHONPATH%

echo DAMASK_ROOT: %DAMASK_ROOT%
echo.

::python -c "import damask as D" -i
ipython -c "import damask as D; import os; os.chdir(os.getenv('DAMASK_ROOT')); print('pwd: '+os.getcwd());" -i
