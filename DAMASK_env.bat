:: sets up an environment for DAMASK on Windows
:: usage:  call DAMASK_env.bat
@echo off
set LOCATION=%~dp0
set DAMASK_ROOT=%LOCATION%\DAMASK
set DAMASK_NUM_THREADS=2
chcp 1252
Title Düsseldorf Advanced Materials Simulation Kit - DAMASK, MPIE Düsseldorf
echo.
echo Düsseldorf Advanced Materials Simulation Kit - DAMASK
echo Max-Planck-Institut für Eisenforschung, Düsseldorf
echo http://damask.mpie.de
echo.
echo Preparing environment ...
echo DAMASK_ROOT=%DAMASK_ROOT%
echo DAMASK_NUM_THREADS=%DAMASK_NUM_THREADS%
set DAMASK_BIN=%DAMASK_ROOT%\bin
set PATH=%PATH%;%DAMASK_BIN%
set PYTHONPATH=%PYTHONPATH%;%DAMASK_ROOT%\lib
