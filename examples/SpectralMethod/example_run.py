#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

from subprocess import call
call('DAMASK_spectral -l tensionX.load -g 20grains16x16x16.geom', shell=True)