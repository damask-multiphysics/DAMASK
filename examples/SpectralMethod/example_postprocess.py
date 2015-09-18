#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os
import glob
from subprocess import call

geom_name = '20grains16x16x16_tensionX'
postResults = 'postResults --cr f,p --split --separation x,y,z '+geom_name+'.spectralOut'

sts = call(postResults, shell=True)

os.chdir('./postProc/')
ascii_files = glob.glob(geom_name+'_inc*.txt')
print ascii_files

showTable = "showTable -a "
addCauchy = 'addCauchy '
addMises = 'addMises -s Cauchy '
addStrainTensors = "addStrainTensors -0 -v "
visualize3D = "3Dvisualize -s 'Mises(Cauchy)',1_p  Cauchy "


postProc = [addCauchy, addMises, addStrainTensors, visualize3D]


for f in ascii_files:
    print f
    for p in postProc:
        p = p+f
        print p
        sts = call(p,shell=True)

