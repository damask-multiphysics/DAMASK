#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask
import pipes

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', 
                      description ='generate 3D RVE from .ang files of EBSD slices .',
                      version = scriptID)

parser.add_option('--offset',
                  dest='offset',
                  type='float',
                  help='offset of EBSD slices [%default]',
                  metavar='float')
parser.add_option('--outname',
                  dest='outName',
                  type='string',
                  help='output file name [%default]',  metavar='string')    
parser.add_option('--vtr',
                  action="store_true",
                  dest='vtr')
parser.add_option('--geom',
                  action="store_true",
                  dest='geom')
parser.set_defaults(offset = 1.0,
                    outName = 'RVE3D')

(options,filenames) = parser.parse_args()

numFiles = len(filenames)
formatwidth = 1+int(math.log10(numFiles))

# copy original files to tmp files to not alter originals
for i in range(numFiles):
  sliceID = 'slice' + str(i).zfill(formatwidth) + '.tmp'
  strCommand = 'cp ' + pipes.quote(filenames[i]) + ' ' + sliceID
  os.system(strCommand)

# modify tmp files
print('Add z-coordinates')
for i in range(numFiles):
  sliceID = 'slice' + str(i).zfill(formatwidth) + '.tmp'
  strCommand = 'OIMgrainFile_toTable ' + sliceID
  os.system(strCommand)
  strCommand = 'addCalculation --label 3Dpos --formula "np.array(#pos#.tolist()+[' + str(i*options.offset) + '])" ' + sliceID
  os.system(strCommand)

# join temp files into one

print('\n Colocate files')
fileOut = open(options.outName + '.ang','w')

# take header information from 1st file
sliceID = 'slice' + str(0).zfill(formatwidth) + '.tmp'
fileRead = open(sliceID)
data = fileRead.readlines()
fileRead.close()
headerLines = int(data[0].split()[0])
fileOut.write(str(headerLines+1) + '\t header\n')
for line in data[1:headerLines]:
  fileOut.write(line)
fileOut.write(scriptID + '\t' + ' '.join(sys.argv[1:]) + '\n')
for line in data[headerLines:]:
  fileOut.write(line)

# append other files content without header
for i in range(numFiles-1):
  sliceID = 'slice' + str(i+1).zfill(formatwidth) + '.tmp'
  fileRead = open(sliceID)
  data = fileRead.readlines()
  fileRead.close()
  headerLines = int(data[0].split()[0])
  for line in data[headerLines+1:]:
    fileOut.write(line)
fileOut.close()

# tidy up and add phase column 
print('\n Remove temp data and add phase info')
strCommand = 'filterTable --black pos ' + options.outName + '.ang'
os.system(strCommand)
strCommand = 'reLabel --label 3Dpos --substitute pos ' + options.outName + '.ang'
os.system(strCommand)
strCommand = 'addCalculation -l phase -f 1 ' + options.outName + '.ang'
os.system(strCommand)


# create geom file when asked for
if options.geom:
  print('\n Build geometry file')
  strCommand = 'geom_fromTable --phase phase --eulers euler --coordinates pos ' + pipes.quote(options.outName) + '.ang'
  os.system(strCommand)

# create paraview file when asked for

if options.vtr:
  print('\n Build Paraview file')
  strCommand = 'addIPFcolor --eulers euler --pole 0.0 0.0 1.0 ' + options.outName + '.ang'
  os.system(strCommand)
  strCommand = 'vtk_rectilinearGrid ' + pipes.quote(options.outName) + '.ang'
  os.system(strCommand)
  os.rename(pipes.quote(options.outName) + '_pos(cell)'+'.vtr', pipes.quote(options.outName) + '.vtr')
  strCommand = 'vtk_addRectilinearGridData --vtk '+ pipes.quote(options.outName) + '.vtr --color IPF_001_cubic ' + pipes.quote(options.outName) + '.ang'
  os.system(strCommand)
  
# delete tmp files
for i in range(numFiles):
  sliceID = 'slice' + str(i).zfill(formatwidth) + '.tmp'
  os.remove(sliceID)