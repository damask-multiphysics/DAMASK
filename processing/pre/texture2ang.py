#!/usr/bin/env python

import damask
import os,sys,math,re,string
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1]

def integerFactorization(i):
  
  j = int(math.floor(math.sqrt(float(i))))
  while (j>1 and int(i)%j != 0):
    j -= 1
  return j

def positiveRadians(angle):

  angle = math.radians(float(angle))
  while angle < 0.0:
    angle += 2.0*math.pi

  return angle

def getHeader(sizeX,sizeY,step,format):
  if format == 'ang':
    return [ 
    '# TEM_PIXperUM          1.000000', 
    '# x-star                0.509548', 
    '# y-star                0.795272', 
    '# z-star                0.611799', 
    '# WorkingDistance       18.000000', 
    '#', 
    '# Phase                 1', 
    '# MaterialName          Al', 
    '# Formula               Fe', 
    '# Info', 
    '# Symmetry              43', 
    '# LatticeConstants      2.870 2.870 2.870  90.000  90.000  90.000', 
    '# NumberFamilies        4', 
    '# hklFamilies           1  1  0 1 0.000000 1', 
    '# hklFamilies           2  0  0 1 0.000000 1', 
    '# hklFamilies           2  1  1 1 0.000000 1', 
    '# hklFamilies           3  1  0 1 0.000000 1', 
    '# Categories            0 0 0 0 0 ', 
    '#', 
    '# GRID: SquareGrid', 
    '# XSTEP: ' + str(step), 
    '# YSTEP: ' + str(step), 
    '# NCOLS_ODD: ' + str(sizeX), 
    '# NCOLS_EVEN: ' + str(sizeX), 
    '# NROWS: ' + str(sizeY), 
    '#', 
    '# OPERATOR: ODFsammpling', 
    '#', 
    '# SAMPLEID: ', 
    '#', 
    '# SCANID: ', 
    '#'
    ]
  else:  # ctf format
    return [ \
    'Channel Text File',
    'Prj      X:\\xxx\\xxxx.cpr',
    'Author   [MPIE-DAMASK]',
    'JobMode  Grid',
    'XCells   '+ str(sizeX),
    'YCells   '+ str(sizeY),
    'XStep    '+ str(step),
    'YStep    '+ str(step),
    'AcqE1    0',
    'AcqE2    90',
    'AcqE3    0',
    'Euler angles refer to Sample Coordinate system (CS0)! Mag 300  Coverage 100  Device 0  KV 20  TiltAngle 70  TiltAxis 0',
    'Phases   1',
    '4.05;4.05;4.05  90;90;90  Aluminium  11  225  3803863129_5.0.6.3    -2102160418    Cryogenics18,54-55',
    'Phase X      Y Bands Error Euler1     Euler2     Euler3    MAD    BC BS'
    ]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
output the Euler angles in the format supported by TSL (.ang or .ctf).

""", version = scriptID)

parser.add_option("-c", "--column", type="int", dest="column",
                  help="starting column of Euler triplet")
parser.add_option("-s", "--skip",   type="int", dest="skip",
                  help="skip this many lines of heading info [%default]")
parser.add_option("-f", "--format", type="string", dest="format",
                  help="the format of the output file [%default]")

parser.set_defaults (column = 1)
parser.set_defaults (skip = 1)
parser.set_defaults (format = 'ang')

(options,filenames) = parser.parse_args()
options.column -= 1

#--- setup file handles ---------------------------------------------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN','input':sys.stdin,'output':sys.stdout,'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,'input':open(name),'output':open(name+'_tmp','w'),'croak':sys.stdout})
    else:
      print('No such file or directory: '+name)

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  # open texture file and read content
  textureFile = open(file['name'])
  content     = textureFile.readlines()
  textureFile.close()

  m = re.match('(\d+)\s+head',content[0],re.I)
  if m != None and options.skip == 0: options.skip = int(m.group(1))+1

  # extract orientation angles
  if options.format == 'ang': # 'ang' file, radian
    angles = [map(positiveRadians,line.split()[options.column:options.column+3]) for line in content[options.skip:]]
  else:                       # 'ctf' file, degree
    angles = [line.split()[options.column:options.column+3] for line in content[options.skip:]]

  nPoints = len(angles)
  sizeY = integerFactorization(nPoints)
  sizeX = nPoints / sizeY
  file['croak'].write('%s: %i*%i = %i (== %i)\n'%(file['name'],sizeX,sizeY,sizeX*sizeY,nPoints) )

  # write ang/ctf file
  for line in getHeader(sizeX,sizeY,1.0,options.format):
    file['output'].write(line + '\n')
  for counter,point in enumerate(angles):
    if options.format == 'ang': 
      file['output'].write(''.join(['%10.5f'%angle for angle in point])+
                           ''.join(['%10.5f'%coord for coord in [counter%sizeX,counter//sizeX]])+
                           ' 100.0 1.0 0 1 1.0\n')
    else:
      file['output'].write('1  '+
                           ' '.join(['%6.1f'%coord for coord in [counter%sizeX,counter//sizeX]])+
                           '  5  0 '+
                           ' '.join(['%10.5f'%float(angle) for angle in point])+
                           '  0.5000  100  0\n')
  if file['name'] != 'STDIN':
     file['output'].close()
     os.rename(file['name']+'_tmp', os.path.splitext(file['name'])[0]+'.'+options.format )