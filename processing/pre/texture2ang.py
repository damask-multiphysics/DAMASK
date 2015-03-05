#!/usr/bin/env python

import os,sys,math,re
from optparse import OptionParser

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


def getHeader(sizeX,sizeY,step):
  
  return [ \
  '# TEM_PIXperUM          1.000000', \
  '# x-star                0.509548', \
  '# y-star                0.795272', \
  '# z-star                0.611799', \
  '# WorkingDistance       18.000000', \
  '#', \
  '# Phase                 1', \
  '# MaterialName          Al', \
  '# Formula               Fe', \
  '# Info', \
  '# Symmetry              43', \
  '# LatticeConstants      2.870 2.870 2.870  90.000  90.000  90.000', \
  '# NumberFamilies        4', \
  '# hklFamilies           1  1  0 1 0.000000 1', \
  '# hklFamilies           2  0  0 1 0.000000 1', \
  '# hklFamilies           2  1  1 1 0.000000 1', \
  '# hklFamilies           3  1  0 1 0.000000 1', \
  '# Categories            0 0 0 0 0 ', \
  '#', \
  '# GRID: SquareGrid', \
  '# XSTEP: ' + str(step), \
  '# YSTEP: ' + str(step), \
  '# NCOLS_ODD: ' + str(sizeX), \
  '# NCOLS_EVEN: ' + str(sizeX), \
  '# NROWS: ' + str(sizeY), \
  '#', \
  '# OPERATOR: ODFsammpling', \
  '#', \
  '# SAMPLEID: ', \
  '#', \
  '# SCANID: ', \
  '#', \
  ]


parser = OptionParser(usage='%prog [options] datafile(s)')
parser.add_option("-c", "--column", type="int",\
                  dest="column",\
                  help="starting column of Euler triplet")
parser.add_option("-s", "--skip", type="int",\
                  dest="skip",\
                  help="skip this many lines of heading info [%default]")

parser.set_defaults (column = 1)
parser.set_defaults (skip = 0)

(options, files) = parser.parse_args()
options.column -= 1

if files == []:
  parser.error('no input file specified...')
  sys.exit(1)
if options.column < 0:
  parser.error('column needs to be 1,2,...')
  sys.exit(1)

while len(files) > 0:
  textureFilename = files.pop()
  baseName = os.path.splitext(textureFilename)[0]

  # open texture file and read content
  textureFile = open(textureFilename)
  content     = textureFile.readlines()
  textureFile.close()

  m = re.match('(\d+)\s+head',content[0],re.I)
  if m != None and options.skip == 0:
    options.skip = int(m.group(1))+1

  # extract orientation angles
  angles = [map(positiveRadians,line.split()[options.column:options.column+3]) for line in content[options.skip:]]

  nPoints = len(angles)
  sizeY = integerFactorization(nPoints)
  sizeX = nPoints / sizeY
  print '%s: %i * %i = %i (== %i)'%(baseName,sizeX,sizeY,sizeX*sizeY,nPoints) 
  # write ang file

  try:
    
    # write header
    angFile = open(baseName + '.ang','w')
    for line in getHeader(sizeX,sizeY,1.0):
      angFile.write(line + '\n')
    
    # write data
    counter = 0
    for point in angles:
      angFile.write(''.join(['%10.5f'%angle for angle in point])+
            ''.join(['%10.5f'%coord for coord in [counter%sizeX,counter//sizeX]])+
            ' 100.0 1.0 0 1 1.0\n')
      counter += 1
    
    angFile.close()

  except: 
    print 'unable to write',baseName
