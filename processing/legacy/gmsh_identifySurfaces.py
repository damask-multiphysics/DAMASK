#!/usr/bin/env python3

import os
import re
from optparse import OptionParser

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def all_same(items,a):
  return all(x == a for x in items)


def func(seq):
  for x in seq:
    try:
      yield float(x)
    except ValueError:
      yield x


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]',
                      description =' Recognize bounding surfaces and append them as physical sufaces in the geo file.',
                      version = scriptID)


parser.add_option('-n','--numvol', dest = 'N', type='int', metavar='int',
                                   help='number of physical volumes' )
parser.add_option('-f','--faces',  dest = 'surfaces', action = 'extend', type = 'string', metavar = '<string LIST>',
                                   help = 'surfaces to tag {x, y, z}')
parser.add_option('-s','--size',   dest = 'size', type='float', metavar='float',
                                   help='mesh size [%default]' )
parser.add_option('-d','--dim',    dest = 'dimension', type='int', metavar='int',
                                   help='dimension of geometry [%default]' )

(options, filename) = parser.parse_args()
parser.set_defaults(size = 0.1,
                    dimension = 3)

my_geofile   = filename[0]
numVol       = options.N

PointCount        = 0
LineCount         = 0
LineLoopCount     = 0
point    = []
line     = []
lineloop = []

f = open(my_geofile,'r')
lines = f.readlines()
f.close()
for eachline in lines:
  if eachline.startswith('Point', 0, 5):
    PointCount +=  1
    r = re.compile('{(.*?)}')
    m = r.search(eachline)
    if m:
      a = m.group(1)
      point.append(list(func(a.split(","))))
  
  elif eachline.startswith('Line (', 0, 6):
    LineCount +=  1
    r = re.compile('{(.*?)}')
    m = r.search(eachline)
    if m:
      a = m.group(1)
      line.append(list(func(a.split(","))))

  elif eachline.startswith('Line Loop', 0, 9):
    LineLoopCount +=  1
    r = re.compile('{(.*?)}')
    m = r.search(eachline)
    if m:
      a = m.group(1)
      lineloop.append(list(func(a.split(","))))

x_coord = []; y_coord = []; z_coord = []
xp = []; xm = []
yp = []; ym = []
zp = []; zm = []
xmin = min([x[0] for x in point]); xmax = max([x[0] for x in point])
ymin = min([x[1] for x in point]); ymax = max([x[1] for x in point])
zmin = min([x[2] for x in point]); zmax = max([x[2] for x in point])


if (options.dimension == 3):
  for i,l in enumerate(lineloop):
    for lines in l:
      for pts in line[int(abs(lines)-1)]:
          x_coord.append(point[int(pts)-1][0])
          y_coord.append(point[int(pts)-1][1])
          z_coord.append(point[int(pts)-1][2])
    if   all_same(x_coord,xmax) and any([surface == 'x' for surface in options.surfaces]):
      xp.append(int(i+1))
    elif all_same(x_coord,xmin) and any([surface == 'x' for surface in options.surfaces]):
      xm.append(int(i+1))
    elif all_same(y_coord,ymax) and any([surface == 'y' for surface in options.surfaces]):
      yp.append(int(i+1))
    elif all_same(y_coord,ymin) and any([surface == 'y' for surface in options.surfaces]):
      ym.append(int(i+1))
    elif all_same(z_coord,zmax) and any([surface == 'z' for surface in options.surfaces]):
      zp.append(int(i+1))
    elif all_same(z_coord,zmin) and any([surface == 'z' for surface in options.surfaces]):
      zm.append(int(i+1))
    x_coord = []
    y_coord = []
    z_coord = []
  
  with open(my_geofile,'a') as f:
    f.write('Delete Physicals; \n')
    if any([surface == 'x' for surface in options.surfaces]):
      f.write('%s%s%s\n' %('Physical Surface(1) = {',','.join(map(str, xp)),'};'))
      f.write('%s%s%s\n' %('Physical Surface(2) = {',','.join(map(str, xm)),'};'))
  
    if any([surface == 'y' for surface in options.surfaces]):
      f.write('%s%s%s\n' %('Physical Surface(3) = {',','.join(map(str, yp)),'};'))
      f.write('%s%s%s\n' %('Physical Surface(4) = {',','.join(map(str, ym)),'};'))
  
    if any([surface == 'z' for surface in options.surfaces]):
      f.write('%s%s%s\n' %('Physical Surface(5) = {',','.join(map(str, zp)),'};'))
      f.write('%s%s%s\n' %('Physical Surface(6) = {',','.join(map(str, zm)),'};'))
  
    for i in range(numVol):
      f.write('%s%d%s%d%s\n' %('Physical Volume (', i+1,') = {',i+1,'};'))
      
    f.write('Field[1] = Box;\n')  
    f.write('%s%f%s\n' %('Field[1].VIn = ', options.size,';')) 
    f.write('%s%f%s\n' %('Field[1].VOut = ',options.size,';')) 
    f.write('%s%f%s\n' %('Field[1].XMin = ',xmin,';')) 
    f.write('%s%f%s\n' %('Field[1].XMax = ',xmax,';'))  
    f.write('%s%f%s\n' %('Field[1].YMin = ',ymin,';')) 
    f.write('%s%f%s\n' %('Field[1].YMax = ',ymax,';'))  
    f.write('%s%f%s\n' %('Field[1].ZMin = ',zmin,';')) 
    f.write('%s%f%s\n' %('Field[1].ZMax = ',zmax,';'))  
    f.write('Background Field = 1;\n')     
  
  f.close()

elif (options.dimension == 2):
  for i,l in enumerate(line):
    #  for pts in line[int(abs(lines)-1)]:
    for pts in l:
      x_coord.append(point[int(pts)-1][0])
      y_coord.append(point[int(pts)-1][1])
      z_coord.append(point[int(pts)-1][2])
 
    if   all_same(x_coord,xmax) and any([surface == 'x' for surface in options.surfaces]):
      xp.append(int(i+1))
    elif all_same(x_coord,xmin) and any([surface == 'x' for surface in options.surfaces]):
      xm.append(int(i+1))
    elif all_same(y_coord,ymax) and any([surface == 'y' for surface in options.surfaces]):
      yp.append(int(i+1))
    elif all_same(y_coord,ymin) and any([surface == 'y' for surface in options.surfaces]):
      ym.append(int(i+1))
    elif all_same(z_coord,zmax) and any([surface == 'z' for surface in options.surfaces]):
      zp.append(int(i+1))
    elif all_same(z_coord,zmin) and any([surface == 'z' for surface in options.surfaces]):
      zm.append(int(i+1))
    x_coord = []
    y_coord = []
    z_coord = []


  with open(my_geofile,'a') as f:
    f.write('Delete Physicals; \n')
    if any([surface == 'x' for surface in options.surfaces]):
      f.write('%s%s%s\n' %('Physical Line(1) = {',','.join(map(str, xp)),'};'))
      f.write('%s%s%s\n' %('Physical Line(2) = {',','.join(map(str, xm)),'};'))
  
    if any([surface == 'y' for surface in options.surfaces]):
      f.write('%s%s%s\n' %('Physical Line(3) = {',','.join(map(str, yp)),'};'))
      f.write('%s%s%s\n' %('Physical Line(4) = {',','.join(map(str, ym)),'};'))
  
    if any([surface == 'z' for surface in options.surfaces]):
      f.write('%s%s%s\n' %('Physical Line(5) = {',','.join(map(str, zp)),'};'))
      f.write('%s%s%s\n' %('Physical Line(6) = {',','.join(map(str, zm)),'};'))
  
    for i in range(numVol):
      f.write('%s%d%s%d%s\n' %('Physical Surface (', i+1,') = {',i+1,'};'))
      
    f.write('Field[1] = Box;\n')  
    f.write('%s%f%s\n' %('Field[1].VIn = ', options.size,';')) 
    f.write('%s%f%s\n' %('Field[1].VOut = ',options.size,';')) 
    f.write('%s%f%s\n' %('Field[1].XMin = ',xmin,';')) 
    f.write('%s%f%s\n' %('Field[1].XMax = ',xmax,';'))  
    f.write('%s%f%s\n' %('Field[1].YMin = ',ymin,';')) 
    f.write('%s%f%s\n' %('Field[1].YMax = ',ymax,';'))  
    f.write('%s%f%s\n' %('Field[1].ZMin = ',zmin,';')) 
    f.write('%s%f%s\n' %('Field[1].ZMax = ',zmax,';'))  
    f.write('Background Field = 1;\n')     
 
  f.close()
