#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import numpy as np
import re
def all_same(items,a):
  return all(x == a for x in items)
def func(seq):
  for x in seq:
    try:
      yield float(x)
    except ValueError:
      yield x

my_geofile   = 'polyXtal_5grains.geo'
PointCount                = 0
LineCount                 = 0
LineLoopCount             = 0
PlaneSurfaceCount         = 0
SurfaceLoopCount          = 0
point    = []
line     = []
lineloop = []
plane    = []
surface  = []

f = open(my_geofile,'r')
lines = f.readlines()
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

  elif eachline.startswith('Plane', 0, len(eachline)):
    PlaneSurfaceCount +=  1
    r = re.compile('{(.*?)}')
    m = r.search(eachline)
    if m:
      a = m.group(1)
      plane.append(list(func(a.split(","))))

  elif eachline.startswith('Surface', 0,len(eachline)):
    SurfaceLoopCount +=  1
    r = re.compile('{(.*?)}')
    m = r.search(eachline)
    if m:
      a = m.group(1)
      surface.append(list(func(a.split(","))))
x_coord = []
y_coord = []
z_coord = []
xp = []
xm = []
yp = []
ym = []
zp = []
zm = []

for i,l in enumerate(lineloop):
  for lines in l:
    for pts in line[int(abs(lines)-1)]:
        x_coord.append(point[int(pts)-1][0])
        y_coord.append(point[int(pts)-1][1])
        z_coord.append(point[int(pts)-1][2])

  if all_same(x_coord,1.):
    xp.append(int(i+1))
  elif all_same(x_coord,0.):
    xm.append(int(i+1))
  elif all_same(y_coord,1.):
    yp.append(int(i+1))
  elif all_same(y_coord,0.):
    ym.append(int(i+1))
  elif all_same(z_coord,1.):
    zp.append(int(i+1))
  elif all_same(z_coord,0.):
    zm.append(int(i+1))
  x_coord = []
  y_coord = []
  z_coord = []

print 'suraces on x + are ', xp
print 'suraces on x - are ', xm
print 'suraces on y + are ', yp
print 'suraces on y - are ', ym
print 'suraces on z + are ', zp
print 'suraces on z - are ', zm


