import numpy as np
import re

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


print point
print line
print 'll'
print np.size(lineloop)
print lineloop
print 'ss'
print surface

#for sn, ll in enumerate(lineloop):
#  for l in ll:
#    for points in l[abs(int(l))]:
#      if point[int(points)-1][0] == float(0):
#       print 'abc' 
