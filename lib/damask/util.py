# damask utility functions
import math 
try:
  import numpy
except:
  numpy=None

# Matlab like trigonometric functions that take and return angles in degrees.
for f in ['cos', 'sin', 'tan']:
    if numpy:
      exec('def %sd(deg): return (numpy.%s(numpy.deg2rad(deg)))'%(f,f))
      exec('def a%sd(val): return (numpy.rad2deg(numpy.arc%s(val)))'%(f,f))
    else:  
      exec('def %sd(deg): return (math.%s(deg/180.*math.pi))'%(f,f))
      exec('def a%sd(val): return (math.a%s(val)*180./math.pi)'%(f,f))

