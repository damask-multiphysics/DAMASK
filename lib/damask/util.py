# some utility functions
import math 
try:
  import numpy
except:
  numpy=None

  




for f in ['cos', 'sin', 'tan']:
    exec('def %sd(deg): return ((numpy or math).%s(deg/180.*math.pi))'%(f,f))
    if numpy:
      exec('def a%sd(val): return (numpy.arc%s(val)*180./math.pi)'%(f,f))
    else:  
      exec('def a%sd(val): return (math.a%s(val)*180./math.pi)'%(f,f))
# Matlab like functions to allow for cosd(degree),...
# it is open how these compare to the numpy.cos, which can operate on arrays
# ->full array compatibility can be achieved through numpy.deg2rad|rad2deg

