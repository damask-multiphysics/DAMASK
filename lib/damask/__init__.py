# $Id$

from .environment import Environment      # only one class
from .asciitable  import ASCIItable       # only one class
from .config      import Material         # will be extended to debug and numerics
from .result      import Result           # one class with subclasses
from .geometry    import Geometry         # one class with subclasses
from .solver      import Solver           # one class with subclasses
from .test        import Test
try:
  from .          import core
except:
  print 'Core module (Fortran code) not available, try to run setup_processing.py'
