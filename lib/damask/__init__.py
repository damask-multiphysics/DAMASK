# -*- coding: UTF-8 no BOM -*-

"""Main aggregator"""
import os,sys,time

with open(os.path.join(os.path.dirname(__file__),'../../VERSION')) as f:
  version = f.readline()[:-1]

from .environment import Environment      # noqa
from .asciitable  import ASCIItable       # noqa
    
from .config      import Material         # noqa
from .colormaps   import Colormap, Color  # noqa
from .orientation import Quaternion, Rodrigues, Symmetry, Orientation # noqa

#from .block       import Block           # only one class
from .result      import Result           # noqa
from .geometry    import Geometry         # noqa
from .solver      import Solver           # noqa
from .test        import Test             # noqa
from .util        import extendableOption # noqa
