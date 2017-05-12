# -*- coding: UTF-8 no BOM -*-

"""Main aggregator"""
import os,sys,time

h5py_flag = os.path.join(os.path.dirname(__file__),'../../.noH5py')
h5py_grace = 7200                         # only complain once every 7200 sec (2 hours)
h5py_msg = "h5py module not found."

now = time.time()

with open(os.path.join(os.path.dirname(__file__),'../../VERSION')) as f:
  version = f.readline()[:-1]

from .environment import Environment      # noqa
from .asciitable  import ASCIItable       # noqa
try:
  from .h5table   import H5Table          # noqa
  if os.path.exists(h5py_flag): os.remove(h5py_flag)  # delete flagging file on success
except ImportError:
  if os.path.exists(h5py_flag):
    if now - os.path.getmtime(h5py_flag) > h5py_grace: # complain (again) every so-and-so often
      sys.stderr.write(h5py_msg+'\n')
      with open(h5py_flag, 'a'):
        os.utime(h5py_flag,(now,now))     # update flag modification time to "now"
  else:
    open(h5py_flag, 'a').close()          # create flagging file
    sys.stderr.write(h5py_msg+'\n')       # complain for the first time
    
from .config      import Material         # noqa
from .colormaps   import Colormap, Color  # noqa
from .orientation import Quaternion, Rodrigues, Symmetry, Orientation # noqa

#from .block       import Block           # only one class
from .result      import Result           # noqa
from .geometry    import Geometry         # noqa
from .solver      import Solver           # noqa
from .test        import Test             # noqa
from .util        import extendableOption # noqa
