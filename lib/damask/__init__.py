# -*- coding: UTF-8 no BOM -*-

# $Id$
import sys, os

from .environment import Environment      # only one class
from .asciitable  import ASCIItable       # only one class
from .config      import Material         # will be extended to debug and numerics
from .colormaps   import Colormap, Color
try:
    from .corientation import Quaternion, Rodrigues, Symmetry, Orientation
    print "Import Cython version of Orientation module"
except:
    from .orientation import Quaternion, Rodrigues, Symmetry, Orientation
#from .block       import Block           # only one class
from .result      import Result           # only one class
from .geometry    import Geometry         # one class with subclasses
from .solver      import Solver           # one class with subclasses
from .test        import Test
from .util        import extendableOption

try:
  from .          import core
# cleaning up namespace
###################################################################################################
# capitalize according to convention
  core.IO                            = core.io
  core.FEsolving                     = core.fesolving
  core.DAMASK_interface              = core.damask_interface
# remove modulePrefix_
  core.prec.init                     = core.prec.prec_init
  core.DAMASK_interface.init         = core.DAMASK_interface.DAMASK_interface_init
  core.IO.init                       = core.IO.IO_init
  core.numerics.init                 = core.numerics.numerics_init
  core.debug.init                    = core.debug.debug_init
  core.math.init                     = core.math.math_init
  core.math.periodicNearestNeighbor  = core.math.math_periodicNearestNeighbor
#  core.math.barycentricInterpolate   = core.math.math_barycentricInterpolate
  core.math.periodicNearestNeighborDistances  = core.math.math_periodicNearestNeighborDistances
  core.math.tensorAvg                = core.math.math_tensorAvg
  core.FEsolving.init                = core.FEsolving.FE_init
  core.mesh.init                     = core.mesh.mesh_init
  core.mesh.regrid                   = core.mesh.mesh_regrid
  core.mesh.nodesAroundCentres       = core.mesh.mesh_nodesAroundCentres
  core.mesh.deformedCoordsFFT        = core.mesh.mesh_deformedCoordsFFT
  core.mesh.volumeMismatch           = core.mesh.mesh_volumeMismatch
  core.mesh.shapeMismatch            = core.mesh.mesh_shapeMismatch

except (ImportError,AttributeError) as e:
  core = None # from http://www.python.org/dev/peps/pep-0008/
  if os.path.split(sys.argv[0])[1] not in ('symLink_Processing.py',
                                           'compile_CoreModule.py',
                                          ):
    sys.stderr.write('\nWARNING: Core module (Fortran code) not available, \n'\
                     "try to run 'make processing'\n"\
                     'Error message when importing core.so: %s\n\n'%e)
