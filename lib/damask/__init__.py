# $Id$
import sys, os

from .environment import Environment      # only one class
from .asciitable  import ASCIItable       # only one class
from .config      import Material         # will be extended to debug and numerics
from .colormaps   import Colormap, Color
from .orientation import Quaternion, Rodrigues, Symmetry, Orientation
#from .block       import Block           # only one class
from .result      import Result           # one class with subclasses
from .geometry    import Geometry         # one class with subclasses
from .solver      import Solver           # one class with subclasses
from .test        import Test

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
  core.math.curlFFT                  = core.math.math_curlFFT
  core.math.divergenceFFT            = core.math.math_divergenceFFT
  core.math.divergenceFDM            = core.math.math_divergenceFDM
  core.math.periodicNearestNeighbor  = core.math.math_periodicNearestNeighbor
  core.math.periodicNearestNeighborDistances  = core.math.math_periodicNearestNeighborDistances
  core.math.tensorAvg                = core.math.math_tensorAvg
  core.math.logstrainSpat            = core.math.math_logstrainSpat
  core.math.logstrainMat             = core.math.math_logstrainMat
  core.math.cauchy                   = core.math.math_cauchy
  core.FEsolving.init                = core.FEsolving.FE_init
  core.mesh.init                     = core.mesh.mesh_init
  core.mesh.regrid                   = core.mesh.mesh_regrid
  core.mesh.nodesAroundCentres       = core.mesh.mesh_nodesAroundCentres
  core.mesh.deformedCoordsLinear     = core.mesh.mesh_deformedCoordsLinear
  core.mesh.deformedCoordsFFT        = core.mesh.mesh_deformedCoordsFFT
  core.mesh.volumeMismatch           = core.mesh.mesh_volumeMismatch
  core.mesh.shapeMismatch            = core.mesh.mesh_shapeMismatch

except (ImportError,AttributeError) as e:
  core = None # from http://www.python.org/dev/peps/pep-0008/
  if os.path.split(sys.argv[0])[1] not in ('symLink_Processing.py',
                                           'compile_CoreModule.py',
                                          ):
    sys.stderr.write('\nWARNING: Core module (Fortran code) not available, \n'\
                     'try to run setup_processing.sh or compile_CoreModule.py\n'\
                     'Error message when importing core.so: %s\n\n'%e)
