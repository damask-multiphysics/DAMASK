"""Tools for pre and post processing of DAMASK simulations."""
from pathlib import Path as _Path
import re as _re

name = 'damask'
with open(_Path(__file__).parent/_Path('VERSION')) as _f:
    version = _re.sub(r'^v','',_f.readline().strip())
    __version__ = version

# make classes directly accessible as damask.Class
from ._environment import Environment as _ # noqa
environment = _()
from .                 import util             # noqa
from .                 import seeds            # noqa
from .                 import mechanics        # noqa
from .                 import solver           # noqa
from .                 import grid_filters     # noqa
from .                 import lattice          # noqa
from ._rotation        import Rotation         # noqa
from ._orientation     import Orientation      # noqa
from ._table           import Table            # noqa
from ._vtk             import VTK              # noqa
from ._colormap        import Colormap         # noqa
from ._config          import Config           # noqa
from ._configmaterial  import ConfigMaterial   # noqa
from ._geom            import Geom             # noqa
from ._result          import Result           # noqa



# deprecated
Environment = _
from ._asciitable  import ASCIItable       # noqa
from ._test        import Test             # noqa
from .util         import extendableOption # noqa
