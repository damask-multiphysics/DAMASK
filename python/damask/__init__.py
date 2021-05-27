"""Tools for managing DAMASK simulations."""

from pathlib import Path as _Path
import re as _re

name = 'damask'
with open(_Path(__file__).parent/_Path('VERSION')) as _f:
    version = _re.sub(r'^v','',_f.readline().strip())
    __version__ = version

from .                 import util             # noqa
from .                 import seeds            # noqa
from .                 import tensor           # noqa
from .                 import mechanics        # noqa
from .                 import solver           # noqa
from .                 import grid_filters     # noqa
from .                 import lattice          # noqa
#Modules that contain only one class (of the same name), are prefixed by a '_'.
#For example, '_colormap' containsa class called 'Colormap' which is imported as 'damask.Colormap'.
from ._rotation        import Rotation         # noqa
from ._orientation     import Orientation      # noqa
from ._table           import Table            # noqa
from ._vtk             import VTK              # noqa
from ._colormap        import Colormap         # noqa
from ._config          import Config           # noqa
from ._configmaterial  import ConfigMaterial   # noqa
from ._grid            import Grid             # noqa
from ._result          import Result           # noqa

# deprecated
from ._test        import Test             # noqa
