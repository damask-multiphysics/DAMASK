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
from ._table       import Table            # noqa
from ._vtk         import VTK              # noqa
from ._colormap    import Colormap         # noqa
from ._rotation    import Rotation         # noqa
from ._lattice     import Symmetry, Lattice# noqa
from ._orientation import Orientation      # noqa
from ._result      import Result           # noqa
from ._geom        import Geom             # noqa
from .             import solver           # noqa

# deprecated
Environment = _
from ._asciitable  import ASCIItable       # noqa
from ._test        import Test             # noqa
from .config       import Material         # noqa
from .util         import extendableOption # noqa
