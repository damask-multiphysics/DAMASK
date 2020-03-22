"""Main aggregator."""
import os as _os
import re as _re

name = 'damask'
with open(_os.path.join(_os.path.dirname(__file__),'VERSION')) as _f:
    version = _re.sub(r'^v','',_f.readline().strip())

# classes
from ._environment import Environment      # noqa
from ._table       import Table            # noqa
from ._vtk         import VTK              # noqa
from ._colormaps   import Colormap, Color  # noqa
from ._rotation    import Rotation         # noqa
from ._lattice     import Symmetry, Lattice# noqa
from ._orientation import Orientation      # noqa
from ._result      import Result           # noqa
from ._geom        import Geom             # noqa
from .             import solver           # noqa

# deprecated
from ._asciitable  import ASCIItable       # noqa
from ._test        import Test             # noqa
from .config       import Material         # noqa
from .util         import extendableOption # noqa
