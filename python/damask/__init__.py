"""Main aggregator."""
import os as _os
import re as _re

name = 'damask'
with open(_os.path.join(_os.path.dirname(__file__),'VERSION')) as _f:
    version = _re.sub(r'^v','',_f.readline().strip())

# classes
from .environment import Environment      # noqa
from .table       import Table            # noqa
from .ktv         import VTK              # noqa
from .colormaps   import Colormap, Color  # noqa
from .rotation    import Rotation         # noqa
from .lattice     import Symmetry, Lattice# noqa
from .orientation import Orientation      # noqa
from .result      import Result           # noqa
from .geom        import Geom             # noqa
from .solver      import Solver           # noqa

# deprecated
from .asciitable  import ASCIItable       # noqa
from .util        import extendableOption # noqa
from .config      import Material         # noqa
from .test        import Test             # noqa
