"""Main aggregator."""
import os
import re

name = 'damask'
with open(os.path.join(os.path.dirname(__file__),'VERSION')) as f:
    version = re.sub(r'^v','',f.readline().strip())

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

# compatibility hack
from .result      import Result as DADF5  # noqa

# deprecated
from .asciitable  import ASCIItable       # noqa
from .util        import extendableOption # noqa
from .config      import Material         # noqa
from .test        import Test             # noqa

# functions in modules
from .            import mechanics        # noqa
from .            import grid_filters     # noqa

# clean temporary variables
del os
del re
del f
