"""Main aggregator."""
import os

name = 'damask'
with open(os.path.join(os.path.dirname(__file__),'VERSION')) as f:
    version = f.readline().strip()

# classes
from .table       import Table            # noqa
from .asciitable  import ASCIItable       # noqa
    
from .config      import Material         # noqa
from .colormaps   import Colormap, Color  # noqa
from .orientation import Symmetry, Lattice, Rotation, Orientation # noqa
from .dadf5       import DADF5 # noqa

from .geom        import Geom             # noqa
from .test        import Test             # noqa
from .util        import extendableOption # noqa

# functions in modules
from .            import mechanics        # noqa

# clean temporary variables
del os
del f
