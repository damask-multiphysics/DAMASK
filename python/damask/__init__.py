"""Main aggregator."""
import os

with open(os.path.join(os.path.dirname(__file__),'VERSION')) as f:
    version = f.readline()[1:-1]

name = 'damask'

# classes
from .environment import Environment      # noqa
from .asciitable  import ASCIItable       # noqa
    
from .config      import Material         # noqa
from .colormaps   import Colormap, Color  # noqa
from .orientation import Symmetry, Lattice, Rotation, Orientation # noqa
from .dadf5       import DADF5 # noqa

from .geom        import Geom             # noqa
from .solver      import Solver           # noqa
from .test        import Test             # noqa
from .util        import extendableOption # noqa

# functions in modules
from . import mechanics                   # noqa

