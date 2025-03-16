"""
Pre- and Post-Processing of DAMASK Simulations.

References
----------
D. Otto de Mentock et al., Journal of Open Source Software, 10(105):7164, 2025
https://doi.org/10.21105/joss.07164
"""

from pathlib import Path as _Path
import re as _re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

name = 'damask'
with open(_Path(__file__).parent/_Path('VERSION')) as _f:
    version = _re.sub(r'^v','',_f.readline().strip())
    __version__ = version

from .                 import _typehints       # noqa
from .                 import util             # noqa
from .                 import seeds            # noqa
from .                 import tensor           # noqa
from .                 import mechanics        # noqa
from .                 import solver           # noqa
from .                 import grid_filters     # noqa
# Modules that contain only one class (of the same name), are prefixed by a '_'.
# For example, '_colormap' contains a class called 'Colormap' which is imported as 'damask.Colormap'.
from ._rotation        import Rotation         # noqa
from ._crystal         import Crystal          # noqa
from ._orientation     import Orientation      # noqa
from ._table           import Table            # noqa
from ._colormap        import Colormap         # noqa
from ._vtk             import VTK              # noqa
from ._yaml            import YAML             # noqa
from ._configmaterial  import ConfigMaterial   # noqa
from ._loadcasegrid    import LoadcaseGrid     # noqa
from ._geomgrid        import GeomGrid         # noqa
from ._result          import Result           # noqa
