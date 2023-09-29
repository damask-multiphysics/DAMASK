"""Functionality for typehints."""

from typing import Sequence, Union, Literal, TextIO
from pathlib import Path

import numpy as np


FloatSequence = Union[np.ndarray,Sequence[float]]
IntSequence = Union[np.ndarray,Sequence[int]]
StrSequence = Union[np.ndarray,Sequence[str]]
FileHandle = Union[TextIO, str, Path]
CrystalFamily = Union[None,Literal['triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic']]
BravaisLattice = Union[None,Literal['aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF']]
CrystalKinematics = Literal['slip', 'twin']
NumpyRngSeed = Union[int, IntSequence, np.random.SeedSequence, np.random.Generator]
# BitGenerator does not exists in older numpy versions
#NumpyRngSeed = Union[int, IntSequence, np.random.SeedSequence, np.random.BitGenerator, np.random.Generator]
