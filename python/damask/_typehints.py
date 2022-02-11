"""Functionality for typehints."""

from typing import Sequence, Union, TextIO
from pathlib import Path

import numpy as np


FloatSequence = Union[np.ndarray,Sequence[float]]
IntSequence = Union[np.ndarray,Sequence[int]]
FileHandle = Union[TextIO, str, Path]
NumpyRngSeed = Union[int, IntSequence, np.random.SeedSequence, np.random.Generator]
# BitGenerator does not exists in older numpy versions
#NumpyRngSeed = Union[int, IntSequence, np.random.SeedSequence, np.random.BitGenerator, np.random.Generator]
