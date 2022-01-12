"""Functionality for typehints."""

from typing import Sequence, Union, TextIO
from pathlib import Path

import numpy as np


FloatSequence = Union[np.ndarray,Sequence[float]]
IntSequence = Union[np.ndarray,Sequence[int]]
FileHandle = Union[TextIO, str, Path]
