# SPDX-License-Identifier: AGPL-3.0-or-later
"""Functionality for type hints."""

from os import PathLike
from collections.abc import Sequence
from typing import Literal, BinaryIO, TextIO, TypedDict

import numpy as np


FloatSequence = np.ndarray | Sequence[float]
IntSequence = np.ndarray | Sequence[int]
StrSequence = np.ndarray | Sequence[str]

FileHandleText = (
    TextIO
    | str
    | PathLike[str]
)

FileHandleBinary = (
    BinaryIO
    | str
    | PathLike[str]
)

CrystalFamily = Literal[
    'triclinic',
    'monoclinic',
    'orthorhombic',
    'tetragonal',
    'hexagonal',
    'cubic',
]

BravaisLattice = Literal[
    'aP',
    'mP',
    'mS',
    'oP',
    'oS',
    'oI',
    'oF',
    'tP',
    'tI',
    'hP',
    'cP',
    'cI',
    'cF',
]

CrystalKinematics = Literal[
    'slip',
    'twin',
]

# https://scientific-python.org/specs/spec-0007/
SeedLike = (
    int
    | np.integer
    | IntSequence
    | np.random.SeedSequence
)
RNGLike = (
    np.random.Generator
    | np.random.BitGenerator
)

# https://peps.python.org/pep-0655/
# Metadata = TypedDict(
#     'Metadata',
#     {
#         'unit': str,
#         'description': str,
#         'creator': str,
#         'lattice': NotRequired[str],
#     },
# )

_Metadata = TypedDict(
    '_Metadata',
    {
        'lattice': str,
        'c/a': float,
        'systems': list[str],
    },
    total=False,
)


class Metadata(_Metadata):
    unit: str
    description: str
    creator: str


DADF5Dataset = TypedDict(
    'DADF5Dataset',
    {
        'data': np.ndarray,
        'label': str,
        'meta': Metadata,
    },
)
