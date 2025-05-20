from typing import Optional, Union, Literal
import math

import numpy as np

from ._typehints import FloatSequence, IntSequence, CrystalFamily, BravaisLattice, CrystalKinematics
from . import util
from . import Rotation


_kinematics: dict[BravaisLattice, dict[CrystalKinematics, list[np.ndarray]]] = {
    'cF': {
        'slip': [np.array([
                   [ 0,+1,-1, +1,+1,+1],
                   [-1, 0,+1, +1,+1,+1],
                   [+1,-1, 0, +1,+1,+1],
                   [ 0,-1,-1, -1,-1,+1],
                   [+1, 0,+1, -1,-1,+1],
                   [-1,+1, 0, -1,-1,+1],
                   [ 0,-1,+1, +1,-1,-1],
                   [-1, 0,-1, +1,-1,-1],
                   [+1,+1, 0, +1,-1,-1],
                   [ 0,+1,+1, -1,+1,-1],
                   [+1, 0,-1, -1,+1,-1],
                   [-1,-1, 0, -1,+1,-1]]),
                 np.array([
                   [+1,+1, 0, +1,-1, 0],
                   [+1,-1, 0, +1,+1, 0],
                   [+1, 0,+1, +1, 0,-1],
                   [+1, 0,-1, +1, 0,+1],
                   [ 0,+1,+1,  0,+1,-1],
                   [ 0,+1,-1,  0,+1,+1]])],
        'twin': [np.array([
                   [-2, 1, 1,  1, 1, 1],
                   [ 1,-2, 1,  1, 1, 1],
                   [ 1, 1,-2,  1, 1, 1],
                   [ 2,-1, 1, -1,-1, 1],
                   [-1, 2, 1, -1,-1, 1],
                   [-1,-1,-2, -1,-1, 1],
                   [-2,-1,-1,  1,-1,-1],
                   [ 1, 2,-1,  1,-1,-1],
                   [ 1,-1, 2,  1,-1,-1],
                   [ 2, 1,-1, -1, 1,-1],
                   [-1,-2,-1, -1, 1,-1],
                   [-1, 1, 2, -1, 1,-1]])]
    },
    'cI': {
        'slip': [np.array([
                   [+1,-1,+1, +0,+1,+1],
                   [+1,-1,+1, +1,+0,-1],
                   [+1,-1,+1, -1,-1,+0],
                   [-1,-1,+1, +0,-1,-1],
                   [-1,-1,+1, +1,+0,+1],
                   [-1,-1,+1, -1,+1,+0],
                   [+1,+1,+1, +0,+1,-1],
                   [+1,+1,+1, -1,+0,+1],
                   [+1,+1,+1, +1,-1,+0],
                   [-1,+1,+1, +0,-1,+1],
                   [-1,+1,+1, -1,+0,-1],
                   [-1,+1,+1, +1,+1,+0]]),
                 np.array([
                   [+1,-1,+1, +2,+1,-1],
                   [+1,-1,+1, -1,+1,+2],
                   [+1,-1,+1, +1,+2,+1],
                   [-1,-1,+1, +2,-1,+1],
                   [-1,-1,+1, +1,+1,+2],
                   [-1,-1,+1, -1,+2,+1],
                   [+1,+1,+1, +1,+1,-2],
                   [+1,+1,+1, +1,-2,+1],
                   [+1,+1,+1, -2,+1,+1],
                   [-1,+1,+1, +1,-1,+2],
                   [-1,+1,+1, +1,+2,-1],
                   [-1,+1,+1, +2,+1,+1]]),
                 np.array([
                   [+1,-1,+1, -1,+2,+3],
                   [+1,-1,+1, +1,+3,+2],
                   [+1,-1,+1, -2,+1,+3],
                   [+1,-1,+1, +2,+3,+1],
                   [+1,-1,+1, +3,+1,-2],
                   [+1,-1,+1, +3,+2,-1],
                   [-1,-1,+1, +1,+2,+3],
                   [-1,-1,+1, -1,+3,+2],
                   [-1,-1,+1, +2,+1,+3],
                   [-1,-1,+1, -2,+3,+1],
                   [-1,-1,+1, +3,-1,+2],
                   [-1,-1,+1, +3,-2,+1],
                   [+1,+1,+1, +1,+2,-3],
                   [+1,+1,+1, +1,-3,+2],
                   [+1,+1,+1, +2,+1,-3],
                   [+1,+1,+1, +2,-3,+1],
                   [+1,+1,+1, -3,+1,+2],
                   [+1,+1,+1, -3,+2,+1],
                   [-1,+1,+1, +1,-2, 3],
                   [-1,+1,+1, +1,+3,-2],
                   [-1,+1,+1, +2,-1,+3],
                   [-1,+1,+1, +2,+3,-1],
                   [-1,+1,+1, +3,+1,+2],
                   [-1,+1,+1, +3,+2,+1]])],
        'twin': [np.array([
                   [+1,-1,+1, +2,+1,-1],
                   [+1,-1,+1, -1,+1,+2],
                   [+1,-1,+1, +1,+2,+1],
                   [-1,-1,+1, +2,-1,+1],
                   [-1,-1,+1, +1,+1,+2],
                   [-1,-1,+1, -1,+2,+1],
                   [+1,+1,+1, +1,+1,-2],
                   [+1,+1,+1, +1,-2,+1],
                   [+1,+1,+1, -2,+1,+1],
                   [-1,+1,+1, +1,-1,+2],
                   [-1,+1,+1, +1,+2,-1],
                   [-1,+1,+1, +2,+1,+1]])]
    },
    'hP': {
        'slip': [np.array([
                   [+2,-1,-1, 0,  0, 0, 0,+1],
                   [-1,+2,-1, 0,  0, 0, 0,+1],
                   [-1,-1,+2, 0,  0, 0, 0,+1]]),
                 np.array([
                   [+2,-1,-1, 0,  0,+1,-1, 0],
                   [-1,+2,-1, 0, -1, 0,+1, 0],
                   [-1,-1,+2, 0, +1,-1, 0, 0]]),
                 np.array([
                   [-1,+2,-1, 0, +1, 0,-1,+1],
                   [-2,+1,+1, 0,  0,+1,-1,+1],
                   [-1,-1,+2, 0, -1,+1, 0,+1],
                   [+1,-2,+1, 0, -1, 0,+1,+1],
                   [+2,-1,-1, 0,  0,-1,+1,+1],
                   [+1,+1,-2, 0, +1,-1, 0,+1]]),
                 np.array([
                   [-2,+1,+1,+3, +1, 0,-1,+1],
                   [-1,-1,+2,+3, +1, 0,-1,+1],
                   [-1,-1,+2,+3,  0,+1,-1,+1],
                   [+1,-2,+1,+3,  0,+1,-1,+1],
                   [+1,-2,+1,+3, -1,+1, 0,+1],
                   [+2,-1,-1,+3, -1,+1, 0,+1],
                   [+2,-1,-1,+3, -1, 0,+1,+1],
                   [+1,+1,-2,+3, -1, 0,+1,+1],
                   [+1,+1,-2,+3,  0,-1,+1,+1],
                   [-1,+2,-1,+3,  0,-1,+1,+1],
                   [-1,+2,-1,+3, +1,-1, 0,+1],
                   [-2,+1,+1,+3, +1,-1, 0,+1]]),
                 np.array([
                   [-1,-1,+2,+3, +1,+1,-2,+2],
                   [+1,-2,+1,+3, -1,+2,-1,+2],
                   [+2,-1,-1,+3, -2,+1,+1,+2],
                   [+1,+1,-2,+3, -1,-1,+2,+2],
                   [-1,+2,-1,+3, +1,-2,+1,+2],
                   [-2,+1,+1,+3, +2,-1,-1,+2]])],
        'twin': [np.array([
                   [-1, 0, 1, 1,  1, 0,-1, 2],   # <-10.1>{10.2}; shear = (3-(c/a)^2)/(sqrt(3) c/a)
                   [ 0,-1, 1, 1,  0, 1,-1, 2],
                   [ 1,-1, 0, 1, -1, 1, 0, 2],
                   [ 1, 0,-1, 1, -1, 0, 1, 2],
                   [ 0, 1,-1, 1,  0,-1, 1, 2],
                   [-1, 1, 0, 1,  1,-1, 0, 2]]),
                 np.array([
                   [-1,-1, 2, 6,  1, 1,-2, 1],   # <11.6>{-1-1.1}; shear = 1/(c/a)
                   [ 1,-2, 1, 6, -1, 2,-1, 1],
                   [ 2,-1,-1, 6, -2, 1, 1, 1],
                   [ 1, 1,-2, 6, -1,-1, 2, 1],
                   [-1, 2,-1, 6,  1,-2, 1, 1],
                   [-2, 1, 1, 6,  2,-1,-1, 1]]),
                 np.array([
                   [ 1, 0,-1,-2, -1,-0, 1,-1],   # <10.-2>{10.1}; shear = (9-4(c/a)^2)/(4 sqrt(3) c/a)
                   [ 0, 1,-1,-2, -0,-1, 1,-1],
                   [-1, 1, 0,-2,  1,-1,-0,-1],
                   [-1, 0, 1,-2,  1,-0,-1,-1],
                   [ 0,-1, 1,-2, -0, 1,-1,-1],
                   [ 1,-1, 0,-2, -1, 1,-0,-1]]),
                 np.array([
                   [ 1, 1,-2,-3, -1,-1, 2,-2],   # <11.-3>{11.2}; shear = 2(2-(c/a)^2)/(3 c/a)
                   [-1, 2,-1,-3,  1,-2, 1,-2],
                   [-2, 1, 1,-3,  2,-1,-1,-2],
                   [-1,-1, 2,-3,  1, 1,-2,-2],
                   [ 1,-2, 1,-3, -1, 2,-1,-2],
                   [ 2,-1,-1,-3, -2, 1, 1,-2]])]
    },
    'tI': {
        'slip': [np.array([
                   [ 0, 0,+1, +1, 0, 0],
                   [ 0, 0,+1,  0,+1, 0]]),
                 np.array([
                   [ 0, 0,+1, +1,+1, 0],
                   [ 0, 0,+1, -1,+1, 0]]),
                 np.array([
                   [ 0,+1, 0, +1, 0, 0],
                   [+1, 0, 0,  0,+1, 0]]),
                 np.array([
                   [+1,-1,+1, +1,+1, 0],
                   [+1,-1,-1, +1,+1, 0],
                   [-1,-1,-1, -1,+1, 0],
                   [-1,-1,+1, -1,+1, 0]]),
                 np.array([
                   [+1,-1, 0, +1,+1, 0],
                   [+1,+1, 0, +1,-1, 0]]),
                 np.array([
                   [ 0,+1,+1, +1, 0, 0],
                   [ 0,-1,+1, +1, 0, 0],
                   [-1, 0,+1,  0,+1, 0],
                   [+1, 0,+1,  0,+1, 0]]),
                 np.array([
                   [ 0,+1, 0,  0, 0,+1],
                   [+1, 0, 0,  0, 0,+1]]),
                 np.array([
                   [+1,+1, 0,  0, 0,+1],
                   [-1,+1, 0,  0, 0,+1]]),
                 np.array([
                   [ 0,+1,-1,  0,+1,+1],
                   [ 0,-1,-1,  0,-1,+1],
                   [-1, 0,-1, -1, 0,+1],
                   [+1, 0,-1, +1, 0,+1]]),
                 np.array([
                   [+1,-1,+1,  0,+1,+1],
                   [+1,+1,-1,  0,+1,+1],
                   [+1,+1,+1,  0,+1,-1],
                   [-1,+1,+1,  0,+1,-1],
                   [+1,-1,-1, +1, 0,+1],
                   [-1,-1,+1, +1, 0,+1],
                   [+1,+1,+1, +1, 0,-1],
                   [+1,-1,+1, +1, 0,-1]]),
                 np.array([
                   [+1, 0, 0,  0,+1,+1],
                   [+1, 0, 0,  0,+1,-1],
                   [ 0,+1, 0, +1, 0,+1],
                   [ 0,+1, 0, +1, 0,-1]]),
                 np.array([
                   [ 0,+1,-1, +2,+1,+1],
                   [ 0,-1,-1, +2,-1,+1],
                   [+1, 0,-1, +1,+2,+1],
                   [-1, 0,-1, -1,+2,+1],
                   [ 0,+1,-1, -2,+1,+1],
                   [ 0,-1,-1, -2,-1,+1],
                   [-1, 0,-1, -1,-2,+1],
                   [+1, 0,-1, +1,-2,+1]]),
                 np.array([
                   [-1,+1,+1, +2,+1,+1],
                   [-1,-1,+1, +2,-1,+1],
                   [+1,-1,+1, +1,+2,+1],
                   [-1,-1,+1, -1,+2,+1],
                   [+1,+1,+1, -2,+1,+1],
                   [+1,-1,+1, -2,-1,+1],
                   [-1,+1,+1, -1,-2,+1],
                   [+1,+1,+1, +1,-2,+1]])]
        }
}


lattice_symmetries: dict[Optional[BravaisLattice], CrystalFamily] = {
                'aP': 'triclinic',

                'mP': 'monoclinic',
                'mS': 'monoclinic',

                'oP': 'orthorhombic',
                'oS': 'orthorhombic',
                'oI': 'orthorhombic',
                'oF': 'orthorhombic',

                'tP': 'tetragonal',
                'tI': 'tetragonal',

                'hP': 'hexagonal',

                'cP': 'cubic',
                'cI': 'cubic',
                'cF': 'cubic',
               }

orientation_relationships: dict[str, dict[str,list[np.ndarray]]] = {
  'KS': { # https://doi.org/10.1016/j.jallcom.2012.02.004
    'cF-->cI' : [
        np.repeat(np.array([
        [[-1, 0, 1],[ 1, 1, 1]],
        [[ 0, 1,-1],[ 1, 1, 1]],
        [[ 1,-1, 0],[ 1, 1, 1]],

        [[ 1, 0,-1],[ 1,-1, 1]],
        [[-1,-1, 0],[ 1,-1, 1]],
        [[ 0, 1, 1],[ 1,-1, 1]],

        [[ 0,-1, 1],[-1, 1, 1]],
        [[-1, 0,-1],[-1, 1, 1]],
        [[ 1, 1, 0],[-1, 1, 1]],

        [[-1, 1, 0],[ 1, 1,-1]],
        [[ 0,-1,-1],[ 1, 1,-1]],
        [[ 1, 0, 1],[ 1, 1,-1]],
        ]),
        2,axis=0),
        np.tile(np.array([[[-1,-1, 1],[ 0, 1, 1]],
                          [[-1, 1,-1],[ 0, 1, 1]]]),
                (12,1,1)),
    ],
    'cI-->cF' : [
        np.repeat(np.array([
        [[ 1, 1,-1],[ 0, 1, 1]],
        [[ 1,-1, 1],[ 0, 1, 1]],

        [[ 1, 1, 1],[ 0, 1,-1]],
        [[-1, 1, 1],[ 0, 1,-1]],

        [[ 1, 1,-1],[ 1, 0, 1]],
        [[ 1,-1,-1],[ 1, 0, 1]],

        [[ 1, 1, 1],[ 1, 0,-1]],
        [[ 1,-1, 1],[ 1, 0,-1]],

        [[ 1,-1, 1],[ 1, 1, 0]],
        [[ 1,-1,-1],[ 1, 1, 0]],

        [[ 1, 1, 1],[ 1,-1, 0]],
        [[ 1, 1,-1],[ 1,-1, 0]],
        ]),
        2,axis=0),
        np.tile(np.array([[[ 0, 1,-1],[ 1, 1, 1]],
                          [[ 0,-1, 1],[ 1, 1, 1]]]),
                (12,1,1)),
    ],
  },
  'GT': { # https://doi.org/10.1107/S0021889805038276
    'cF-->cI' : [
        np.array([
        [[ -5,-12, 17],[  1,  1,  1]],
        [[ 17, -5,-12],[  1,  1,  1]],
        [[-12, 17, -5],[  1,  1,  1]],
        [[  5, 12, 17],[ -1, -1,  1]],
        [[-17,  5,-12],[ -1, -1,  1]],
        [[ 12,-17, -5],[ -1, -1,  1]],
        [[ -5, 12,-17],[ -1,  1,  1]],
        [[ 17,  5, 12],[ -1,  1,  1]],
        [[-12,-17,  5],[ -1,  1,  1]],
        [[  5,-12,-17],[  1, -1,  1]],
        [[-17, -5, 12],[  1, -1,  1]],
        [[ 12, 17,  5],[  1, -1,  1]],
        [[ -5, 17,-12],[  1,  1,  1]],
        [[-12, -5, 17],[  1,  1,  1]],
        [[ 17,-12, -5],[  1,  1,  1]],
        [[  5,-17,-12],[ -1, -1,  1]],
        [[ 12,  5, 17],[ -1, -1,  1]],
        [[-17, 12, -5],[ -1, -1,  1]],
        [[ -5,-17, 12],[ -1,  1,  1]],
        [[-12,  5,-17],[ -1,  1,  1]],
        [[ 17, 12,  5],[ -1,  1,  1]],
        [[  5, 17, 12],[  1, -1,  1]],
        [[ 12, -5,-17],[  1, -1,  1]],
        [[-17,-12,  5],[  1, -1,  1]],
        ]),
        np.array([
        [[-17, -7, 17],[  1,  0,  1]],
        [[ 17,-17, -7],[  1,  1,  0]],
        [[ -7, 17,-17],[  0,  1,  1]],
        [[ 17,  7, 17],[ -1,  0,  1]],
        [[-17, 17, -7],[ -1, -1,  0]],
        [[  7,-17,-17],[  0, -1,  1]],
        [[-17,  7,-17],[ -1,  0,  1]],
        [[ 17, 17,  7],[ -1,  1,  0]],
        [[ -7,-17, 17],[  0,  1,  1]],
        [[ 17, -7,-17],[  1,  0,  1]],
        [[-17,-17,  7],[  1, -1,  0]],
        [[  7, 17, 17],[  0, -1,  1]],
        [[-17, 17, -7],[  1,  1,  0]],
        [[ -7,-17, 17],[  0,  1,  1]],
        [[ 17, -7,-17],[  1,  0,  1]],
        [[ 17,-17, -7],[ -1, -1,  0]],
        [[  7, 17, 17],[  0, -1,  1]],
        [[-17,  7,-17],[ -1,  0,  1]],
        [[-17,-17,  7],[ -1,  1,  0]],
        [[ -7, 17,-17],[  0,  1,  1]],
        [[ 17,  7, 17],[ -1,  0,  1]],
        [[ 17, 17,  7],[  1, -1,  0]],
        [[  7,-17,-17],[  0, -1,  1]],
        [[-17, -7, 17],[  1,  0,  1]],
        ]),
    ],
    'cI-->cF' : [
        np.array([
        [[-17, -7, 17],[  1,  0,  1]],
        [[ 17,-17, -7],[  1,  1,  0]],
        [[ -7, 17,-17],[  0,  1,  1]],
        [[ 17,  7, 17],[ -1,  0,  1]],
        [[-17, 17, -7],[ -1, -1,  0]],
        [[  7,-17,-17],[  0, -1,  1]],
        [[-17,  7,-17],[ -1,  0,  1]],
        [[ 17, 17,  7],[ -1,  1,  0]],
        [[ -7,-17, 17],[  0,  1,  1]],
        [[ 17, -7,-17],[  1,  0,  1]],
        [[-17,-17,  7],[  1, -1,  0]],
        [[  7, 17, 17],[  0, -1,  1]],
        [[-17, 17, -7],[  1,  1,  0]],
        [[ -7,-17, 17],[  0,  1,  1]],
        [[ 17, -7,-17],[  1,  0,  1]],
        [[ 17,-17, -7],[ -1, -1,  0]],
        [[  7, 17, 17],[  0, -1,  1]],
        [[-17,  7,-17],[ -1,  0,  1]],
        [[-17,-17,  7],[ -1,  1,  0]],
        [[ -7, 17,-17],[  0,  1,  1]],
        [[ 17,  7, 17],[ -1,  0,  1]],
        [[ 17, 17,  7],[  1, -1,  0]],
        [[  7,-17,-17],[  0, -1,  1]],
        [[-17, -7, 17],[  1,  0,  1]],
        ]),
        np.array([
        [[ -5,-12, 17],[  1,  1,  1]],
        [[ 17, -5,-12],[  1,  1,  1]],
        [[-12, 17, -5],[  1,  1,  1]],
        [[  5, 12, 17],[ -1, -1,  1]],
        [[-17,  5,-12],[ -1, -1,  1]],
        [[ 12,-17, -5],[ -1, -1,  1]],
        [[ -5, 12,-17],[ -1,  1,  1]],
        [[ 17,  5, 12],[ -1,  1,  1]],
        [[-12,-17,  5],[ -1,  1,  1]],
        [[  5,-12,-17],[  1, -1,  1]],
        [[-17, -5, 12],[  1, -1,  1]],
        [[ 12, 17,  5],[  1, -1,  1]],
        [[ -5, 17,-12],[  1,  1,  1]],
        [[-12, -5, 17],[  1,  1,  1]],
        [[ 17,-12, -5],[  1,  1,  1]],
        [[  5,-17,-12],[ -1, -1,  1]],
        [[ 12,  5, 17],[ -1, -1,  1]],
        [[-17, 12, -5],[ -1, -1,  1]],
        [[ -5,-17, 12],[ -1,  1,  1]],
        [[-12,  5,-17],[ -1,  1,  1]],
        [[ 17, 12,  5],[ -1,  1,  1]],
        [[  5, 17, 12],[  1, -1,  1]],
        [[ 12, -5,-17],[  1, -1,  1]],
        [[-17,-12,  5],[  1, -1,  1]],
        ]),
    ],
  },
  'GT_prime': { # https://doi.org/10.1107/S0021889805038276
    'cF-->cI' : [
        np.array([
        [[  0,  1, -1],[  7, 17, 17]],
        [[ -1,  0,  1],[ 17,  7, 17]],
        [[  1, -1,  0],[ 17, 17,  7]],
        [[  0, -1, -1],[ -7,-17, 17]],
        [[  1,  0,  1],[-17, -7, 17]],
        [[  1, -1,  0],[-17,-17,  7]],
        [[  0,  1, -1],[  7,-17,-17]],
        [[  1,  0,  1],[ 17, -7,-17]],
        [[ -1, -1,  0],[ 17,-17, -7]],
        [[  0, -1, -1],[ -7, 17,-17]],
        [[ -1,  0,  1],[-17,  7,-17]],
        [[ -1, -1,  0],[-17, 17, -7]],
        [[  0, -1,  1],[  7, 17, 17]],
        [[  1,  0, -1],[ 17,  7, 17]],
        [[ -1,  1,  0],[ 17, 17,  7]],
        [[  0,  1,  1],[ -7,-17, 17]],
        [[ -1,  0, -1],[-17, -7, 17]],
        [[ -1,  1,  0],[-17,-17,  7]],
        [[  0, -1,  1],[  7,-17,-17]],
        [[ -1,  0, -1],[ 17, -7,-17]],
        [[  1,  1,  0],[ 17,-17, -7]],
        [[  0,  1,  1],[ -7, 17,-17]],
        [[  1,  0, -1],[-17,  7,-17]],
        [[  1,  1,  0],[-17, 17, -7]],
        ]),
        np.array([
        [[  1,  1, -1],[ 12,  5, 17]],
        [[ -1,  1,  1],[ 17, 12,  5]],
        [[  1, -1,  1],[  5, 17, 12]],
        [[ -1, -1, -1],[-12, -5, 17]],
        [[  1, -1,  1],[-17,-12,  5]],
        [[  1, -1, -1],[ -5,-17, 12]],
        [[ -1,  1, -1],[ 12, -5,-17]],
        [[  1,  1,  1],[ 17,-12, -5]],
        [[ -1, -1,  1],[  5,-17,-12]],
        [[  1, -1, -1],[-12,  5,-17]],
        [[ -1, -1,  1],[-17, 12, -5]],
        [[ -1, -1, -1],[ -5, 17,-12]],
        [[  1, -1,  1],[ 12, 17,  5]],
        [[  1,  1, -1],[  5, 12, 17]],
        [[ -1,  1,  1],[ 17,  5, 12]],
        [[ -1,  1,  1],[-12,-17,  5]],
        [[ -1, -1, -1],[ -5,-12, 17]],
        [[ -1,  1, -1],[-17, -5, 12]],
        [[ -1, -1,  1],[ 12,-17, -5]],
        [[ -1,  1, -1],[  5,-12,-17]],
        [[  1,  1,  1],[ 17, -5,-12]],
        [[  1,  1,  1],[-12, 17, -5]],
        [[  1, -1, -1],[ -5, 12,-17]],
        [[  1,  1, -1],[-17,  5,-12]],
        ]),
    ],
    'cI-->cF' : [
        np.array([
        [[  1,  1, -1],[ 12,  5, 17]],
        [[ -1,  1,  1],[ 17, 12,  5]],
        [[  1, -1,  1],[  5, 17, 12]],
        [[ -1, -1, -1],[-12, -5, 17]],
        [[  1, -1,  1],[-17,-12,  5]],
        [[  1, -1, -1],[ -5,-17, 12]],
        [[ -1,  1, -1],[ 12, -5,-17]],
        [[  1,  1,  1],[ 17,-12, -5]],
        [[ -1, -1,  1],[  5,-17,-12]],
        [[  1, -1, -1],[-12,  5,-17]],
        [[ -1, -1,  1],[-17, 12, -5]],
        [[ -1, -1, -1],[ -5, 17,-12]],
        [[  1, -1,  1],[ 12, 17,  5]],
        [[  1,  1, -1],[  5, 12, 17]],
        [[ -1,  1,  1],[ 17,  5, 12]],
        [[ -1,  1,  1],[-12,-17,  5]],
        [[ -1, -1, -1],[ -5,-12, 17]],
        [[ -1,  1, -1],[-17, -5, 12]],
        [[ -1, -1,  1],[ 12,-17, -5]],
        [[ -1,  1, -1],[  5,-12,-17]],
        [[  1,  1,  1],[ 17, -5,-12]],
        [[  1,  1,  1],[-12, 17, -5]],
        [[  1, -1, -1],[ -5, 12,-17]],
        [[  1,  1, -1],[-17,  5,-12]],
        ]),
        np.array([
        [[  0,  1, -1],[  7, 17, 17]],
        [[ -1,  0,  1],[ 17,  7, 17]],
        [[  1, -1,  0],[ 17, 17,  7]],
        [[  0, -1, -1],[ -7,-17, 17]],
        [[  1,  0,  1],[-17, -7, 17]],
        [[  1, -1,  0],[-17,-17,  7]],
        [[  0,  1, -1],[  7,-17,-17]],
        [[  1,  0,  1],[ 17, -7,-17]],
        [[ -1, -1,  0],[ 17,-17, -7]],
        [[  0, -1, -1],[ -7, 17,-17]],
        [[ -1,  0,  1],[-17,  7,-17]],
        [[ -1, -1,  0],[-17, 17, -7]],
        [[  0, -1,  1],[  7, 17, 17]],
        [[  1,  0, -1],[ 17,  7, 17]],
        [[ -1,  1,  0],[ 17, 17,  7]],
        [[  0,  1,  1],[ -7,-17, 17]],
        [[ -1,  0, -1],[-17, -7, 17]],
        [[ -1,  1,  0],[-17,-17,  7]],
        [[  0, -1,  1],[  7,-17,-17]],
        [[ -1,  0, -1],[ 17, -7,-17]],
        [[  1,  1,  0],[ 17,-17, -7]],
        [[  0,  1,  1],[ -7, 17,-17]],
        [[  1,  0, -1],[-17,  7,-17]],
        [[  1,  1,  0],[-17, 17, -7]],
        ]),
    ],
  },
  'NW': { # https://doi.org/10.1016/j.matchar.2004.12.015
    'cF-->cI' : [
        np.array([
        [[ 2,-1,-1],[ 1, 1, 1]],
        [[-1, 2,-1],[ 1, 1, 1]],
        [[-1,-1, 2],[ 1, 1, 1]],

        [[-2,-1,-1],[-1, 1, 1]],
        [[ 1, 2,-1],[-1, 1, 1]],
        [[ 1,-1, 2],[-1, 1, 1]],

        [[ 2, 1,-1],[ 1,-1, 1]],
        [[-1,-2,-1],[ 1,-1, 1]],
        [[-1, 1, 2],[ 1,-1, 1]],

        [[ 2,-1, 1],[ 1, 1,-1]],
        [[-1, 2, 1],[ 1, 1,-1]],
        [[-1,-1,-2],[ 1, 1,-1]],
        ]),
        np.broadcast_to(np.array([[ 0,-1, 1],[ 0, 1, 1]]),
                        (12,2,3)),
    ],
    'cI-->cF' : [
        np.repeat(np.array([
            [[ 0, 1,-1],[ 0, 1, 1]],
            [[ 0, 1, 1],[ 0, 1,-1]],
            [[ 1, 0,-1],[ 1, 0, 1]],
            [[ 1, 0, 1],[ 1, 0,-1]],
            [[ 1,-1, 0],[ 1, 1, 0]],
            [[ 1, 1, 0],[ 1,-1, 0]],
            ]),
            2,axis=0),
        np.tile(np.array([
            [[ 2,-1,-1],[ 1, 1, 1]],
            [[-2, 1, 1],[ 1, 1, 1]],
            ]),
            (6,1,1)),
    ],
  },
  'Pitsch': { # https://doi.org/10.1080/14786435908238253
    'cF-->cI' : [
        np.repeat(np.array([
        [[ 0, 1, 1],[ 1, 0, 0]],
        [[ 0, 1,-1],[ 1, 0, 0]],
        [[ 1, 0, 1],[ 0, 1, 0]],
        [[ 1, 0,-1],[ 0, 1, 0]],
        [[ 1, 1, 0],[ 0, 0, 1]],
        [[ 1,-1, 0],[ 0, 0, 1]],
        ]),
        2,axis=0),
        np.tile(np.array([
        [[ 1, 1,-1],[ 0, 1, 1]],
        [[-1, 1,-1],[ 0, 1, 1]],
        ]),
        (6,1,1)),
    ],
    'cI-->cF' : [
        np.array([
        [[ 1, 1,-1],[ 0, 1, 1]],
        [[ 1,-1, 1],[ 0, 1, 1]],
        [[ 1, 1, 1],[ 0, 1,-1]],
        [[-1, 1, 1],[ 0, 1,-1]],
        [[ 1, 1,-1],[ 1, 0, 1]],
        [[ 1,-1,-1],[ 1, 0, 1]],
        [[ 1, 1, 1],[ 1, 0,-1]],
        [[ 1,-1, 1],[ 1, 0,-1]],
        [[ 1,-1, 1],[ 1, 1, 0]],
        [[ 1,-1,-1],[ 1, 1, 0]],
        [[ 1, 1, 1],[ 1,-1, 0]],
        [[ 1, 1,-1],[ 1,-1, 0]],
        ]),
        np.broadcast_to(np.array([[ 1, 1, 0],[ 0, 0, 1]]),
                        (12,2,3)),
    ],
  },
  'Bain': { # https://doi.org/10.1107/S0021889805038276
    'cF-->cI' : [
        np.array([
        [[ 0, 1, 0],[ 1, 0, 0]],
        [[ 0, 0, 1],[ 0, 1, 0]],
        [[ 1, 0, 0],[ 0, 0, 1]],
        ]),
        np.broadcast_to(np.array([[ 1, 1, 0],[ 0, 0, 1]]),
                        (3,2,3)),
    ],
    'cI-->cF' : [
        np.array([
        [[ 0, 1, 1],[ 1, 0, 0]],
        [[ 1, 0, 1],[ 0, 1, 0]],
        [[ 1, 1, 0],[ 0, 0, 1]],
        ]),
        np.broadcast_to(np.array([[ 1, 0, 0],[ 0, 0, 1]]),
                        (3,2,3)),
    ]
  },
  'Burgers' : { # https://doi.org/10.1016/S0031-8914(34)80244-3
    'cI-->hP' : [
        np.array([
        [[ 1, 1,-1],[ 0, 1, 1]],
        [[ 1,-1, 1],[ 0, 1, 1]],
        [[ 1, 1, 1],[ 0, 1,-1]],
        [[-1, 1, 1],[ 0, 1,-1]],
        [[ 1, 1,-1],[ 1, 0, 1]],
        [[ 1,-1,-1],[ 1, 0, 1]],
        [[ 1, 1, 1],[ 1, 0,-1]],
        [[ 1,-1, 1],[ 1, 0,-1]],
        [[ 1,-1, 1],[ 1, 1, 0]],
        [[ 1,-1,-1],[ 1, 1, 0]],
        [[ 1, 1, 1],[ 1,-1, 0]],
        [[ 1, 1,-1],[ 1,-1, 0]],
        ]),
        np.broadcast_to(np.array([[ 2,-1,-1, 0],[ 0, 0, 0, 1]]),
                        (12,2,4)),
    ],
    'hP-->cI' : [
        np.repeat(np.array([
        [[ 2,-1,-1, 0],[ 0, 0, 0, 1]],
        [[-1, 2,-1, 0],[ 0, 0, 0, 1]],
        [[-1,-1, 2, 0],[ 0, 0, 0, 1]],
        ]),
        2,axis=0),
        np.tile(np.array([
        [[ 1, 1,-1],[ 0, 1, 1]],
        [[-1, 1,-1],[ 0, 1, 1]],
        ]),
        (3,1,1)),
    ]
  },
}

class Crystal():
    """
    Representation of a crystal as (general) crystal family or (more specific) as a scaled Bravais lattice.

    Attributes
    ----------
    family : str
        Name of the crystal family.
    lattice : str, optional
        Name of the Bravais lattice in Pearson notation.
    a : float, optional
        Length of lattice parameter 'a'.
    b : float, optional
        Length of lattice parameter 'b'.
    c : float, optional
        Length of lattice parameter 'c'.
    alpha : float, optional
        Angle between 'b' and 'c' lattice basis.
    beta : float, optional
        Angle between 'c' and 'a' lattice basis.
    gamma : float, optional
        Angle between 'a' and 'b' lattice basis.

    Examples
    --------
    Cubic crystal family:

    >>> import damask
    >>> (cubic := damask.Crystal(family='cubic'))
    Crystal family: cubic

    Body-centered cubic Bravais lattice with parameters of iron:

    >>> import damask
    >>> (Fe := damask.Crystal(lattice='cI', a=287e-12))
    Crystal family: cubic
    Bravais lattice: cI
    a=2.87e-10 m, b=2.87e-10 m, c=2.87e-10 m
    α=90°, β=90°, γ=90°
    """

    def __init__(self, *,
                 family: Optional[CrystalFamily] = None,
                 lattice: Optional[BravaisLattice] = None,
                 a: Optional[float] = None, b: Optional[float] = None, c: Optional[float] = None,
                 alpha: Optional[float] = None, beta: Optional[float] = None, gamma: Optional[float] = None,
                 degrees: bool = False):
        """
        New representation of a crystal.

        Parameters
        ----------
        family : {'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic'}, optional
            Name of the crystal family.
            Will be inferred if 'lattice' is given.
        lattice : {'aP', 'mP', 'mS', 'oP', 'oS', 'oI', 'oF', 'tP', 'tI', 'hP', 'cP', 'cI', 'cF'}, optional
            Name of the Bravais lattice in Pearson notation.
        a : float, optional
            Length of lattice parameter 'a'.
        b : float, optional
            Length of lattice parameter 'b'.
        c : float, optional
            Length of lattice parameter 'c'.
        alpha : float, optional
            Angle between b and c lattice basis.
        beta : float, optional
            Angle between c and a lattice basis.
        gamma : float, optional
            Angle between a and b lattice basis.
        degrees : bool, optional
            Angles are given in degrees. Defaults to False.
        """
        if not lattice and not family:
            raise KeyError('Crystal initialization requires either lattice or family information')
        if family is not None and family not in list(lattice_symmetries.values()):
            raise KeyError(f'invalid crystal family "{family}"')
        if lattice is not None and family is not None and family != lattice_symmetries[lattice]:
            raise KeyError(f'incompatible family "{family}" for lattice "{lattice}"')

        self.family  = lattice_symmetries[lattice] if family is None else family
        self.lattice = lattice

        if self.lattice is not None:
            self.a = 1. if a is None else a
            self.b = b
            self.c = c
            self.a = float(self.a) if self.a is not None else \
                     (self.b / self.ratio['b'] if self.b is not None and self.ratio['b'] is not None else
                      self.c / self.ratio['c'] if self.c is not None and self.ratio['c'] is not None else None)
            self.b = float(self.b) if self.b is not None else \
                     (self.a * self.ratio['b'] if self.a is not None and self.ratio['b'] is not None else
                      self.c / self.ratio['c'] * self.ratio['b']
                      if self.c is not None and self.ratio['b'] is not None and self.ratio['c'] is not None else None)
            self.c = float(self.c) if self.c is not None else \
                     (self.a * self.ratio['c'] if self.a is not None and self.ratio['c'] is not None else
                      self.b / self.ratio['b'] * self.ratio['c']
                      if self.c is not None and self.ratio['b'] is not None and self.ratio['c'] is not None else None)

            self.alpha = math.radians(alpha) if degrees and alpha is not None else alpha
            self.beta  = math.radians(beta)  if degrees and beta  is not None else beta
            self.gamma = math.radians(gamma) if degrees and gamma is not None else gamma
            if self.alpha is None and 'alpha' in self.immutable: self.alpha = self.immutable['alpha']
            if self.beta  is None and 'beta'  in self.immutable: self.beta  = self.immutable['beta']
            if self.gamma is None and 'gamma' in self.immutable: self.gamma = self.immutable['gamma']

            if \
                (self.a     is None) \
             or (self.b     is None or ('b'     in self.immutable and self.b     != self.immutable['b'] * self.a)) \
             or (self.c     is None or ('c'     in self.immutable and self.c     != self.immutable['c'] * self.b)) \
             or (self.alpha is None or ('alpha' in self.immutable and self.alpha != self.immutable['alpha'])) \
             or (self.beta  is None or ('beta'  in self.immutable and self.beta  != self.immutable['beta'])) \
             or (self.gamma is None or ('gamma' in self.immutable and self.gamma != self.immutable['gamma'])):
                raise ValueError (f'incompatible parameters {self.parameters} for crystal family {self.family}')

            if np.any(np.array([self.alpha,self.beta,self.gamma]) <= 0):
                raise ValueError ('lattice angles must be positive')
            if np.any([np.roll([self.alpha,self.beta,self.gamma],r)[0]
              >= np.sum(np.roll([self.alpha,self.beta,self.gamma],r)[1:]) for r in range(3)]):
                raise ValueError ('each lattice angle must be less than sum of others')


    def __repr__(self):
        """
        Return repr(self).

        Give short, human-readable summary.
        """
        family = f'Crystal family: {self.family}'
        return family if self.lattice is None else \
               util.srepr([family,
                           f'Bravais lattice: {self.lattice}',
                           'a={a:.5g} m, b={b:.5g} m, c={c:.5g} m'.format(**self.parameters),
                           'α={alpha:.5g}°, β={beta:.5g}°, γ={gamma:.5g}°'
                           .format(**dict(map(lambda kv: (kv[0], np.degrees(kv[1])), self.parameters.items()))),
                           ])


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : Crystal
            Crystal to check for equality.

        Returns
        -------
        equal : bool
            Whether both arguments are equal.
        """
        return (NotImplemented if not isinstance(other, Crystal) else
                self.lattice == other.lattice and
                self.parameters == other.parameters and
                self.family == other.family)

    @property
    def parameters(self) -> Optional[dict]:
        """
        Return lattice parameters.

        Returns
        -------
        parameters : dict
            Lattice parameters a, b, c, alpha, beta, gamma.
        """
        has_parameters = all([hasattr(self,p) for p in ['a','b','c','alpha','beta','gamma']])
        return dict(a=self.a,b=self.b,c=self.c,
                    alpha=self.alpha,beta=self.beta,gamma=self.gamma) if has_parameters else None

    @property
    def immutable(self) -> dict[str, float]:
        """
        Return immutable lattice parameters.

        Returns
        -------
        immutable : dict
            Lattice parameters a, b, c, alpha, beta, gamma
            that are fixed for the given crystal family.
        """
        match self.family:
            case 'cubic':
                return {
                         'b': 1.0,
                         'c': 1.0,
                         'alpha': math.pi/2.,
                         'beta':  math.pi/2.,
                         'gamma': math.pi/2.,
                       }
            case 'hexagonal':
                return {
                         'b': 1.0,
                         'alpha': math.pi/2.,
                         'beta':  math.pi/2.,
                         'gamma': 2.*math.pi/3.,
                       }
            case 'tetragonal':
                return {
                         'b': 1.0,
                         'alpha': math.pi/2.,
                         'beta':  math.pi/2.,
                         'gamma': math.pi/2.,
                       }
            case 'orthorhombic':
                return {
                         'alpha': math.pi/2.,
                         'beta':  math.pi/2.,
                         'gamma': math.pi/2.,
                       }
            case 'monoclinic':
                return {
                         'alpha': math.pi/2.,
                         'gamma': math.pi/2.,
                       }
            case 'triclinic':
                return {}

    @property
    def orientation_relationships(self) -> list[str]:
        """
        Return labels of orientation relationships.

        Returns
        -------
        labels : list of str
            Labels of the applicable orientation relationships.
        """
        return [k for k,v in orientation_relationships.items() if np.any([m.startswith(str(self.lattice)) for m in v])]


    @property
    def standard_triangle(self) -> Union[dict[str, np.ndarray], None]:
        """
        Return corners of the standard triangle.

        Returns
        -------
        standard_triangle : dict
            Proper and improper corners of the standard triangle
            for the given crystal family.

        Notes
        -----
        Not yet defined for monoclinic.

        References
        ----------
        Bases are computed from

        >>> basis = {
        ...    'cubic' :       np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,1.]/np.sqrt(2.),               #              green
        ...                                            [1.,1.,1.]/np.sqrt(3.)]).T),          #              blue
        ...    'hexagonal' :   np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,0.],                           #              green
        ...                                            [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T), #              blue
        ...    'tetragonal' :  np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,0.],                           #              green
        ...                                            [1.,1.,0.]/np.sqrt(2.)]).T),          #              blue
        ...    'orthorhombic': np.linalg.inv(np.array([[0.,0.,1.],                           # direction of red
        ...                                            [1.,0.,0.],                           #              green
        ...                                            [0.,1.,0.]]).T),                      #              blue
        ...    }
        """
        match self.family:
            case 'cubic':
                return {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                              [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                              [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                          'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                              [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                              [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
                       }
            case 'hexagonal':
                return {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                              [ 1.            , -np.sqrt(3.)   ,  0. ],
                                              [ 0.            ,  2.            ,  0. ] ]),
                          'proper':np.array([ [ 0.            ,  0.            ,  1. ],
                                              [-1.            ,  np.sqrt(3.)   ,  0. ],
                                              [ np.sqrt(3.)   , -1.            ,  0. ] ]),
                       }
            case 'tetragonal':
                return {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                              [ 1.            , -1.            ,  0. ],
                                              [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
                          'proper':np.array([ [ 0.            ,  0.            ,  1. ],
                                              [-1.            ,  1.            ,  0. ],
                                              [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
                       }
            case 'orthorhombic':
                return {'improper':np.array([ [ 0., 0., 1.],
                                              [ 1., 0., 0.],
                                              [ 0., 1., 0.] ]),
                          'proper':np.array([ [ 0., 0., 1.],
                                              [-1., 0., 0.],
                                              [ 0., 1., 0.] ]),
                       }
            case _:
                return None


    @property
    def symmetry_operations(self) -> Rotation:
        """
        Return symmetry operations.

        Returns
        -------
        symmetry_operations : damask.Rotation
            Symmetry operations for given crystal family.

        Notes
        -----
        The symmetry operations defined here only consider Rotations.
        More specifically, for each crystal family, an enantiomorphic
        point symmetry is selected. In case that there are multiple
        point groups with enantiomorphic point symmetry, the one with
        the highest order is chosen:

        Overview of crystal classes and point group in Hermann-Mauguin
        notation used for definition of symmetry operations.
        - tricinic: 1
        - monoclinic: 2
        - orthorhombic: 222
        - tetragonal: 422
        - hexagonal: 622
        - cubic: 432

        References
        ----------
        U.F. Kocks et al.,
        Texture and Anisotropy: Preferred Orientations in Polycrystals
        and their Effect on Materials Properties.
        Cambridge University Press 1998. Table II

        https://en.wikipedia.org/wiki/Crystal_system#Crystal_classes
        """
        match self.family:
            case 'cubic': # 432
                ops = [
                        [ 1.0,            0.0,            0.0,            0.0            ],
                        [ 0.0,            1.0,            0.0,            0.0            ],
                        [ 0.0,            0.0,            1.0,            0.0            ],
                        [ 0.0,            0.0,            0.0,            1.0            ],
                        [ 0.0,            0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2) ],
                        [ 0.0,            0.0,            0.5*np.sqrt(2),-0.5*np.sqrt(2) ],
                        [ 0.0,            0.5*np.sqrt(2), 0.0,            0.5*np.sqrt(2) ],
                        [ 0.0,            0.5*np.sqrt(2), 0.0,           -0.5*np.sqrt(2) ],
                        [ 0.0,            0.5*np.sqrt(2),-0.5*np.sqrt(2), 0.0            ],
                        [ 0.0,           -0.5*np.sqrt(2),-0.5*np.sqrt(2), 0.0            ],
                        [ 0.5,            0.5,            0.5,            0.5            ],
                        [-0.5,            0.5,            0.5,            0.5            ],
                        [-0.5,            0.5,            0.5,           -0.5            ],
                        [-0.5,            0.5,           -0.5,            0.5            ],
                        [-0.5,           -0.5,            0.5,            0.5            ],
                        [-0.5,           -0.5,            0.5,           -0.5            ],
                        [-0.5,           -0.5,           -0.5,            0.5            ],
                        [-0.5,            0.5,           -0.5,           -0.5            ],
                        [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                        [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                        [-0.5*np.sqrt(2), 0.0,            0.5*np.sqrt(2), 0.0            ],
                        [-0.5*np.sqrt(2), 0.0,           -0.5*np.sqrt(2), 0.0            ],
                        [-0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0,            0.0            ],
                        [-0.5*np.sqrt(2),-0.5*np.sqrt(2), 0.0,            0.0            ],
                      ]
            case 'hexagonal': # 622
                ops = [
                        [ 1.0,            0.0,            0.0,            0.0            ],
                        [-0.5*np.sqrt(3), 0.0,            0.0,           -0.5            ],
                        [ 0.5,            0.0,            0.0,            0.5*np.sqrt(3) ],
                        [ 0.0,            0.0,            0.0,            1.0            ],
                        [-0.5,            0.0,            0.0,            0.5*np.sqrt(3) ],
                        [-0.5*np.sqrt(3), 0.0,            0.0,            0.5            ],
                        [ 0.0,            1.0,            0.0,            0.0            ],
                        [ 0.0,           -0.5*np.sqrt(3), 0.5,            0.0            ],
                        [ 0.0,            0.5,           -0.5*np.sqrt(3), 0.0            ],
                        [ 0.0,            0.0,            1.0,            0.0            ],
                        [ 0.0,           -0.5,           -0.5*np.sqrt(3), 0.0            ],
                        [ 0.0,            0.5*np.sqrt(3), 0.5,            0.0            ],
                      ]
            case 'tetragonal': # 422
                ops = [
                        [ 1.0,            0.0,            0.0,            0.0            ],
                        [ 0.0,            1.0,            0.0,            0.0            ],
                        [ 0.0,            0.0,            1.0,            0.0            ],
                        [ 0.0,            0.0,            0.0,            1.0            ],
                        [ 0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                        [ 0.0,           -0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                        [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                        [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                      ]
            case 'orthorhombic': # 222
                ops = [
                        [ 1.0,0.0,0.0,0.0 ],
                        [ 0.0,1.0,0.0,0.0 ],
                        [ 0.0,0.0,1.0,0.0 ],
                        [ 0.0,0.0,0.0,1.0 ],
                      ]
            case 'monoclinic':
                ops = [ # 2
                        [ 1.0,0.0,0.0,0.0 ],
                        [ 0.0,0.0,1.0,0.0 ],
                      ]
            case 'triclinic': # 1
                ops = [
                        [ 1.0,0.0,0.0,0.0 ],
                      ]
        return Rotation.from_quaternion(ops,accept_homomorph=True)


    @property
    def ratio(self):
        """
        Return axes ratios.

        Returns
        -------
        ratio : dict
            Ratio of lattice parameters 'b' and 'c' with respect to 'a'.
        """
        _ratio = { 'hexagonal': {'c': math.sqrt(8./3.)}}

        return dict(b = self.immutable['b']
                        if 'b' in self.immutable else
                        _ratio[self.family]['b'] if self.family in _ratio and 'b' in _ratio[self.family] else None,
                    c = self.immutable['c']
                        if 'c' in self.immutable else
                        _ratio[self.family]['c'] if self.family in _ratio and 'c' in _ratio[self.family] else None,
                   )


    @property
    def basis_real(self) -> np.ndarray:
        """
        Return orthogonal real space crystal basis.

        Returns
        -------
        basis_real : numpy.ndarray, shape(3)
            Orthogonal real space crystal basis.

        References
        ----------
        C.T. Young and J.L. Lytton, Journal of Applied Physics 43:1408–1417, 1972
        https://doi.org/10.1063/1.1661333
        """
        if (p := self.parameters) is not None:
            return np.array([
                              [1,0,0],
                              [np.cos(p['gamma']),np.sin(p['gamma']),0],
                              [np.cos(p['beta']),
                               (np.cos(p['alpha'])-np.cos(p['beta'])*np.cos(p['gamma']))                     /np.sin(p['gamma']),
                               np.sqrt(1 - np.cos(p['alpha'])**2 - np.cos(p['beta'])**2 - np.cos(p['gamma'])**2
                                     + 2 * np.cos(p['alpha'])    * np.cos(p['beta'])    * np.cos(p['gamma']))/np.sin(p['gamma'])],
                             ]).T \
                 * np.array([p['a'],p['b'],p['c']])
        else:
            raise KeyError('missing crystal lattice parameters')


    @property
    def basis_reciprocal(self) -> np.ndarray:
        """
        Return reciprocal (dual) crystal basis.

        Returns
        -------
        basis_reciprocal : numpy.ndarray, shape(3)
            Reciprocal (dual) crystal basis.
        """
        return np.linalg.inv(self.basis_real.T)


    @property
    def lattice_points(self) -> np.ndarray:                                                         # type: ignore[return]
        """
        Return lattice points.

        Returns
        -------
        lattice_points : numpy.ndarray, shape(:,3)
            Positions of atoms.
        """
        if self.lattice is None: raise KeyError('no lattice type specified')

        origin = [0.,0.,0.]
        match list(self.lattice):
            case ['h','P']:
                return np.array([origin] + [ [2./3.,1./3.,0.5] ])
            case [_,'P']:
                return np.array([origin])
            case [_,'S']:
                return np.array([origin] + [ [0.5,0.5,0.0] ])
            case [_,'I']:
                return np.array([origin] + [ [0.5,0.5,0.5] ])
            case [_,'F']:
                return np.array([origin] + [
                                             [0.0,0.5,0.5],
                                             [0.5,0.0,0.5],
                                             [0.5,0.5,0.0],
                                           ])

    def to_lattice(self, *,
                   direction: Optional[FloatSequence] = None,
                   plane: Optional[FloatSequence] = None) -> np.ndarray:                            # numpydoc ignore=PR01,PR02
        """
        Calculate lattice vector corresponding to crystal frame direction or plane normal.

        Parameters
        ----------
        direction|plane : numpy.ndarray, shape (...,3)
            Real space vector along direction or
            reciprocal space vector along plane normal.

        Returns
        -------
        Miller : numpy.ndarray, shape (...,3)
            Lattice vector of direction or plane.
            Use util.scale_to_coprime to convert to (integer) Miller indices.
        """
        if (direction is not None) ^ (plane is None):
            raise KeyError('specify either "direction" or "plane"')
        basis,axis = (self.basis_reciprocal,np.asarray(direction)) \
                     if plane is None else \
                     (self.basis_real,np.asarray(plane))
        return np.einsum('li,...l',basis,axis)


    def to_frame(self, *,
                 uvw: Optional[IntSequence] = None,
                 hkl: Optional[IntSequence] = None,
                 uvtw: Optional[IntSequence] = None,
                 hkil: Optional[IntSequence] = None) -> np.ndarray:                                 # numpydoc ignore=PR01,PR02
        """
        Calculate crystal frame vector corresponding to lattice direction [uvw]/[uvtw] or plane normal (hkl)/(hkil).

        Parameters
        ----------
        uvw|hkl|uvtw|hkil : numpy.ndarray, shape (...,3) or shape (...,4)
            Miller(–Bravais) indices of crystallographic direction or plane normal.

        Returns
        -------
        vector : numpy.ndarray, shape (...,3)
            Crystal frame vector in real space along [uvw]/[uvtw] direction or
            in reciprocal space along (hkl)/(hkil) plane normal.

        Examples
        --------
        Crystal frame vector (real space) of Magnesium corresponding to [1,1,0] direction:

        >>> import damask
        >>> Mg = damask.Crystal(lattice='hP', a=321e-12, c=521e-12)
        >>> Mg.to_frame(uvw=[1, 1, 0])
        array([1.60500000e-10, 2.77994155e-10, 0.00000000e+00])

        Crystal frame vector (reciprocal space) of Titanium along (1,0,0) plane normal:

        >>> import damask
        >>> Ti = damask.Crystal(lattice='hP', a=295e-12, c=469e-12)
        >>> Ti.to_frame(hkl=(1, 0, 0))
        array([ 3.38983051e+09,  1.95711956e+09, -4.15134508e-07])
        """
        if sum(arg is not None for arg in (uvw,hkl, uvtw,hkil)) != 1:
            raise KeyError('specify either "uvw", "hkl", "uvtw", or "hkil"')
        basis,axis = (self.basis_real,np.asarray(uvw if uvtw is None else util.Bravais_to_Miller(uvtw=uvtw))) \
                     if hkl is None and hkil is None else \
                     (self.basis_reciprocal,np.asarray(hkl if hkil is None else util.Bravais_to_Miller(hkil=hkil)))
        return np.einsum('il,...l',basis,axis)


    def kinematics(self,
                   mode: CrystalKinematics) -> dict[str, list[np.ndarray]]:
        """
        Return crystal kinematics systems.

        Parameters
        ----------
        mode : {'slip','twin'}
            Deformation mode.

        Returns
        -------
        direction_plane : dictionary
            Directions and planes of deformation mode families.

        Notes
        -----
        Kinematics of slip systems are bidirectional, i.e. the
        shape change equals the kinematics multiplied by the (signed) shear.
        In contrast, twin kinematics are unidirectional and already consider
        the distinction between extension and compression twins such that the
        shape change equals the kinematics multiplied by (absolute) twin shear.
        """
        if self.lattice is None: raise KeyError('no lattice type specified')
        master = _kinematics[self.lattice][mode]
        kinematics = {'direction':[util.Bravais_to_Miller(uvtw=m[:,0:4]) if self.lattice == 'hP'
                                                    else m[:,0:3] for m in master],
                      'plane':    [util.Bravais_to_Miller(hkil=m[:,4:8]) if self.lattice == 'hP'
                                                    else m[:,3:6] for m in master]}
        if mode == 'twin':
            gamma_char = self.characteristic_shear_twin()
            kinematics['plane'] = [np.sign(gamma_char[i]).astype(int).reshape(-1,1)*k
                                                         for i,k in enumerate(kinematics['plane'])]

        return kinematics


    def characteristic_shear_twin(self,
                                  N_twin: Union[list[int], Literal['*']] = '*') -> np.ndarray:
        """
        Return characteristic shear for twinning.

        A positive value indicates a tension twin, a negative value a
        compression twin.

        Parameters
        ----------
        N_twin : '*' or sequence of int
            Number of twin systems per twin family.
            Use '*' to select all.

        Returns
        -------
        s : numpy.ndarray, shape (...)
            Characteristic shear for twinning.

        References
        ----------
        J.W. Christian and S. Mahajan, Progress in Materials Science 39(1-2):1-157, 1995
        https://doi.org/10.1016/0079-6425(94)00007-7
        """
        if self.lattice in ['cI', 'cF']:
            N_twin_ = [len(a) for a in _kinematics[self.lattice]['twin']] if N_twin == '*' else N_twin
            return np.array([[0.5*np.sqrt(2.0)]*N_twin_[0]])
        elif self.lattice == 'hP':
            N_twin_ = [len(a) for a in _kinematics[self.lattice]['twin']] if N_twin == '*' else N_twin
            c_a = self.c/self.a                                                                     # type: ignore[operator]
            return np.array([[(3.0-c_a**2)/np.sqrt(3.0)/c_a]*N_twin_[0],
                             [1.0/c_a]*N_twin_[1],
                             [(9.0-4.0*c_a**2)/np.sqrt(48.0)/c_a]*N_twin_[2],
                             [2.0*(2.0-c_a**2)/3.0/c_a]*N_twin_[3]]
                           )
        else:
            raise KeyError(f'twin systems not defined for lattice "{self.lattice}"')


    def relation_operations(self,
                            model: str,
                            target = None) -> tuple[BravaisLattice, Rotation]:
        """
        Crystallographic orientation relationships for phase transformations.

        Parameters
        ----------
        model : str
            Name of orientation relationship.
        target : Crystal, optional
            Crystal to transform to.
            Providing this parameter allows specification of non-standard lattice parameters.
            Default is inferred from selected model and uses standard lattice parameters.

        Returns
        -------
        operations : (string, damask.Rotation)
            Resulting lattice and rotations characterizing the orientation relationship.

        References
        ----------
        S. Morito et al., Journal of Alloys and Compounds 577:s587-s592, 2013
        https://doi.org/10.1016/j.jallcom.2012.02.004

        K. Kitahara et al., Acta Materialia 54(5):1279-1288, 2006
        https://doi.org/10.1016/j.actamat.2005.11.001

        Y. He et al., Journal of Applied Crystallography 39:72-81, 2006
        https://doi.org/10.1107/S0021889805038276

        H. Kitahara et al., Materials Characterization 54(4-5):378-386, 2005
        https://doi.org/10.1016/j.matchar.2004.12.015

        Y. He et al., Acta Materialia 53(4):1179-1190, 2005
        https://doi.org/10.1016/j.actamat.2004.11.021
        """
        m_l: BravaisLattice
        o_l: BravaisLattice

        if model not in self.orientation_relationships:
            raise KeyError(f'unknown orientation relationship "{model}"')

        sep = '-->'
        search = self.lattice+sep+('' if target is None else target.lattice)                        # type: ignore[operator]
        transform = [t for t in orientation_relationships[model].keys() if t.startswith(search)]

        if len(transform) != 1:
            raise ValueError(f'invalid target lattice "{search.split(sep)[1]}"')

        m_l,o_l = transform[0].split(sep)                                                           # type: ignore[assignment]
        m_p,o_p = orientation_relationships[model][m_l+sep+o_l]
        m = Crystal(lattice=m_l) if self.parameters is None else Crystal(lattice=m_l,**self.parameters)
        o = Crystal(lattice=o_l) if target is None else target
        m_p = np.stack((m.to_frame(uvw=m_p[:,0] if m_l != 'hP' else util.Bravais_to_Miller(uvtw=m_p[:,0])),
                        m.to_frame(hkl=m_p[:,1] if m_l != 'hP' else util.Bravais_to_Miller(hkil=m_p[:,1]))),
                        axis=-2)
        o_p = np.stack((o.to_frame(uvw=o_p[:,0] if o_l != 'hP' else util.Bravais_to_Miller(uvtw=o_p[:,0])),
                        o.to_frame(hkl=o_p[:,1] if o_l != 'hP' else util.Bravais_to_Miller(hkil=o_p[:,1]))),
                        axis=-2)
        return (o_l,Rotation.from_parallel(source=m_p,target=o_p,active=True))


    def Schmid(self, *,
               N_slip: Optional[Union[IntSequence, Literal['*']]] = None,
               N_twin: Optional[Union[IntSequence, Literal['*']]] = None) -> np.ndarray:            # numpydoc ignore=PR01,PR02
        u"""
        Calculate Schmid matrix P = d ⨂ n for selected deformation systems.

        Parameters
        ----------
        N_slip|N_twin : '*' or sequence of int
            Number of deformation systems per family of the deformation system.
            Use '*' to select all.

        Returns
        -------
        P : numpy.ndarray, shape (N,3,3)
            Schmid matrix for each of the N deformation systems.

        Examples
        --------
        Schmid matrix of third octahedral slip system of a face-centered
        cubic crystal.

        >>> import numpy as np
        >>> import damask
        >>> C = damask.Crystal(lattice='cF')
        >>> np.round(C.Schmid(N_slip=[12])[3],3)
        array([[ 0.   ,  0.   , -0.   ],
               [ 0.408,  0.408, -0.408],
               [ 0.408,  0.408, -0.408]])

        Schmid matrix of the first 2nd order pyramidal <c+a> slip system of
        hexagonal crystals in dependence of the c/a ratio.

        >>> import numpy as np
        >>> import damask
        >>> Mg = damask.Crystal(lattice='hP',a=320.91e-12,c=521.03e-12)
        >>> Ti = damask.Crystal(lattice='hP',a=295.05e-12,c=468.33e-12)
        >>> np.round(Mg.Schmid(N_slip=[0,0,0,0,1]),3)
        array([[[-0.112, -0.193, -0.138],
                [-0.193, -0.335, -0.238],
                [ 0.362,  0.628,  0.447]]])
        >>> np.round(Ti.Schmid(N_slip=[0,0,0,0,1]),3)
        array([[[-0.113, -0.195, -0.142],
                [-0.195, -0.338, -0.246],
                [ 0.358,  0.62 ,  0.451]]])
        """
        if (N_slip is not None) ^ (N_twin is None):
            raise KeyError('specify either "N_slip" or "N_twin"')

        kinematics,active = (self.kinematics('slip'),N_slip) if N_twin is None else \
                            (self.kinematics('twin'),N_twin)
        everylen = list(map(len,kinematics['direction']))

        if active == '*': active = everylen
        if not active or (np.array(active) > everylen[:len(active)]).any():
            raise ValueError('Invalid number of slip/twin systems')

        d = Crystal.to_frame(self,uvw=np.vstack([kinematics['direction'][i][:n] for i,n in enumerate(active)]))
        p = Crystal.to_frame(self,hkl=np.vstack([kinematics['plane'][i][:n] for i,n in enumerate(active)]))
        return np.einsum('...i,...j',d/np.linalg.norm(d,axis=-1,keepdims=True),
                                     p/np.linalg.norm(p,axis=-1,keepdims=True))
