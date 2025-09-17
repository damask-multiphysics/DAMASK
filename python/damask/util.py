"""Miscellaneous helper functionality."""

import sys as _sys
import datetime as _datetime
import os as _os
import subprocess as _subprocess
import shlex as _shlex
import re as _re
import signal as _signal
import fractions as _fractions
import contextlib as _contextlib
from collections import abc as _abc
from functools import reduce as _reduce, partial as _partial
from pathlib import Path as _Path
import logging
from typing import Optional as _Optional, Union as _Union, Iterable as _Iterable, \
                   Literal as _Literal, NamedTuple as _NamedTuple, \
                   Any as _Any, TextIO as _TextIO, Generator as _Generator

import numpy as _np
import h5py as _h5py

from . import version as _version
from ._typehints import FloatSequence as _FloatSequence, IntSequence as _IntSequence, \
                        NumpyRngSeed as _NumpyRngSeed, FileHandle as _FileHandle

class stdioTuple(_NamedTuple):
    stdout: str
    stderr: str


logger = logging.getLogger(__name__)

# https://svn.blender.org/svnroot/bf-blender/trunk/blender/build_files/scons/tools/bcolors.py
# https://stackoverflow.com/questions/287871
_colors = {
           'header' :   '\033[95m',
           'OK_blue':   '\033[94m',
           'OK_green':  '\033[92m',
           'warning':   '\033[93m',
           'fail':      '\033[91m',
           'end_color': '\033[0m',
           'bold':      '\033[1m',
           'dim':       '\033[2m',
           'underline': '\033[4m',
           'crossout':  '\033[9m'
          }

####################################################################################################
# Functions
####################################################################################################
def srepr(msg,
          glue: str = '\n',
          quote: bool = False) -> str:
    r"""
    Join (quoted) items with glue string.

    Parameters
    ----------
    msg : (sequence of) object with __repr__
        Items to join.
    glue : str, optional
        Glue used for joining operation. Defaults to '\n'.
    quote : bool, optional
        Quote items. Defaults to False.

    Returns
    -------
    joined : str
        String representation of the joined and quoted items.
    """
    q = '"' if quote else ''
    if (not hasattr(msg, 'strip') and
           (hasattr(msg, '__getitem__') or
            hasattr(msg, '__iter__'))):
        return glue.join(q+str(x)+q for x in msg)
    else:
        return q+(msg if isinstance(msg,str) else repr(msg))+q


def emph(msg) -> str:
    """
    Format with emphasis.

    Parameters
    ----------
    msg : (sequence of) object with __repr__
        Message to format.

    Returns
    -------
    formatted : str
        Formatted string representation of the joined items.
    """
    return _colors['bold']+srepr(msg)+_colors['end_color']

def deemph(msg) -> str:
    """
    Format with deemphasis.

    Parameters
    ----------
    msg : (sequence of) object with __repr__
        Message to format.

    Returns
    -------
    formatted : str
        Formatted string representation of the joined items.
    """
    return _colors['dim']+srepr(msg)+_colors['end_color']

def warn(msg) -> str:
    """
    Format for warning.

    Parameters
    ----------
    msg : (sequence of) object with __repr__
        Message to format.

    Returns
    -------
    formatted : str
        Formatted string representation of the joined items.
    """
    return _colors['warning']+emph(msg)+_colors['end_color']

def strikeout(msg) -> str:
    """
    Format as strikeout.

    Parameters
    ----------
    msg : (iterable of) object with __repr__
        Message to format.

    Returns
    -------
    formatted : str
        Formatted string representation of the joined items.
    """
    return _colors['crossout']+srepr(msg)+_colors['end_color']


def run(cmd: str,
        wd: str = './',
        env: _Optional[dict[str, str]] = None,
        timeout: _Optional[int] = None) -> stdioTuple:
    """
    Run a command.

    Parameters
    ----------
    cmd : str
        Command to be executed.
    wd : str, optional
        Working directory of process. Defaults to './'.
    env : dict, optional
        Environment for execution.
    timeout : int, optional
        Timeout in seconds.

    Returns
    -------
    stdout, stderr : (str, str)
        Output of the executed command.
    """
    def pass_signal(sig,_,proc,default):
        proc.send_signal(sig)
        _signal.signal(sig,default)
        _signal.raise_signal(sig)

    signals = [_signal.SIGINT,_signal.SIGTERM]

    logger.info(f"running '{cmd}' in '{wd}'")
    process = _subprocess.Popen(_shlex.split(cmd),
                                stdout = _subprocess.PIPE,
                                stderr = _subprocess.PIPE,
                                env = _os.environ if env is None else env,
                                cwd = wd,
                                encoding = 'utf-8')
    # ensure that process is terminated (https://stackoverflow.com/questions/22916783)
    sig_states = [_signal.signal(sig,_partial(pass_signal,proc=process,default=_signal.getsignal(sig))) for sig in signals]

    try:
        stdout,stderr = process.communicate(timeout=timeout)
    finally:
        for sig,state in zip(signals,sig_states):
            _signal.signal(sig,state)

    if process.returncode != 0:
        logger.error(stdout)
        logger.error(stderr)
        raise RuntimeError(f"'{cmd}' failed with returncode {process.returncode}")

    return stdioTuple(stdout, stderr)


@_contextlib.contextmanager
def open_text(fname: _FileHandle,
              mode: _Literal['r','w'] = 'r') -> _Generator[_TextIO, None, None]:                    # noqa
    """
    Open a text file with Unix line endings.

    If a path or string is given, a context manager ensures that
    the file handle is closed.
    If a file handle is given, it remains unmodified.

    Parameters
    ----------
    fname : file, str, or pathlib.Path
        Name or handle of file.
    mode : {'r','w'}, optional
        Access mode: 'r'ead or 'w'rite, defaults to 'r'.

    Returns
    -------
    f : file handle
        File handle for a text file.
    """
    if isinstance(fname, (str,_Path)):
        fhandle = open(_Path(fname).expanduser(),mode,newline=('\n' if mode == 'w' else None))
        yield fhandle
        fhandle.close()
    else:
        yield fname


def time_stamp() -> str:
    """
    Provide current time as formatted string.

    Returns
    -------
    time_stamp : str
        Current time as string in %Y-%m-%d %H:%M:%S%z format.
    """
    return _datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S%z')

def execution_stamp(class_name: str,
                    function_name: _Optional[str] = None) -> str:
    """
    Timestamp the execution of a (function within a) class.

    execution_stamp : str
        Fingerprint of an operation: Class, (function), version, and
        current time.
    """
    _function_name = '' if function_name is None else f'.{function_name}'
    return f'damask.{class_name}{_function_name} v{_version} ({time_stamp()})'


def natural_sort(key: str) -> list[_Union[int, str]]:
    """
    Natural sort.

    For use in python's 'sorted'.

    References
    ----------
    https://en.wikipedia.org/wiki/Natural_sort_order
    """
    return [int(c) if c.isdigit() else c for c in _re.split('([0-9]+)', key)]


def show_progress(iterable: _Iterable,
                  N_iter: _Optional[int] = None,
                  prefix: str = '',
                  bar_length: int = 50) -> _Any:
    """
    Decorate a loop with a progress bar.

    Use similar like enumerate.

    Parameters
    ----------
    iterable : iterable
        Iterable to be decorated.
    N_iter : int, optional
        Total number of iterations. Required if iterable is not a sequence.
    prefix : str, optional
        Prefix string. Defaults to ''.
    bar_length : int, optional
        Length of progress bar in characters. Defaults to 50.
    """
    if isinstance(iterable,_abc.Sequence):
        if N_iter is None:
            N = len(iterable)
        else:
            raise ValueError('N_iter given for sequence')
    else:
        if N_iter is None:
            raise ValueError('N_iter not given')

        N = N_iter

    if N <= 1:
        for item in iterable:
            yield item
    else:
        status = ProgressBar(N,prefix,bar_length)
        for i,item in enumerate(iterable):
            yield item
            status.update(i)


def scale_to_coprime(v: _FloatSequence,
                     N_significant: int = 9) -> _np.ndarray:
    """
    Scale vector to co-prime (relatively prime) integers.

    Parameters
    ----------
    v : sequence of float, len (:)
        Vector to scale.
    N_significant : int, optional
        Number of significant digits to consider. Defaults to 9.

    Returns
    -------
    m : numpy.ndarray, shape (:)
        Vector scaled to co-prime numbers.
    """

    def get_square_denominator(x,max_denominator):
        """Denominator of the square of a number."""
        return _fractions.Fraction(x ** 2).limit_denominator(max_denominator).denominator

    def abs_lcm(a,b):
        """Absolute value of least common multiple."""
        return _np.abs(_np.lcm(a,b))

    max_denominator = int(10**(N_significant-1))

    v_ = _np.asarray(v)
    if _np.issubdtype(v_.dtype,_np.inexact):
        v_ = _np.round(_np.asarray(v,_np.float64)/_np.max(_np.abs(v)),N_significant)
    m = (v_ * _reduce(abs_lcm, map(lambda x: int(get_square_denominator(x,max_denominator)),v_))**0.5).astype(_np.int64)
    m = m//_reduce(_np.gcd,m)

    if not _np.allclose(m/_np.max(_np.abs(m)),v/_np.max(_np.abs(v)),atol=1e-2,rtol=0):
        raise ValueError(f'invalid result "{m}" for input "{v}"')

    return m


def project_equal_angle(vector: _np.ndarray,
                        direction: _Literal['x', 'y', 'z'] = 'z',                                   # noqa
                        normalize: bool = True,
                        keepdims: bool = False) -> _np.ndarray:
    """
    Apply equal-angle projection to vector.

    Parameters
    ----------
    vector : numpy.ndarray, shape (...,3)
        Vector coordinates to be projected.
    direction : {'x', 'y', 'z'}
        Projection direction. Defaults to 'z'.
    normalize : bool
        Ensure unit length of input vector. Defaults to True.
    keepdims : bool
        Maintain three-dimensional output coordinates.
        Defaults to False.

    Returns
    -------
    coordinates : numpy.ndarray, shape (...,2 | 3)
        Projected coordinates.

    Notes
    -----
    Two-dimensional output uses right-handed frame spanned by
    the next and next-next axis relative to the projection direction,
    e.g. x-y when projecting along z and z-x when projecting along y.

    Examples
    --------
    >>> import damask
    >>> import numpy as np
    >>> project_equal_angle(np.ones(3))
    array([0.3660, 0.3660])
    >>> project_equal_angle(np.ones(3),direction='x',normalize=False,keepdims=True)
    array([0. , 0.5, 0.5])
    >>> project_equal_angle([0,1,1],direction='y',normalize=True,keepdims=False)
    array([0.4142, 0. ])
    """
    shift = 'zyx'.index(direction)
    v = _np.roll(vector/_np.linalg.norm(vector,axis=-1,keepdims=True) if normalize else vector,
                 shift,axis=-1)
    return _np.roll(_np.block([v[...,:2]/(1.0+_np.abs(v[...,2:3])),_np.zeros_like(v[...,2:3])]),
                    -shift if keepdims else 0,axis=-1)[...,:3 if keepdims else 2]

def project_equal_area(vector: _np.ndarray,
                       direction: _Literal['x', 'y', 'z'] = 'z',                                    # noqa
                       normalize: bool = True,
                       keepdims: bool = False) -> _np.ndarray:
    """
    Apply equal-area projection to vector.

    Parameters
    ----------
    vector : numpy.ndarray, shape (...,3)
        Vector coordinates to be projected.
    direction : {'x', 'y', 'z'}
        Projection direction. Defaults to 'z'.
    normalize : bool
        Ensure unit length of input vector. Defaults to True.
    keepdims : bool
        Maintain three-dimensional output coordinates.
        Defaults to False.

    Returns
    -------
    coordinates : numpy.ndarray, shape (...,2 | 3)
        Projected coordinates.

    Notes
    -----
    Two-dimensional output uses right-handed frame spanned by
    the next and next-next axis relative to the projection direction,
    e.g. x-y when projecting along z and z-x when projecting along y.

    Examples
    --------
    >>> import damask
    >>> import numpy as np
    >>> project_equal_area(np.ones(3))
    array([0.4597, 0.4597])
    >>> project_equal_area(np.ones(3),direction='x',normalize=False,keepdims=True)
    array([0. , 0.7071, 0.7071])
    >>> project_equal_area([0,1,1],direction='y',normalize=True,keepdims=False)
    array([0.5412, 0. ])
    """
    shift = 'zyx'.index(direction)
    v = _np.roll(vector/_np.linalg.norm(vector,axis=-1,keepdims=True) if normalize else vector,
                 shift,axis=-1)
    return _np.roll(_np.block([v[...,:2]/_np.sqrt(1.0+_np.abs(v[...,2:3])),_np.zeros_like(v[...,2:3])]),
                    -shift if keepdims else 0,axis=-1)[...,:3 if keepdims else 2]


def hybrid_IA(dist: _FloatSequence,
              N: int,
              rng_seed: _Optional[_NumpyRngSeed] = None) -> _np.ndarray:
    """
    Hybrid integer approximation.

    Parameters
    ----------
    dist : numpy.ndarray
        Distribution to be approximated.
    N : int
        Number of samples to draw.
    rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
        A seed to initialize the BitGenerator. Defaults to None.
        If None, then fresh, unpredictable entropy will be pulled from the OS.

    Returns
    -------
    hist : numpy.ndarray, shape (N)
        Integer approximation of the distribution.
    """
    N_opt_samples = _np.maximum(_np.count_nonzero(dist),N)                                          # random subsampling if too little samples requested
    N_inv_samples = _np.int_(0)

    scale_,scale,inc_factor = (0.0,float(N_opt_samples),1.0)
    while (not _np.isclose(scale, scale_)) and (N_inv_samples != N_opt_samples):
        repeats = _np.rint(scale*_np.array(dist)).astype(_np.int64)
        N_inv_samples = _np.sum(repeats)
        scale_,scale,inc_factor = (scale,scale+inc_factor*0.5*(scale - scale_), inc_factor*2.0) \
                                   if N_inv_samples < N_opt_samples else \
                                  (scale_,0.5*(scale_ + scale), 1.0)

    return _np.repeat(_np.arange(len(dist)),repeats)[_np.random.default_rng(rng_seed).permutation(N_inv_samples)[:N]]


def shapeshifter(fro: tuple[int, ...],
                 to: tuple[int, ...],
                 mode: _Literal['left','right'] = 'left',                                           # noqa
                 keep_ones: bool = False) -> tuple[int, ...]:
    """
    Return dimensions that reshape 'fro' to become broadcastable to 'to'.

    Parameters
    ----------
    fro : tuple
        Original shape of array.
    to : tuple
        Target shape of array after broadcasting.
        len(to) cannot be less than len(fro).
    mode : {'left', 'right'}, optional
        Indicates whether new axes are preferably added to
        either left or right of the original shape.
        Defaults to 'left'.
    keep_ones : bool, optional
        Treat '1' in fro as literal value instead of dimensional placeholder.
        Defaults to False.

    Returns
    -------
    new_dims : tuple
        Dimensions for reshape.

    Examples
    --------
    >>> import numpy as np
    >>> from damask import util
    >>> a = np.ones((3,4,2))
    >>> b = np.ones(4)
    >>> b_extended = b.reshape(util.shapeshifter(b.shape,a.shape))
    >>> (a * np.broadcast_to(b_extended,a.shape)).shape
    (3, 4, 2)
    """
    if len(fro) == 0 and len(to) == 0: return tuple()
    _fro = [1] if len(fro) == 0 else list(fro)[::-1 if mode=='left' else 1]
    _to  = [1] if len(to)  == 0 else list(to) [::-1 if mode=='left' else 1]

    final_shape: list[int] = []
    index = 0
    for i,item in enumerate(_to):
        if item == _fro[index]:
            final_shape.append(item)
            index+=1
        else:
            final_shape.append(1)
            if _fro[index] == 1 and not keep_ones:
                index+=1
        if index == len(_fro):
            final_shape = final_shape+[1]*(len(_to)-i-1)
            break
    if index != len(_fro): raise ValueError(f'shapes cannot be shifted {fro} --> {to}')
    return tuple(final_shape[::-1] if mode == 'left' else final_shape)

def shapeblender(a: tuple[int, ...],
                 b: tuple[int, ...],
                 keep_ones: bool = False) -> tuple[int, ...]:
    """
    Calculate shape that overlaps the rightmost entries of 'a' with the leftmost of 'b'.

    Parameters
    ----------
    a : tuple
        Shape of first ("left") array.
    b : tuple
        Shape of second ("right") array.
    keep_ones : bool, optional
        Treat innermost '1's as literal value instead of dimensional placeholder.
        Defaults to False.

    Returns
    -------
    shape : tuple
        Shape that overlaps the rightmost entries of 'a' with the leftmost of 'b'.

    Examples
    --------
    >>> shapeblender((3,2),(3,2))
    (3, 2)
    >>> shapeblender((4,3),(3,2))
    (4, 3, 2)
    >>> shapeblender((4,4),(3,2))
    (4, 4, 3, 2)
    >>> shapeblender((1,2),(1,2,3))
    (1, 2, 3)
    >>> shapeblender((),(2,2,1))
    (2, 2, 1)
    >>> shapeblender((1,),(2,2,1))
    (2, 2, 1)
    >>> shapeblender((1,),(2,2,1),True)
    (1, 2, 2, 1)
    """
    def is_broadcastable(a,b):
        try:
            _np.broadcast_shapes(a,b)
            return True
        except ValueError:
            return False

    a_,_b = a,b
    if keep_ones:
        i = min(len(a_),len(_b))
        while i > 0 and a_[-i:] != _b[:i]: i -= 1
        return a_ + _b[i:]
    else:
        a_ += max(0,len(_b)-len(a_))*(1,)
        while not is_broadcastable(a_,_b):
            a_ = a_ + ((1,) if len(a_)<=len(_b) else ())
            _b = ((1,) if len(_b)<len(a_) else ()) + _b
        return _np.broadcast_shapes(a_,_b)


def DREAM3D_base_group(fname: _Union[str, _Path, _h5py.File]) -> str:
    """
    Determine the base group of a DREAM.3D file.

    The base group is defined as the group (folder) that contains
    a 'SPACING' dataset in a '_SIMPL_GEOMETRY' group.

    Parameters
    ----------
    fname : str, pathlib.Path, or _h5py.File
        Filename of the DREAM.3D (HDF5) file.

    Returns
    -------
    path : str
        Path to the base group.
    """
    def get_base_group(f: _h5py.File) -> str:
        base_group = f.visit(lambda path: path.rsplit('/',2)[0] if '_SIMPL_GEOMETRY/SPACING' in path else None)
        if base_group is None:
            raise ValueError(f'could not determine base group in file "{fname}"')
        return base_group

    if isinstance(fname,_h5py.File):
        return get_base_group(fname)

    with _h5py.File(_Path(fname).expanduser(),'r') as f:
        return get_base_group(f)

def DREAM3D_cell_data_group(fname: _Union[str, _Path, _h5py.File]) -> str:
    """
    Determine the cell data group of a DREAM.3D file.

    The cell data group is defined as the group (folder) that contains
    a dataset in the base group whose length matches the total number
    of points as specified in '_SIMPL_GEOMETRY/DIMENSIONS'.

    Parameters
    ----------
    fname : str, pathlib.Path, or h5py.File
        Filename of the DREAM.3D (HDF5) file.

    Returns
    -------
    path : str
        Path to the cell data group.
    """
    def get_cell_data_group(f: _h5py.File) -> str:
        base_group = DREAM3D_base_group(f)
        cells = tuple(f['/'.join([base_group,'_SIMPL_GEOMETRY','DIMENSIONS'])][()][::-1])
        cell_data_group = f[base_group].visititems(lambda path,obj: path.split('/')[0] \
                                                   if isinstance(obj,_h5py._hl.dataset.Dataset) and _np.shape(obj)[:-1] == cells \
                                                   else None)
        if cell_data_group is None:
            raise ValueError(f'could not determine cell-data group in file "{fname}/{base_group}"')
        return cell_data_group

    if isinstance(fname,_h5py.File):
        return get_cell_data_group(fname)

    with _h5py.File(_Path(fname).expanduser(),'r') as f:
        return get_cell_data_group(f)


def _standardize_MillerBravais(idx: _IntSequence) -> _np.ndarray:
    """
    Convert Miller-Bravais indices with missing component to standard (full) form.

    Parameters
    ----------
    idx : numpy.ndarray, shape (...,4) or (...,3)
        Miller–Bravais indices of crystallographic direction [uvtw] or plane normal (hkil).
        The third index (t or i) can be omitted completely or given as "..." (Ellipsis).

    Returns
    -------
    uvtw|hkil : numpy.ndarray, shape (...,4)
        Miller-Bravais indices of [uvtw] direction or (hkil) plane normal.
    """
    def expand(v: _np.ndarray) -> _np.ndarray:
        """Expand from 3 to 4 indices."""
        return _np.block([v[...,:2], -_np.sum(v[...,:2],axis=-1,keepdims=True), v[...,2:]])

    a = _np.asarray(idx)
    if _np.issubdtype(a.dtype,_np.signedinteger):
        if a.shape[-1] == 4:
            if (_np.sum(a[...,:3],axis=-1) != 0).any(): raise ValueError(rf'u+v+t≠0 | h+k+i≠0: {a}')
            return a
        elif a.shape[-1] == 3:
            return expand(a)
    else:
        if a.shape[-1] == 4:
            b = (_np.block([a[...,:2],
                            _np.where(a[...,2:3] == ..., -_np.sum(a[...,:2],axis=-1,keepdims=True),a[...,2:3]),
                            a[...,3:]]))
            if (_np.sum(b[...,:3].astype(int),axis=-1) != 0).any(): raise ValueError(rf'u+v+t≠0 | h+k+i≠0: {b}')
        elif a.shape[-1] == 3:
            b = expand(a)

        if (b != (c := b.astype(int))).any():
            raise ValueError(f'"uvtw" | "hkil" are not (castable to) signed integers: {a}')
        return c

    raise ValueError(f'invalid Miller-Bravais indices {a}')


def Bravais_to_Miller(*,
                      uvtw: _Optional[_IntSequence] = None,
                      hkil: _Optional[_IntSequence] = None) -> _np.ndarray:                         # numpydoc ignore=PR01,PR02
    """
    Transform 4 Miller–Bravais indices to 3 Miller indices of crystal direction [uvw] or plane normal (hkl).

    Parameters
    ----------
    uvtw|hkil : numpy.ndarray, shape (...,4) or (...,3)
        Miller–Bravais indices of crystallographic direction [uvtw] or plane normal (hkil).
        The third index (t or i) can be omitted completely or given as "..." (Ellipsis).

    Returns
    -------
    uvw|hkl : numpy.ndarray, shape (...,3)
        Miller indices of [uvw] direction or (hkl) plane normal.
    """
    if (uvtw is not None) ^ (hkil is None):
        raise KeyError('specify either "uvtw" or "hkil"')
    elif uvtw is not None:
        axis,basis = _standardize_MillerBravais(uvtw),_np.array([[2,1,0,0],
                                                                 [1,2,0,0],
                                                                 [0,0,0,1]])
    elif hkil is not None:
        axis,basis = _standardize_MillerBravais(hkil),_np.array([[1,0,0,0],
                                                                 [0,1,0,0],
                                                                 [0,0,0,1]])
    uvw_hkl = _np.einsum('il,...l',basis,axis)

    return uvw_hkl//_np.gcd.reduce(uvw_hkl,axis=-1,keepdims=True)


MillerBravais_to_Miller = Bravais_to_Miller


def Miller_to_Bravais(*,
                      uvw: _Optional[_IntSequence] = None,
                      hkl: _Optional[_IntSequence] = None) -> _np.ndarray:                          # numpydoc ignore=PR01,PR02
    """
    Transform 3 Miller indices to 4 Miller–Bravais indices of crystal direction [uvtw] or plane normal (hkil).

    Parameters
    ----------
    uvw|hkl : numpy.ndarray, shape (...,3)
        Miller indices of crystallographic direction [uvw] or plane normal (hkl).

    Returns
    -------
    uvtw|hkil : numpy.ndarray, shape (...,4)
        Miller–Bravais indices of [uvtw] direction or (hkil) plane normal.
    """
    if (uvw is not None) ^ (hkl is None):
        raise KeyError('specify either "uvw" or "hkl"')
    axis,basis = (_np.asarray(uvw),_np.array([[ 2,-1, 0],
                                              [-1, 2, 0],
                                              [-1,-1, 0],
                                              [ 0, 0, 3]])) \
                 if hkl is None else \
                 (_np.asarray(hkl),_np.array([[ 1, 0, 0],
                                              [ 0, 1, 0],
                                              [-1,-1, 0],
                                              [ 0, 0, 1]]))
    if (axis != axis.astype(int)).any():
        raise ValueError(f'"uvt" | "hki" are not (castable to) signed integers: {axis}')
    uvtw_hkil = _np.einsum('il,...l',basis,axis.astype(int))

    return uvtw_hkil//_np.gcd.reduce(uvtw_hkil,axis=-1,keepdims=True)


Miller_to_MillerBravais = Miller_to_Bravais


def dict_prune(d: dict) -> dict:
    """
    Recursively remove empty dictionaries.

    Parameters
    ----------
    d : dict
        Dictionary to prune.

    Returns
    -------
    pruned : dict
        Pruned dictionary.
    """
    # https://stackoverflow.com/questions/48151953
    new = {}
    for k,v in d.items():
        if isinstance(v, dict):
            v = dict_prune(v)
        if not isinstance(v,dict) or v != {}:
            new[k] = v

    return new

def dict_flatten(d: dict) -> dict:
    """
    Recursively remove keys of single-entry dictionaries.

    Parameters
    ----------
    d : dict
        Dictionary to flatten.

    Returns
    -------
    flattened : dict
        Flattened dictionary.
    """
    if isinstance(d,dict) and len(d) == 1:
        entry = d[list(d.keys())[0]]
        new = dict_flatten(entry.copy()) if isinstance(entry,dict) else entry
    else:
        new = {k: (dict_flatten(v) if isinstance(v, dict) else v) for k,v in d.items()}

    return new


def to_list(a: _Any) -> list:
    """
    Put into list.

    Parameters
    ----------
    a : any
        Variable to put into list or convert to list.

    Returns
    -------
    l : list
        Data in list.
    """
    return [a] if not hasattr(a,'__iter__') or isinstance(a,str) else list(a)


####################################################################################################
# Classes
####################################################################################################
class ProgressBar:
    """
    Report progress of an interation as a status bar.

    Works for 0-based loops, ETA is estimated by linear extrapolation.
    """

    def __init__(self,
                 total: int,
                 prefix: str,
                 bar_length: int):
        """
        New progress bar.

        Parameters
        ----------
        total : int
            Total # of iterations.
        prefix : str
            Prefix string.
        bar_length : int
            Character length of bar.
        """
        self.total = total
        self.prefix = prefix
        self.bar_length = bar_length
        self.time_start = self.time_last_update = _datetime.datetime.now()
        self.fraction_last = 0.0

        if _sys.stdout.isatty():
            _sys.stdout.write(f"{self.prefix} {'░'*self.bar_length}   0% ETA n/a")

    def update(self,
               iteration: int) -> None:

        fraction = (iteration+1) / self.total

        if (filled_length := int(self.bar_length * fraction)) > int(self.bar_length * self.fraction_last) or \
            _datetime.datetime.now() - self.time_last_update > _datetime.timedelta(seconds=10):
            self.time_last_update = _datetime.datetime.now()
            bar = '█' * filled_length + '░' * (self.bar_length - filled_length)
            remaining_time = (_datetime.datetime.now() - self.time_start) \
                           * (self.total - (iteration+1)) / (iteration+1)
            remaining_time -= _datetime.timedelta(microseconds=remaining_time.microseconds)         # remove μs
            if _sys.stdout.isatty():
                _sys.stdout.write(f'\r{self.prefix} {bar} {fraction:>4.0%} ETA {remaining_time}')

        self.fraction_last = fraction

        if iteration == self.total - 1 and _sys.stdout.isatty():
            _sys.stdout.write('\n')
