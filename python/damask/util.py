"""Miscellaneous helper functionality."""

import sys as _sys
import datetime as _datetime
import os as _os
import subprocess as _subprocess
import shlex as _shlex
import re as _re
import signal as _signal
import fractions as _fractions
from collections import abc as _abc
from functools import reduce as _reduce, partial as _partial
from typing import Callable as _Callable, Union as _Union, Iterable as _Iterable, Sequence as _Sequence, Dict as _Dict,  \
                   List as _List, Tuple as _Tuple, Literal as _Literal, Any as _Any, Collection as _Collection, TextIO as _TextIO
from pathlib import Path as _Path

import numpy as _np
import h5py as _h5py

from . import version as _version
from ._typehints import FloatSequence as _FloatSequence, NumpyRngSeed as _NumpyRngSeed, IntCollection as _IntCollection, \
                        FileHandle as _FileHandle

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
          glue: str = '\n') -> str:
    r"""
    Join items with glue string.

    Parameters
    ----------
    msg : object with __repr__ or sequence of objects with __repr__
        Items to join.
    glue : str, optional
        Glue used for joining operation. Defaults to '\n'.

    Returns
    -------
    joined : str
        String representation of the joined items.

    """
    if (not hasattr(msg, 'strip') and
           (hasattr(msg, '__getitem__') or
            hasattr(msg, '__iter__'))):
        return glue.join(str(x) for x in msg)
    else:
        return msg if isinstance(msg,str) else repr(msg)


def emph(msg) -> str:
    """
    Format with emphasis.

    Parameters
    ----------
    msg : object with __repr__ or sequence of objects with __repr__
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
    msg : object with __repr__ or sequence of objects with __repr__
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
    msg : object with __repr__ or sequence of objects with __repr__
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
    msg : object with __repr__ or iterable of objects with __repr__
        Message to format.

    Returns
    -------
    formatted : str
        Formatted string representation of the joined items.

    """
    return _colors['crossout']+srepr(msg)+_colors['end_color']


def run(cmd: str,
        wd: str = './',
        env: _Dict[str, str] = None,
        timeout: int = None) -> _Tuple[str, str]:
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
    timeout : integer, optional
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

    print(f"running '{cmd}' in '{wd}'")
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
        print(stdout)
        print(stderr)
        raise RuntimeError(f"'{cmd}' failed with returncode {process.returncode}")

    return stdout, stderr


def open_text(fname: _FileHandle,
              mode: _Literal['r','w'] = 'r') -> _TextIO:                                            # noqa
    """
    Open a text file.

    Parameters
    ----------
    fname : file, str, or pathlib.Path
        Name or handle of file.
    mode: {'r','w'}, optional
        Access mode: 'r'ead or 'w'rite, defaults to 'r'.

    Returns
    -------
    f : file handle

    """
    return fname if not isinstance(fname, (str,_Path)) else \
           open(_Path(fname).expanduser(),mode,newline=('\n' if mode == 'w' else None))


def natural_sort(key: str) -> _List[_Union[int, str]]:
    """
    Natural sort.

    For use in python's 'sorted'.

    References
    ----------
    https://en.wikipedia.org/wiki/Natural_sort_order

    """
    convert = lambda text: int(text) if text.isdigit() else text
    return [ convert(c) for c in _re.split('([0-9]+)', key) ]


def show_progress(iterable: _Iterable,
                  N_iter: int = None,
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
        Prefix string.
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


def scale_to_coprime(v: _FloatSequence) -> _np.ndarray:
    """
    Scale vector to co-prime (relatively prime) integers.

    Parameters
    ----------
    v : sequence of float, len (:)
        Vector to scale.

    Returns
    -------
    m : numpy.ndarray, shape (:)
        Vector scaled to co-prime numbers.

    """
    MAX_DENOMINATOR = 1000000

    def get_square_denominator(x):
        """Denominator of the square of a number."""
        return _fractions.Fraction(x ** 2).limit_denominator(MAX_DENOMINATOR).denominator

    def lcm(a,b):
        """Least common multiple."""
        try:
            return _np.lcm(a,b)                                                                      # numpy > 1.18
        except AttributeError:
            return a * b // _np.gcd(a, b)

    v_ = _np.array(v)
    m = (v_ * _reduce(lcm, map(lambda x: int(get_square_denominator(x)),v_))**0.5).astype(_np.int64)
    m = m//_reduce(_np.gcd,m)

    with _np.errstate(invalid='ignore'):
        if not _np.allclose(_np.ma.masked_invalid(v_/m),v_[_np.argmax(abs(v_))]/m[_np.argmax(abs(v_))]):
            raise ValueError(f'invalid result "{m}" for input "{v_}"')

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
        [0.3660254, 0.3660254]
    >>> project_equal_angle(np.ones(3),direction='x',normalize=False,keepdims=True)
        [0, 0.5, 0.5]
    >>> project_equal_angle([0,1,1],direction='y',normalize=True,keepdims=False)
        [0.41421356, 0]

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
        [0.45970084, 0.45970084]
    >>> project_equal_area(np.ones(3),direction='x',normalize=False,keepdims=True)
        [0.0, 0.70710678, 0.70710678]
    >>> project_equal_area([0,1,1],direction='y',normalize=True,keepdims=False)
        [0.5411961, 0.0]

    """
    shift = 'zyx'.index(direction)
    v = _np.roll(vector/_np.linalg.norm(vector,axis=-1,keepdims=True) if normalize else vector,
                 shift,axis=-1)
    return _np.roll(_np.block([v[...,:2]/_np.sqrt(1.0+_np.abs(v[...,2:3])),_np.zeros_like(v[...,2:3])]),
                    -shift if keepdims else 0,axis=-1)[...,:3 if keepdims else 2]

def execution_stamp(class_name: str,
                    function_name: str = None) -> str:
    """Timestamp the execution of a (function within a) class."""
    now = _datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S%z')
    _function_name = '' if function_name is None else f'.{function_name}'
    return f'damask.{class_name}{_function_name} v{_version} ({now})'


def hybrid_IA(dist: _np.ndarray,
              N: int,
              rng_seed: _NumpyRngSeed = None) -> _np.ndarray:
    """
    Hybrid integer approximation.

    Parameters
    ----------
    dist : numpy.ndarray
        Distribution to be approximated
    N : int
        Number of samples to draw.
    rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
        A seed to initialize the BitGenerator. Defaults to None.
        If None, then fresh, unpredictable entropy will be pulled from the OS.

    """
    N_opt_samples,N_inv_samples = (max(_np.count_nonzero(dist),N),0)                                 # random subsampling if too little samples requested

    scale_,scale,inc_factor = (0.0,float(N_opt_samples),1.0)
    while (not _np.isclose(scale, scale_)) and (N_inv_samples != N_opt_samples):
        repeats = _np.rint(scale*dist).astype(_np.int64)
        N_inv_samples = _np.sum(repeats)
        scale_,scale,inc_factor = (scale,scale+inc_factor*0.5*(scale - scale_), inc_factor*2.0) \
                                   if N_inv_samples < N_opt_samples else \
                                  (scale_,0.5*(scale_ + scale), 1.0)

    return _np.repeat(_np.arange(len(dist)),repeats)[_np.random.default_rng(rng_seed).permutation(N_inv_samples)[:N]]


def shapeshifter(fro: _Tuple[int, ...],
                 to: _Tuple[int, ...],
                 mode: _Literal['left','right'] = 'left',                                           # noqa
                 keep_ones: bool = False) -> _Tuple[int, ...]:
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

    Example
    -------
    >>> import numpy as np
    >>> from damask import util
    >>> a = np.ones((3,4,2))
    >>> b = np.ones(4)
    >>> b_extended = b.reshape(util.shapeshifter(b.shape,a.shape))
    >>> (a * np.broadcast_to(b_extended,a.shape)).shape
    (3,4,2)


    """
    if len(fro) == 0 and len(to) == 0: return ()

    beg = dict(left ='(^.*\\b)',
               right='(^.*?\\b)')
    sep = dict(left ='(.*\\b)',
               right='(.*?\\b)')
    end = dict(left ='(.*?$)',
               right='(.*$)')
    fro = (1,) if len(fro) == 0 else fro
    to  = (1,) if len(to) == 0 else to
    try:
        match = _re.match(beg[mode]
                      +f',{sep[mode]}'.join(map(lambda x: f'{x}'
                                                          if x>1 or (keep_ones and len(fro)>1) else
                                                          '\\d+',fro))
                      +f',{end[mode]}',','.join(map(str,to))+',')
        assert match
        grp = match.groups()
    except AssertionError:
        raise ValueError(f'shapes cannot be shifted {fro} --> {to}')
    fill: _Any = ()
    for g,d in zip(grp,fro+(None,)):
        fill += (1,)*g.count(',')+(d,)
    return fill[:-1]


def shapeblender(a: _Tuple[int, ...],
                 b: _Tuple[int, ...]) -> _Tuple[int, ...]:
    """
    Return a shape that overlaps the rightmost entries of 'a' with the leftmost of 'b'.

    Parameters
    ----------
    a : tuple
        Shape of first array.
    b : tuple
        Shape of second array.

    Examples
    --------
    >>> shapeblender((4,4,3),(3,2,1))
        (4,4,3,2,1)
    >>> shapeblender((1,2),(1,2,3))
        (1,2,3)
    >>> shapeblender((1,),(2,2,1))
        (1,2,2,1)
    >>> shapeblender((3,2),(3,2))
        (3,2)

    """
    i = min(len(a),len(b))
    while i > 0 and a[-i:] != b[:i]: i -= 1
    return a + b[i:]


def extend_docstring(extra_docstring: str) -> _Callable:
    """
    Decorator: Append to function's docstring.

    Parameters
    ----------
    extra_docstring : str
       Docstring to append.

    """
    def _decorator(func):
        func.__doc__ += extra_docstring
        return func
    return _decorator


def extended_docstring(f: _Callable,
                       extra_docstring: str) -> _Callable:
    """
    Decorator: Combine another function's docstring with a given docstring.

    Parameters
    ----------
    f : function
       Function of which the docstring is taken.
    extra_docstring : str
       Docstring to append.

    """
    def _decorator(func):
        func.__doc__ = f.__doc__ + extra_docstring
        return func
    return _decorator


def DREAM3D_base_group(fname: _Union[str, _Path]) -> str:
    """
    Determine the base group of a DREAM.3D file.

    The base group is defined as the group (folder) that contains
    a 'SPACING' dataset in a '_SIMPL_GEOMETRY' group.

    Parameters
    ----------
    fname : str or pathlib.Path
        Filename of the DREAM.3D (HDF5) file.

    Returns
    -------
    path : str
        Path to the base group.

    """
    with _h5py.File(_Path(fname).expanduser(),'r') as f:
        base_group = f.visit(lambda path: path.rsplit('/',2)[0] if '_SIMPL_GEOMETRY/SPACING' in path else None)

    if base_group is None:
        raise ValueError(f'could not determine base group in file "{fname}"')

    return base_group

def DREAM3D_cell_data_group(fname: _Union[str, _Path]) -> str:
    """
    Determine the cell data group of a DREAM.3D file.

    The cell data group is defined as the group (folder) that contains
    a dataset in the base group whose length matches the total number
    of points as specified in '_SIMPL_GEOMETRY/DIMENSIONS'.

    Parameters
    ----------
    fname : str or pathlib.Path
        Filename of the DREAM.3D (HDF5) file.

    Returns
    -------
    path : str
        Path to the cell data group.

    """
    base_group = DREAM3D_base_group(fname)
    with _h5py.File(_Path(fname).expanduser(),'r') as f:
        cells = tuple(f['/'.join([base_group,'_SIMPL_GEOMETRY','DIMENSIONS'])][()][::-1])
        cell_data_group = f[base_group].visititems(lambda path,obj: path.split('/')[0] \
                                                   if isinstance(obj,_h5py._hl.dataset.Dataset) and _np.shape(obj)[:-1] == cells \
                                                   else None)

    if cell_data_group is None:
        raise ValueError(f'could not determine cell-data group in file "{fname}/{base_group}"')

    return cell_data_group


def Bravais_to_Miller(*,
                      uvtw: _np.ndarray = None,
                      hkil: _np.ndarray = None) -> _np.ndarray:
    """
    Transform 4 Miller–Bravais indices to 3 Miller indices of crystal direction [uvw] or plane normal (hkl).

    Parameters
    ----------
    uvtw|hkil : numpy.ndarray, shape (...,4)
        Miller–Bravais indices of crystallographic direction [uvtw] or plane normal (hkil).

    Returns
    -------
    uvw|hkl : numpy.ndarray, shape (...,3)
        Miller indices of [uvw] direction or (hkl) plane normal.

    """
    if (uvtw is not None) ^ (hkil is None):
        raise KeyError('specify either "uvtw" or "hkil"')
    axis,basis  = (_np.array(uvtw),_np.array([[1,0,-1,0],
                                              [0,1,-1,0],
                                              [0,0, 0,1]])) \
                  if hkil is None else \
                  (_np.array(hkil),_np.array([[1,0,0,0],
                                              [0,1,0,0],
                                              [0,0,0,1]]))
    return _np.einsum('il,...l',basis,axis)


def Miller_to_Bravais(*,
                      uvw: _np.ndarray = None,
                      hkl: _np.ndarray = None) -> _np.ndarray:
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
    axis,basis  = (_np.array(uvw),_np.array([[ 2,-1, 0],
                                             [-1, 2, 0],
                                             [-1,-1, 0],
                                             [ 0, 0, 3]])/3) \
                  if hkl is None else \
                  (_np.array(hkl),_np.array([[ 1, 0, 0],
                                             [ 0, 1, 0],
                                             [-1,-1, 0],
                                             [ 0, 0, 1]]))
    return _np.einsum('il,...l',basis,axis)


def dict_prune(d: _Dict) -> _Dict:
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


def dict_flatten(d: _Dict) -> _Dict:
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


def tail_repack(extended: _Union[str, _Sequence[str]],
                existing: _List[str] = []) -> _List[str]:
    """
    Repack tailing characters into single string if all are new.

    Parameters
    ----------
    extended : str or list of str
        Extended string list with potentially autosplitted tailing string relative to `existing`.
    existing : list of str
        Base string list.

    Returns
    -------
    repacked : list of str
        Repacked version of `extended`.

    Examples
    --------
    >>> tail_repack(['a','new','e','n','t','r','y'],['a','new'])
        ['a','new','entry']
    >>> tail_repack(['a','new','shiny','e','n','t','r','y'],['a','new'])
        ['a','new','shiny','e','n','t','r','y']

    """
    return [extended] if isinstance(extended,str) else existing + \
         ([''.join(extended[len(existing):])] if _np.prod([len(i) for i in extended[len(existing):]]) == 1 else
          list(extended[len(existing):]))


def aslist(arg: _Union[_IntCollection, int, None]) -> _List:
    """
    Transform argument to list.

    Parameters
    ----------
    arg : int or collection of int or None
        Entity to transform into list.

    Returns
    -------
    transformed : list
        Entity transformed into list.

    """
    return [] if arg is None else list(arg) if isinstance(arg,(_np.ndarray,_Collection)) else [arg]


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
        Set current time as basis for ETA estimation.

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

        _sys.stderr.write(f"{self.prefix} {'░'*self.bar_length}   0% ETA n/a")
        _sys.stderr.flush()

    def update(self,
               iteration: int) -> None:

        fraction = (iteration+1) / self.total

        if (filled_length := int(self.bar_length * fraction)) > int(self.bar_length * self.fraction_last) or \
            _datetime.datetime.now() - self.time_last_update > _datetime.timedelta(seconds=10):
            self.time_last_update = _datetime.datetime.now()
            bar = '█' * filled_length + '░' * (self.bar_length - filled_length)
            remaining_time = (_datetime.datetime.now() - self.time_start) \
                           * (self.total - (iteration+1)) / (iteration+1)
            remaining_time -= _datetime.timedelta(microseconds=remaining_time.microseconds)          # remove μs
            _sys.stderr.write(f'\r{self.prefix} {bar} {fraction:>4.0%} ETA {remaining_time}')
            _sys.stderr.flush()

        self.fraction_last = fraction

        if iteration == self.total - 1:
            _sys.stderr.write('\n')
            _sys.stderr.flush()
