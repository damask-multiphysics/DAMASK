import sys
import datetime
import os
import subprocess
import shlex
import re
import fractions
from functools import reduce

import numpy as np
import h5py

from . import version

# limit visibility
__all__=[
         'srepr',
         'emph','deemph','warn','strikeout',
         'execute',
         'natural_sort',
         'show_progress',
         'scale_to_coprime',
         'project_stereographic',
         'hybrid_IA',
         'execution_stamp',
         'shapeshifter', 'shapeblender',
         'extend_docstring', 'extended_docstring',
         'DREAM3D_base_group', 'DREAM3D_cell_data_group',
         'dict_prune', 'dict_flatten'
        ]

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
def srepr(arg,glue = '\n'):
    r"""
    Join arguments with glue string.

    Parameters
    ----------
    arg : iterable
        Items to join.
    glue : str, optional
        Glue used for joining operation. Defaults to \n.

    """
    if (not hasattr(arg, 'strip') and
           (hasattr(arg, '__getitem__') or
            hasattr(arg, '__iter__'))):
        return glue.join(str(x) for x in arg)
    else:
       return arg if isinstance(arg,str) else repr(arg)


def emph(what):
    """Formats string with emphasis."""
    return _colors['bold']+srepr(what)+_colors['end_color']

def deemph(what):
    """Formats string with deemphasis."""
    return _colors['dim']+srepr(what)+_colors['end_color']

def warn(what):
    """Formats string for warning."""
    return _colors['warning']+emph(what)+_colors['end_color']

def strikeout(what):
    """Formats string as strikeout."""
    return _colors['crossout']+srepr(what)+_colors['end_color']


def execute(cmd,wd='./',env=None):
    """
    Execute command.

    Parameters
    ----------
    cmd : str
        Command to be executed.
    wd : str, optional
        Working directory of process. Defaults to ./ .
    env : dict, optional
        Environment for execution.

    """
    print(f"executing '{cmd}' in '{wd}'")
    process = subprocess.run(shlex.split(cmd),
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE,
                             env = os.environ if env is None else env,
                             cwd = wd,
                             encoding = 'utf-8')

    if process.returncode != 0:
        print(process.stdout)
        print(process.stderr)
        raise RuntimeError(f"'{cmd}' failed with returncode {process.returncode}")

    return process.stdout, process.stderr


def natural_sort(key):
    convert = lambda text: int(text) if text.isdigit() else text
    return [ convert(c) for c in re.split('([0-9]+)', key) ]


def show_progress(iterable,N_iter=None,prefix='',bar_length=50):
    """
    Decorate a loop with a status bar.

    Use similar like enumerate.

    Parameters
    ----------
    iterable : iterable/function with yield statement
        Iterable (or function with yield statement) to be decorated.
    N_iter : int
        Total # of iterations. Needed if number of iterations can not be obtained as len(iterable).
    prefix : str, optional.
        Prefix string.
    bar_length : int, optional
        Character length of bar. Defaults to 50.

    """
    if N_iter in [0,1] or (hasattr(iterable,'__len__') and len(iterable) <= 1):
        for item in iterable:
            yield item
    else:
        status = _ProgressBar(N_iter if N_iter is not None else len(iterable),prefix,bar_length)

        for i,item in enumerate(iterable):
            yield item
            status.update(i)


def scale_to_coprime(v):
    """
    Scale vector to co-prime (relatively prime) integers.

    Parameters
    ----------
    v : numpy.ndarray of shape (:)
        Vector to scale.

    """
    MAX_DENOMINATOR = 1000000

    def get_square_denominator(x):
        """Denominator of the square of a number."""
        return fractions.Fraction(x ** 2).limit_denominator(MAX_DENOMINATOR).denominator

    def lcm(a,b):
        """Least common multiple."""
        try:
            return np.lcm(a,b)                                                                      # numpy > 1.18
        except AttributeError:
            return a * b // np.gcd(a, b)

    m = (np.array(v) * reduce(lcm, map(lambda x: int(get_square_denominator(x)),v)) ** 0.5).astype(int)
    m = m//reduce(np.gcd,m)

    with np.errstate(invalid='ignore'):
        if not np.allclose(np.ma.masked_invalid(v/m),v[np.argmax(abs(v))]/m[np.argmax(abs(v))]):
            raise ValueError(f'Invalid result {m} for input {v}. Insufficient precision?')

    return m


def project_stereographic(vector,direction='z',normalize=True,keepdims=False):
    """
    Apply stereographic projection to vector.

    Parameters
    ----------
    vector : numpy.ndarray of shape (...,3)
        Vector coordinates to be projected.
    direction : str
        Projection direction 'x', 'y', or 'z'.
        Defaults to 'z'.
    normalize : bool
        Ensure unit length of input vector. Defaults to True.
    keepdims : bool
        Maintain three-dimensional output coordinates.
        Default two-dimensional output uses right-handed frame spanned by
        the next and next-next axis relative to the projection direction,
        e.g. x-y when projecting along z and z-x when projecting along y.

    Returns
    -------
    coordinates : numpy.ndarray of shape (...,2 | 3)
        Projected coordinates.

    Examples
    --------
    >>> project_stereographic(np.ones(3))
        [0.3660254, 0.3660254]
    >>> project_stereographic(np.ones(3),direction='x',normalize=False,keepdims=True)
        [0, 0.5, 0.5]
    >>> project_stereographic([0,1,1],direction='y',normalize=True,keepdims=False)
        [0.41421356, 0]

    """
    shift = 'zyx'.index(direction)
    v_ = np.roll(vector/np.linalg.norm(vector,axis=-1,keepdims=True) if normalize else vector,
                 shift,axis=-1)
    return np.roll(np.block([v_[...,:2]/(1+np.abs(v_[...,2:3])),np.zeros_like(v_[...,2:3])]),
                   -shift if keepdims else 0,axis=-1)[...,:3 if keepdims else 2]


def execution_stamp(class_name,function_name=None):
    """Timestamp the execution of a (function within a) class."""
    now = datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S%z')
    _function_name = '' if function_name is None else f'.{function_name}'
    return f'damask.{class_name}{_function_name} v{version} ({now})'


def hybrid_IA(dist,N,rng_seed=None):
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
    N_opt_samples,N_inv_samples = (max(np.count_nonzero(dist),N),0)                                 # random subsampling if too little samples requested

    scale_,scale,inc_factor = (0.0,float(N_opt_samples),1.0)
    while (not np.isclose(scale, scale_)) and (N_inv_samples != N_opt_samples):
        repeats = np.rint(scale*dist).astype(int)
        N_inv_samples = np.sum(repeats)
        scale_,scale,inc_factor = (scale,scale+inc_factor*0.5*(scale - scale_), inc_factor*2.0) \
                                   if N_inv_samples < N_opt_samples else \
                                  (scale_,0.5*(scale_ + scale), 1.0)

    return np.repeat(np.arange(len(dist)),repeats)[np.random.default_rng(rng_seed).permutation(N_inv_samples)[:N]]


def shapeshifter(fro,to,mode='left',keep_ones=False):
    """
    Return a tuple that reshapes 'fro' to become broadcastable to 'to'.

    Parameters
    ----------
    fro : tuple
        Original shape of array.
    to : tuple
        Target shape of array after broadcasting.
        len(to) cannot be less than len(fro).
    mode : str, optional
        Indicates whether new axes are preferably added to
        either 'left' or 'right' of the original shape.
        Defaults to 'left'.
    keep_ones : bool, optional
        Treat '1' in fro as literal value instead of dimensional placeholder.
        Defaults to False.

    """
    beg = dict(left ='(^.*\\b)',
               right='(^.*?\\b)')
    sep = dict(left ='(.*\\b)',
               right='(.*?\\b)')
    end = dict(left ='(.*?$)',
               right='(.*$)')
    fro = (1,) if not len(fro) else fro
    to  = (1,) if not len(to)  else to
    try:
        grp = re.match(beg[mode]
                      +f',{sep[mode]}'.join(map(lambda x: f'{x}'
                                                          if x>1 or (keep_ones and len(fro)>1) else
                                                          '\\d+',fro))
                      +f',{end[mode]}',
                       ','.join(map(str,to))+',').groups()
    except AttributeError:
        raise ValueError(f'Shapes can not be shifted {fro} --> {to}')
    fill = ()
    for g,d in zip(grp,fro+(None,)):
        fill += (1,)*g.count(',')+(d,)
    return fill[:-1]


def shapeblender(a,b):
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


def extend_docstring(extra_docstring):
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


def extended_docstring(f,extra_docstring):
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


def DREAM3D_base_group(fname):
    """
    Determine the base group of a DREAM.3D file.

    The base group is defined as the group (folder) that contains
    a 'SPACING' dataset in a '_SIMPL_GEOMETRY' group.

    Parameters
    ----------
    fname : str
        Filename of the DREAM.3D (HDF5) file.

    """
    with h5py.File(fname,'r') as f:
        base_group = f.visit(lambda path: path.rsplit('/',2)[0] if '_SIMPL_GEOMETRY/SPACING' in path else None)

    if base_group is None:
        raise ValueError(f'Could not determine base group in file {fname}.')

    return base_group

def DREAM3D_cell_data_group(fname):
    """
    Determine the cell data group of a DREAM.3D file.

    The cell data group is defined as the group (folder) that contains
    a dataset in the base group whose length matches the total number
    of points as specified in '_SIMPL_GEOMETRY/DIMENSIONS'.

    Parameters
    ----------
    fname : str
        Filename of the DREAM.3D (HDF5) file.

    """
    base_group = DREAM3D_base_group(fname)
    with h5py.File(fname,'r') as f:
        cells = tuple(f['/'.join([base_group,'_SIMPL_GEOMETRY','DIMENSIONS'])][()][::-1])
        cell_data_group = f[base_group].visititems(lambda path,obj: path.split('/')[0] \
                                                   if isinstance(obj,h5py._hl.dataset.Dataset) and np.shape(obj)[:-1] == cells \
                                                   else None)

    if cell_data_group is None:
        raise ValueError(f'Could not determine cell data group in file {fname}/{base_group}.')

    return cell_data_group


def dict_prune(d):
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


def dict_flatten(d):
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



####################################################################################################
# Classes
####################################################################################################
class _ProgressBar:
    """
    Report progress of an interation as a status bar.

    Works for 0-based loops, ETA is estimated by linear extrapolation.
    """

    def __init__(self,total,prefix,bar_length):
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
        self.start_time = datetime.datetime.now()
        self.last_fraction = 0.0

        sys.stderr.write(f"{self.prefix} {'░'*self.bar_length}   0% ETA n/a")
        sys.stderr.flush()

    def update(self,iteration):

        fraction = (iteration+1) / self.total
        filled_length = int(self.bar_length * fraction)

        if filled_length > int(self.bar_length * self.last_fraction):
            bar = '█' * filled_length + '░' * (self.bar_length - filled_length)
            delta_time = datetime.datetime.now() - self.start_time
            remaining_time = (self.total - (iteration+1)) * delta_time / (iteration+1)
            remaining_time -= datetime.timedelta(microseconds=remaining_time.microseconds)          # remove μs
            sys.stderr.write(f'\r{self.prefix} {bar} {fraction:>4.0%} ETA {remaining_time}')
            sys.stderr.flush()

        self.last_fraction = fraction

        if iteration == self.total - 1:
            sys.stderr.write('\n')
            sys.stderr.flush()
