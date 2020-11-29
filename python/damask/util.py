import sys
import datetime
import os
import subprocess
import shlex
import re
import fractions
from functools import reduce
from optparse import Option

import numpy as np

from . import version

# limit visibility
__all__=[
         'srepr',
         'croak',
         'report',
         'emph','deemph','warn','strikeout',
         'execute',
         'show_progress',
         'scale_to_coprime',
         'project_stereographic',
         'hybrid_IA',
         'return_message',
         'extendableOption',
         'execution_stamp',
         'shapeshifter', 'shapeblender',
         'extend_docstring', 'extended_docstring'
        ]

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
    if (not hasattr(arg, "strip") and
           (hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))):
        return glue.join(str(x) for x in arg)
    return arg if isinstance(arg,str) else repr(arg)


def croak(what, newline = True):
    """
    Write formated to stderr.

    DEPRECATED

    Parameters
    ----------
    what : str or iterable
        Content to be displayed.
    newline : bool, optional
        Separate items of what by newline. Defaults to True.

    """
    if what is not None:
        sys.stderr.write(srepr(what,glue = '\n') + ('\n' if newline else ''))
    sys.stderr.flush()


def report(who = None,
           what = None):
    """
    Report script and file name.

    DEPRECATED

    """
    croak( (emph(who)+': ' if who is not None else '') + (what if what is not None else '') + '\n' )


def emph(what):
    """Formats string with emphasis."""
    return bcolors.BOLD+srepr(what)+bcolors.ENDC

def deemph(what):
    """Formats string with deemphasis."""
    return bcolors.DIM+srepr(what)+bcolors.ENDC

def warn(what):
    """Formats string for warning."""
    return bcolors.WARNING+emph(what)+bcolors.ENDC

def strikeout(what):
    """Formats string as strikeout."""
    return bcolors.CROSSOUT+srepr(what)+bcolors.ENDC


def execute(cmd,
            stream_in = None,
            wd = './',
            env = None):
    """
    Execute command.

    Parameters
    ----------
    cmd : str
        Command to be executed.
    stream_in : file object, optional
        Input (via pipe) for executed process.
    wd : str, optional
        Working directory of process. Defaults to ./ .
    env : dict, optional
        Environment for execution.

    """
    initialPath = os.getcwd()
    myEnv = os.environ if env is None else env
    os.chdir(wd)
    print(f"executing '{cmd}' in '{wd}'")
    process = subprocess.Popen(shlex.split(cmd),
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               stdin  = subprocess.PIPE,
                               env = myEnv)
    stdout, stderr = [i for i in (process.communicate() if stream_in is None
                             else process.communicate(stream_in.read().encode('utf-8')))]
    os.chdir(initialPath)
    stdout = stdout.decode('utf-8').replace('\x08','')
    stderr = stderr.decode('utf-8').replace('\x08','')
    if process.returncode != 0:
        raise RuntimeError(f"'{cmd}' failed with returncode {process.returncode}")
    return stdout, stderr


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
    status = _ProgressBar(N_iter if N_iter else len(iterable),prefix,bar_length)

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

    def lcm(a, b):
        """Least common multiple."""
        # Python 3.9 provides math.lcm, see https://stackoverflow.com/questions/51716916.
        return a * b // np.gcd(a, b)

    m = (np.array(v) * reduce(lcm, map(lambda x: int(get_square_denominator(x)),v)) ** 0.5).astype(np.int)
    m = m//reduce(np.gcd,m)

    with np.errstate(invalid='ignore'):
        if not np.allclose(np.ma.masked_invalid(v/m),v[np.argmax(abs(v))]/m[np.argmax(abs(v))]):
            raise ValueError(f'Invalid result {m} for input {v}. Insufficient precision?')

    return m


def project_stereographic(vector,normalize=False):
    """
    Apply stereographic projection to vector.

    Parameters
    ----------
    vector : numpy.ndarray of shape (...,3)
        Vector coordinates to be projected.
    normalize : bool
        Ensure unit length for vector. Defaults to False.

    Returns
    -------
    coordinates : numpy.ndarray of shape (...,2)
        Projected coordinates.

    """
    v_ = vector/np.linalg.norm(vector,axis=-1,keepdims=True) if normalize else vector
    return np.block([v_[...,:2]/(1+np.abs(v_[...,2:3])),
                     np.zeros_like(v_[...,2:3])])


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


####################################################################################################
# Classes
####################################################################################################
class extendableOption(Option):
    """
    Used for definition of new option parser action 'extend', which enables to take multiple option arguments.

    Adopted from online tutorial http://docs.python.org/library/optparse.html
    DEPRECATED
    """

    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            lvalue = value.split(",")
            values.ensure_value(dest, []).extend(lvalue)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)


class _ProgressBar:
    """
    Report progress of an interation as a status bar.

    Works for 0-based loops, ETA is estimated by linear extrapolation.
    """

    def __init__(self,total,prefix,bar_length):
        """
        Inititalize a progress bar to current time as basis for ETA estimation.

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
            remaining_time -= datetime.timedelta(microseconds=remaining_time.microseconds)           # remove μs
            sys.stderr.write(f'\r{self.prefix} {bar} {fraction:>4.0%} ETA {remaining_time}')
            sys.stderr.flush()

        self.last_fraction = fraction

        if iteration == self.total - 1:
            sys.stderr.write('\n')
            sys.stderr.flush()


class bcolors:
    """
    ASCII colors.

    https://svn.blender.org/svnroot/bf-blender/trunk/blender/build_files/scons/tools/bcolors.py
    https://stackoverflow.com/questions/287871
    """

    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    DIM       = '\033[2m'
    UNDERLINE = '\033[4m'
    CROSSOUT  = '\033[9m'


class return_message:
    """Object with formatted return message."""

    def __init__(self,message):
        """
        Sets return message.

        Parameters
        ----------
        message : str or list of str
            message for output to screen

        """
        self.message = message

    def __repr__(self):
        """Return message suitable for interactive shells."""
        return srepr(self.message)
