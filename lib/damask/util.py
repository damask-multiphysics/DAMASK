# -*- coding: UTF-8 no BOM -*-

# damask utility functions
import sys,time,random,threading
import numpy as np
from optparse import OptionParser, Option

# -----------------------------
def emph(what):
# -----------------------------
  return '\033[1m'+str(what)+'\033[0m'

# -----------------------------
# Matlab like trigonometric functions that take and return angles in degrees.
# -----------------------------
for f in ['cos', 'sin', 'tan']:
  exec('def %sd(deg): return (np.%s(np.deg2rad(deg)))'%(f,f))
  exec('def a%sd(val): return (np.rad2deg(np.arc%s(val)))'%(f,f))


# -----------------------------
def gridLocation(idx,res):
# -----------------------------
  return ( idx  % res[0], \
         ( idx // res[0]) % res[1], \
         ( idx // res[0] // res[1]) % res[2] )


# -----------------------------
def gridIndex(location,res):
# -----------------------------
  return ( location[0] % res[0]                   + \
         ( location[1] % res[1]) * res[0]          + \
         ( location[2] % res[2]) * res[1] * res[0]   )


# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
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


# -----------------------------
class backgroundMessage(threading.Thread):
# -----------------------------
  choices = {'bounce':['_','o','O','°','¯','¯','°','O','o',],
             'circle': [u'\u25f4',u'\u25f5',u'\u25f6',u'\u25f7'],
             'hexagon': [u'\u2b22',u'\u2b23'],
             'pentagon': [u'\u2b20',u'\u2b54'],
             'square': [u'\u2596',u'\u2598',u'\u259d',u'\u2597'],
             'triangle': [u'\u140a',u'\u140a',u'\u1403',u'\u1405',u'\u1405',u'\u1403'],
             'amoeba': [u'\u2596',u'\u258f',u'\u2598',u'\u2594',u'\u259d',u'\u2595',u'\u2597',u'\u2582'],
             'beat': [u'\u2581',u'\u2582',u'\u2583',u'\u2584',u'\u2585',u'\u2586',u'\u2587',u'\u2588',u'\u2587',u'\u2586',u'\u2585',u'\u2584',u'\u2583',u'\u2582',],
             'prison': [u'\u168b',u'\u168c',u'\u168d',u'\u168f',u'\u168e',u'\u168d',u'\u168c',u'\u168b',],
             'breath': [u'\u1690',u'\u1691',u'\u1692',u'\u1693',u'\u1694',u'\u1693',u'\u1692',u'\u1691',u'\u1690'],
             'pulse': ['·','•',u'\u25cf',u'\u25cf','•',],
             'ant': [u'\u2801',u'\u2802',u'\u2810',u'\u2820',u'\u2804',u'\u2840',u'\u2880',u'\u2820',u'\u2804',u'\u2802',u'\u2810',u'\u2808'],
             'classic':['-', '\\', '|', '/',],
            }    

  def __init__(self,
               symbol,
               wait = 0.1):
    threading.Thread.__init__(self)
    self.message = ''
    self.new_message = ''
    self.counter = 0
    self.gap     = ' '
    self.symbols = self.choices[symbol if symbol in self.choices else random.choice(self.choices.keys())]
    self.waittime = wait
  
  def __quit__(self):
    length = len(self.symbols[self.counter] + self.gap + self.message)
    sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)
    sys.stderr.write('')
  
  def run(self):
    while not threading.enumerate()[0]._Thread__stopped:
      time.sleep(self.waittime)
      self.update_message()
    self.__quit__()

  def set_message(self, new_message):
    self.new_message = new_message
    self.print_message()
  
  def print_message(self):
    length = len(self.symbols[self.counter] + self.gap + self.message)
    sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)                                # delete former message
    sys.stderr.write(self.symbols[self.counter] + self.gap + self.new_message)                  # print new message
    self.message = self.new_message
      
  def update_message(self):
    self.counter = (self.counter + 1)%len(self.symbols)
    self.print_message()

'''
Non-linear least square fitting (Levenberg-Marquardt method) with 
the bounded parameters.
the codes of transformation between int <-> ext refers to the work of
Jonathan J. Helmus: https://github.com/jjhelmus/leastsqbound-scipy
other codes refers to the source code of minpack.py:
..\Lib\site-packages\scipy\optimize\minpack.py
'''
from numpy import (array, arcsin, asarray, cos, dot, eye, empty_like,
                   isscalar,finfo, take, triu, transpose, sqrt, sin)
from scipy.optimize import _minpack

def _check_func(checker, argname, thefunc, x0, args, numinputs,
                output_shape=None):
    from numpy import atleast_1d, shape, issubdtype, dtype, inexact
    ''' 
    The same as that of minpack.py, 
    '''
    res = atleast_1d(thefunc(*((x0[:numinputs],) + args)))
    if (output_shape is not None) and (shape(res) != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return shape(res)
            msg = "%s: there is a mismatch between the input and output " \
                  "shape of the '%s' argument" % (checker, argname)
            func_name = getattr(thefunc, '__name__', None)
            if func_name:
                msg += " '%s'." % func_name
            else:
                msg += "."
            raise TypeError(msg)
    if issubdtype(res.dtype, inexact):
        dt = res.dtype
    else:
        dt = dtype(float)
    return shape(res), dt

def _int2extGrad(p_int, bounds):
    """
    Calculate the gradients of transforming the internal (unconstrained) 
    to external (constained) parameter.
    """
    grad = empty_like(p_int)
    for i, (x, bound) in enumerate(zip(p_int, bounds)):
        lower, upper = bound
        if lower is None and upper is None:  # No constraints
            grad[i] = 1.0
        elif upper is None:                  # only lower bound
            grad[i] =  x/sqrt(x*x + 1.0)
        elif lower is None:                  # only upper bound
            grad[i] = -x/sqrt(x*x + 1.0)
        else:                                # lower and upper bounds
            grad[i] = (upper - lower)*cos(x)/2.0
    return grad

def _int2extFunc(bounds):
    '''
    transform internal parameters into external parameters.
    '''
    local = [_int2extLocal(b) for b in bounds]
    def _transform_i2e(p_int):
        p_ext = empty_like(p_int)
        p_ext[:] = [i(j) for i, j in zip(local, p_int)]
        return p_ext
    return _transform_i2e

def _ext2intFunc(bounds):
    '''
    transform external parameters into internal parameters.
    '''
    local = [_ext2intLocal(b) for b in bounds]
    def _transform_e2i(p_ext):
        p_int = empty_like(p_ext)
        p_int[:] = [i(j) for i, j in zip(local, p_ext)]
        return p_int
    return _transform_e2i
    
def _int2extLocal(bound):
    '''
    transform a single internal parameter to an external parameter.
    '''
    lower, upper = bound
    if lower is None and upper is None:      # no constraints
        return lambda x: x
    elif upper is None:                      # only lower bound
        return lambda x: lower - 1.0 + sqrt(x*x + 1.0)
    elif lower is None:                      # only upper bound
        return lambda x: upper + 1.0 - sqrt(x*x + 1.0)
    else:
        return lambda x: lower + ((upper - lower)/2.0)*(sin(x) + 1.0)

def _ext2intLocal(bound):
    '''
    transform a single external parameter to an internal parameter.
    '''
    lower, upper = bound
    if lower is None and upper is None:  # no constraints
        return lambda x: x
    elif upper is None:                  # only lower bound
        return lambda x: sqrt((x - lower + 1.0)**2 - 1.0)
    elif lower is None:                  # only upper bound
        return lambda x: sqrt((x - upper - 1.0)**2 - 1.0)
    else:
        return lambda x: arcsin((2.0*(x - lower)/(upper - lower)) - 1.0)

def leastsqBound(func, x0, args=(), bounds=None, Dfun=None, full_output=0,
                 col_deriv=0, ftol=1.49012e-8, xtol=1.49012e-8,
                 gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None):
    '''
    An internal parameter list is used to enforce contraints on the fitting 
    parameters. The transfomation is based on that of MINUIT package. 
    please see: F. James and M. Winkler. MINUIT User's Guide, 2004.

    bounds : list
      (min, max) pairs for each parameter, use None for 'min' or 'max' 
      when there is no bound in that direction. 
      For example: if there are two parameters needed to be fitting, then
      bounds is [(min1,max1), (min2,max2)]
      
    This function is based on 'leastsq' of minpack.py, the annotation of
    other parameters can be found in 'leastsq'.
    ..\Lib\site-packages\scipy\optimize\minpack.py
    '''
    i2e = _int2extFunc(bounds)
    e2i = _ext2intFunc(bounds)
    
    x0 = asarray(x0).flatten()
    n = len(x0)

    if len(bounds) != n:
        raise ValueError('the length of bounds is inconsistent with the number of parameters ')
    
    if not isinstance(args, tuple):
        args = (args,)
    
    shape, dtype = _check_func('leastsq', 'func', func, x0, args, n)
    m = shape[0]

    if n > m:
        raise TypeError('Improper input: N=%s must not exceed M=%s' % (n, m))
    if epsfcn is None:
        epsfcn = finfo(dtype).eps

    # wrapped func
    def funcWarp(x, *args):
        return func(i2e(x), *args)

    xi0 = e2i(x0)
    
    if Dfun is None:
        if maxfev == 0:
            maxfev = 200*(n + 1)
        retval = _minpack._lmdif(funcWarp, xi0, args, full_output, ftol, xtol,
                                 gtol, maxfev, epsfcn, factor, diag)
    else:
        if col_deriv:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (n, m))
        else:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (m, n))
        if maxfev == 0:
            maxfev = 100*(n + 1)

        # wrapped Dfun
        def DfunWarp(x, *args):
            return Dfun(i2e(x), *args)

        retval = _minpack._lmder(funcWarp, DfunWarp, xi0, args, full_output, col_deriv,
                                 ftol, xtol, gtol, maxfev, factor, diag)

    errors = {0: ["Improper input parameters.", TypeError],
              1: ["Both actual and predicted relative reductions "
                  "in the sum of squares\n  are at most %f" % ftol, None],
              2: ["The relative error between two consecutive "
                  "iterates is at most %f" % xtol, None],
              3: ["Both actual and predicted relative reductions in "
                  "the sum of squares\n  are at most %f and the "
                  "relative error between two consecutive "
                  "iterates is at \n  most %f" % (ftol, xtol), None],
              4: ["The cosine of the angle between func(x) and any "
                  "column of the\n  Jacobian is at most %f in "
                  "absolute value" % gtol, None],
              5: ["Number of calls to function has reached "
                  "maxfev = %d." % maxfev, ValueError],
              6: ["ftol=%f is too small, no further reduction "
                  "in the sum of squares\n  is possible.""" % ftol,
                  ValueError],
              7: ["xtol=%f is too small, no further improvement in "
                  "the approximate\n  solution is possible." % xtol,
                  ValueError],
              8: ["gtol=%f is too small, func(x) is orthogonal to the "
                  "columns of\n  the Jacobian to machine "
                  "precision." % gtol, ValueError],
              'unknown': ["Unknown error.", TypeError]}

    info = retval[-1]    # The FORTRAN return value
    
    if info not in [1, 2, 3, 4] and not full_output:
        if info in [5, 6, 7, 8]:
            warnings.warn(errors[info][0], RuntimeWarning)
        else:
            try:
                raise errors[info][1](errors[info][0])
            except KeyError:
                raise errors['unknown'][1](errors['unknown'][0])

    mesg = errors[info][0]
    x = i2e(retval[0])

    if full_output:
        grad = _int2extGrad(retval[0], bounds)
        retval[1]['fjac'] = (retval[1]['fjac'].T / take(grad,
                             retval[1]['ipvt'] - 1)).T
        cov_x = None
        if info in [1, 2, 3, 4]:
            from numpy.dual import inv
            from numpy.linalg import LinAlgError
            perm = take(eye(n), retval[1]['ipvt'] - 1, 0)
            r = triu(transpose(retval[1]['fjac'])[:n, :])
            R = dot(r, perm)
            try:
                cov_x = inv(dot(transpose(R), R))
            except LinAlgError as inverror:
                print inverror
                pass
        return (x, cov_x) + retval[1:-1] + (mesg, info)
    else:
        return (x, info)

def _general_function(params, ydata, xdata, function):
    return  function(xdata, *params) - ydata
def _weighted_general_function(params, ydata, xdata, function, weights):
    return (function(xdata, *params) - ydata)*weights

def curve_fit_bound(f, xdata, ydata, p0=None, sigma=None, bounds=None, **kw):
    ''' Similar as 'curve_fit' in minpack.py'''
    if p0 is None:
        # determine number of parameters by inspecting the function
        import inspect
        args, varargs, varkw, defaults = inspect.getargspec(f)
        if len(args) < 2:
            msg = "Unable to determine number of fit parameters."
            raise ValueError(msg)
        if 'self' in args:
            p0 = [1.0] * (len(args)-2)
        else:
            p0 = [1.0] * (len(args)-1)

    if isscalar(p0):
        p0 = array([p0])

    args = (ydata, xdata, f)
    if sigma is None:
        func = _general_function
    else:
        func = _weighted_general_function
        args += (1.0/asarray(sigma),)

    return_full = kw.pop('full_output', False)
    res = leastsqBound(func, p0, args=args, bounds = bounds, full_output=True, **kw)
    (popt, pcov, infodict, errmsg, ier) = res

    if ier not in [1, 2, 3, 4]:
        msg = "Optimal parameters not found: " + errmsg
        raise RuntimeError(msg)

    if (len(ydata) > len(p0)) and pcov is not None:
        s_sq = (func(popt, *args)**2).sum()/(len(ydata)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = inf

    if return_full:
        return popt, pcov, infodict, errmsg, ier
    else:
        return popt, pcov