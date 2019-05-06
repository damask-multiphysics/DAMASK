# -*- coding: UTF-8 no BOM -*-

import math
import numpy as np
from . import Lambert
from .quaternion import Quaternion
from .quaternion import P


####################################################################################################
class Rotation:
    u"""
    Orientation stored as unit quaternion.
    
    Following: D Rowenhorst et al. Consistent representations of and conversions between 3D rotations
               10.1088/0965-0393/23/8/083501
    Convention 1: coordinate frames are right-handed
    Convention 2: a rotation angle ω is taken to be positive for a counterclockwise rotation
                  when viewing from the end point of the rotation axis towards the origin
    Convention 3: rotations will be interpreted in the passive sense
    Convention 4: Euler angle triplets are implemented using the Bunge convention,
                  with the angular ranges as [0, 2π],[0, π],[0, 2π]
    Convention 5: the rotation angle ω is limited to the interval [0, π]
    Convention 6: P = -1 (as default)
    
    q is the real part, p = (x, y, z) are the imaginary parts.

    Vector "a" (defined in coordinate system "A") is passively rotated
               resulting in new coordinates "b" when expressed in system "B".
    b = Q * a
    b = np.dot(Q.asMatrix(),a)
    """

    __slots__ = ['quaternion']
    
    def __init__(self,quaternion = np.array([1.0,0.0,0.0,0.0])):
      """
      Initializes to identity unless specified
      
      If a quaternion is given, it needs to comply with the convection. Use .fromQuaternion
      to check the input.
      """
      if isinstance(quaternion,Quaternion):
        self.quaternion = quaternion.copy()
      else:
        self.quaternion = Quaternion(q=quaternion[0],p=quaternion[1:4])
        
    def __copy__(self):
      """Copy"""
      return self.__class__(self.quaternion)
      
    copy = __copy__


    def __repr__(self):
      """Value in selected representation"""
      return '\n'.join([
            '{}'.format(self.quaternion),
            'Matrix:\n{}'.format( '\n'.join(['\t'.join(list(map(str,self.asMatrix()[i,:]))) for i in range(3)]) ),
            'Bunge Eulers / deg: {}'.format('\t'.join(list(map(str,self.asEulers(degrees=True)))) ),
              ])
              
    def __mul__(self, other):
      """
      Multiplication
      
      Rotation: Details needed (active/passive), rotation of (3,3,3,3)-matrix should be considered
      """
      if isinstance(other, Rotation):                                                               # rotate a rotation
        return self.__class__(self.quaternion * other.quaternion).standardize()
      elif isinstance(other, np.ndarray):
        if other.shape == (3,):                                                                     # rotate a single (3)-vector
          ( x, y, z) = self.quaternion.p
          (Vx,Vy,Vz) = other[0:3]
          A = self.quaternion.q*self.quaternion.q - np.dot(self.quaternion.p,self.quaternion.p)
          B = 2.0 * (x*Vx + y*Vy + z*Vz)
          C = 2.0 * P*self.quaternion.q

          return np.array([
            A*Vx + B*x + C*(y*Vz - z*Vy),
            A*Vy + B*y + C*(z*Vx - x*Vz),
            A*Vz + B*z + C*(x*Vy - y*Vx),
            ])
        elif other.shape == (3,3,):                                                                 # rotate a single (3x3)-matrix
           return np.dot(self.asMatrix(),np.dot(other,self.asMatrix().T))
        elif other.shape == (3,3,3,3,):
          raise NotImplementedError
        else:
          return NotImplemented
      elif isinstance(other, tuple):                                                                # used to rotate a meshgrid-tuple
          ( x, y, z) = self.quaternion.p
          (Vx,Vy,Vz) = other[0:3]
          A = self.quaternion.q*self.quaternion.q - np.dot(self.quaternion.p,self.quaternion.p)
          B = 2.0 * (x*Vx + y*Vy + z*Vz)
          C = 2.0 * P*self.quaternion.q

          return np.array([
            A*Vx + B*x + C*(y*Vz - z*Vy),
            A*Vy + B*y + C*(z*Vx - x*Vz),
            A*Vz + B*z + C*(x*Vy - y*Vx),
            ])
      else:
        return NotImplemented
        
    
    def inverse(self):
      """In-place inverse rotation/backward rotation"""
      self.quaternion.conjugate()
      return self
      
    def inversed(self):
      """Inverse rotation/backward rotation"""
      return self.copy().inverse()


    def standardize(self):
      """In-place quaternion representation with positive q"""
      if self.quaternion.q < 0.0: self.quaternion.homomorph()
      return self

    def standardized(self):
      """Quaternion representation with positive q"""
      return self.copy().standardize()


    def misorientation(self,other):
      """Misorientation"""
      return other*self.inversed()


    def average(self,other):
      """Calculate the average rotation"""
      return Rotation.fromAverage([self,other])
      

    ################################################################################################
    # convert to different orientation representations (numpy arrays)

    def asQuaternion(self):
      """Unit quaternion: (q, [p_1, p_2, p_3])"""
      return self.quaternion.asArray()
      
    def asEulers(self,
                 degrees = False):
      """Bunge-Euler angles: (φ_1, ϕ, φ_2)"""
      eu = qu2eu(self.quaternion.asArray())
      if degrees: eu = np.degrees(eu)
      return eu
    
    def asAxisAngle(self,
                    degrees = False):
      """Axis-angle pair: ([n_1, n_2, n_3], ω)"""
      ax = qu2ax(self.quaternion.asArray())
      if degrees: ax[3] = np.degrees(ax[3])
      return ax
      
    def asMatrix(self):
      """Rotation matrix"""
      return qu2om(self.quaternion.asArray())

    def asRodrigues(self,vector=False):
      """Rodrigues-Frank vector: ([n_1, n_2, n_3], tan(ω/2))"""
      ro = qu2ro(self.quaternion.asArray())
      return ro[:3]*ro[3] if vector else ro
  
    def asHomochoric(self):
      """Homochoric vector: (h_1, h_2, h_3)"""
      return qu2ho(self.quaternion.asArray())
         
    def asCubochoric(self):
      return qu2cu(self.quaternion.asArray())
      
    def asM(self):
      """Intermediate representation supporting quaternion averaging (see F. Landis Markley et al.)"""
      return self.quaternion.asM()
      

    ################################################################################################
    # static constructors. The input data needs to follow the convention, options allow to
    # relax these convections
    @classmethod
    def fromQuaternion(cls,
                       quaternion,
                       acceptHomomorph = False,
                       P = -1):

      qu = quaternion if isinstance(quaternion, np.ndarray) else np.array(quaternion)
      if P > 0: qu[1:4] *= -1                                                                       # convert from P=1 to P=-1
      if qu[0] < 0.0:
        if acceptHomomorph:
          qu *= -1.
        else:
          raise ValueError('Quaternion has negative first component.\n{}'.format(qu[0]))
      if not np.isclose(np.linalg.norm(qu), 1.0):
        raise ValueError('Quaternion is not of unit length.\n{} {} {} {}'.format(*qu))

      return cls(qu)
      
    @classmethod
    def fromEulers(cls,
                   eulers,
                   degrees = False):

      eu = eulers if isinstance(eulers, np.ndarray) else np.array(eulers)
      eu = np.radians(eu) if degrees else eu
      if np.any(eu < 0.0) or np.any(eu > 2.0*np.pi) or eu[1] > np.pi:
        raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π].\n{} {} {}.'.format(*eu))
        
      return cls(eu2qu(eu))
      
    @classmethod
    def fromAxisAngle(cls,
                      angleAxis,
                      degrees = False,
                      normalise = False,
                      P = -1):
      
      ax = angleAxis if isinstance(angleAxis, np.ndarray) else np.array(angleAxis)
      if P > 0: ax[0:3] *= -1                                                                       # convert from P=1 to P=-1
      if degrees:   ax[3] = np.radians(ax[3])
      if normalise: ax[0:3] /=np.linalg.norm(ax[0:3])
      if ax[3] < 0.0 or ax[3] > np.pi:
        raise ValueError('Axis angle rotation angle outside of [0..π].\n'.format(ax[3]))
      if not np.isclose(np.linalg.norm(ax[0:3]), 1.0):
        raise ValueError('Axis angle rotation axis is not of unit length.\n{} {} {}'.format(*ax[0:3]))
      
      return cls(ax2qu(ax))
      
    @classmethod
    def fromMatrix(cls,
                   matrix,
                   containsStretch = False):                                 #ToDo: better name?
      
      om = matrix if isinstance(matrix, np.ndarray) else np.array(matrix).reshape((3,3))            # ToDo: Reshape here or require explicit?
      if containsStretch:
        (U,S,Vh) = np.linalg.svd(om)                                                                # singular value decomposition
        om = np.dot(U,Vh)     
      if not np.isclose(np.linalg.det(om),1.0):
        raise ValueError('matrix is not a proper rotation.\n{}'.format(om))
      if    not np.isclose(np.dot(om[0],om[1]), 0.0) \
         or not np.isclose(np.dot(om[1],om[2]), 0.0) \
         or not np.isclose(np.dot(om[2],om[0]), 0.0):
        raise ValueError('matrix is not orthogonal.\n{}'.format(om))

      return cls(om2qu(om))
      
    @classmethod
    def fromRodrigues(cls,
                      rodrigues,
                      normalise = False,
                      P = -1):
        
      ro = rodrigues if isinstance(rodrigues, np.ndarray) else np.array(rodrigues) 
      if P > 0: ro[0:3] *= -1                                                                       # convert from P=1 to P=-1
      if normalise: ro[0:3] /=np.linalg.norm(ro[0:3])
      if not np.isclose(np.linalg.norm(ro[0:3]), 1.0):
        raise ValueError('Rodrigues rotation axis is not of unit length.\n{} {} {}'.format(*ro[0:3]))
      if ro[3] < 0.0:
        raise ValueError('Rodriques rotation angle not positive.\n'.format(ro[3]))
        
      return cls(ro2qu(ro))
      
    @classmethod
    def fromHomochoric(cls,
                      homochoric,
                      P = -1):
        
      ho = homochoric if isinstance(homochoric, np.ndarray) else np.array(homochoric) 
      if P > 0: ho *= -1                                                                            # convert from P=1 to P=-1

      return cls(ho2qu(ho))
      
    @classmethod
    def fromCubochoric(cls,
                      cubochoric,
                      P = -1):
      
      cu = cubochoric if isinstance(cubochoric, np.ndarray) else np.array(cubochoric)   
      ho = cu2ho(cu)
      if P > 0: ho *= -1                                                                            # convert from P=1 to P=-1
        
      return cls(ho2qu(ho))


    @classmethod
    def fromAverage(cls,
                    rotations,
                    weights = []):
      """
      Average rotation

      ref: F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
           Averaging Quaternions,
           Journal of Guidance, Control, and Dynamics, Vol. 30, No. 4 (2007), pp. 1193-1197.
           doi: 10.2514/1.28949
           
      usage: input a list of rotations and optional weights
      """
      if not all(isinstance(item, Rotation) for item in rotations):
        raise TypeError("Only instances of Rotation can be averaged.")

      N = len(rotations)
      if weights == [] or not weights:
        weights = np.ones(N,dtype='i')

      for i,(r,n) in enumerate(zip(rotations,weights)):
        M =          r.asM() * n if i == 0 \
            else M + r.asM() * n                                                                    # noqa add (multiples) of this rotation to average noqa
      eig, vec = np.linalg.eig(M/N)

      return cls.fromQuaternion(np.real(vec.T[eig.argmax()]),acceptHomomorph = True)


    @classmethod
    def fromRandom(cls):
      r = np.random.random(3) 
      A = np.sqrt(r[2]) 
      B = np.sqrt(1.0-r[2]) 
      w = np.cos(2.0*np.pi*r[0])*A 
      x = np.sin(2.0*np.pi*r[1])*B 
      y = np.cos(2.0*np.pi*r[1])*B 
      z = np.sin(2.0*np.pi*r[0])*A 
      return cls.fromQuaternion([w,x,y,z],acceptHomomorph=True) 



# ******************************************************************************************
class Symmetry:
  """
  Symmetry operations for lattice systems
  
  https://en.wikipedia.org/wiki/Crystal_system
  """

  lattices = [None,'orthorhombic','tetragonal','hexagonal','cubic',]

  def __init__(self, symmetry = None):
    if isinstance(symmetry, str) and symmetry.lower() in Symmetry.lattices:
      self.lattice = symmetry.lower()
    else:
      self.lattice = None


  def __copy__(self):
    """Copy"""
    return self.__class__(self.lattice)

  copy = __copy__


  def __repr__(self):
    """Readable string"""
    return '{}'.format(self.lattice)


  def __eq__(self, other):
    """Equal to other"""
    return self.lattice == other.lattice

  def __neq__(self, other):
    """Not equal to other"""
    return not self.__eq__(other)

  def __cmp__(self,other):
    """Linear ordering"""
    myOrder    = Symmetry.lattices.index(self.lattice)
    otherOrder = Symmetry.lattices.index(other.lattice)
    return (myOrder > otherOrder) - (myOrder < otherOrder)

  def symmetryOperations(self,members=[]):
    """List (or single element) of symmetry operations as rotations."""
    if self.lattice == 'cubic':
      symQuats =  [
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
    elif self.lattice == 'hexagonal':
      symQuats =  [
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
    elif self.lattice == 'tetragonal':
      symQuats =  [
                    [ 1.0,            0.0,            0.0,            0.0            ],
                    [ 0.0,            1.0,            0.0,            0.0            ],
                    [ 0.0,            0.0,            1.0,            0.0            ],
                    [ 0.0,            0.0,            0.0,            1.0            ],
                    [ 0.0,            0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                    [ 0.0,           -0.5*np.sqrt(2), 0.5*np.sqrt(2), 0.0            ],
                    [ 0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                    [-0.5*np.sqrt(2), 0.0,            0.0,            0.5*np.sqrt(2) ],
                  ]
    elif self.lattice == 'orthorhombic':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,0.0,0.0,1.0 ],
                  ]
    else:
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                  ]

    symOps = list(map(Rotation,
                  np.array(symQuats)[np.atleast_1d(members) if members != [] else range(len(symQuats))]))
    try:
        iter(members)                                                                               # asking for (even empty) list of members?
    except TypeError:
        return symOps[0]                                                                            # no, return rotation object
    else:
        return symOps                                                                               # yes, return list of rotations
    

  def inFZ(self,rodrigues):
    """
    Check whether given Rodrigues vector falls into fundamental zone of own symmetry.
    
    Fundamental zone in Rodrigues space is point symmetric around origin.
    """
    if (len(rodrigues) != 3):
      raise ValueError('Input is not a Rodriques-Frank vector.\n')
    
    if np.any(rodrigues == np.inf): return False

    Rabs = abs(rodrigues)

    if self.lattice == 'cubic':
      return     math.sqrt(2.0)-1.0 >= Rabs[0] \
             and math.sqrt(2.0)-1.0 >= Rabs[1] \
             and math.sqrt(2.0)-1.0 >= Rabs[2] \
             and 1.0 >= Rabs[0] + Rabs[1] + Rabs[2]
    elif self.lattice == 'hexagonal':
      return     1.0 >= Rabs[0] and 1.0 >= Rabs[1] and 1.0 >= Rabs[2] \
             and 2.0 >= math.sqrt(3)*Rabs[0] + Rabs[1] \
             and 2.0 >= math.sqrt(3)*Rabs[1] + Rabs[0] \
             and 2.0 >= math.sqrt(3) + Rabs[2]
    elif self.lattice == 'tetragonal':
      return     1.0 >= Rabs[0] and 1.0 >= Rabs[1] \
             and math.sqrt(2.0) >= Rabs[0] + Rabs[1] \
             and math.sqrt(2.0) >= Rabs[2] + 1.0
    elif self.lattice == 'orthorhombic':
      return     1.0 >= Rabs[0] and 1.0 >= Rabs[1] and 1.0 >= Rabs[2]
    else:
      return True


  def inDisorientationSST(self,rodrigues):
    """
    Check whether given Rodrigues vector (of misorientation) falls into standard stereographic triangle of own symmetry.

    Determination of disorientations follow the work of A. Heinz and P. Neumann:
    Representation of Orientation and Disorientation Data for Cubic, Hexagonal, Tetragonal and Orthorhombic Crystals
    Acta Cryst. (1991). A47, 780-789
    """
    if (len(rodrigues) != 3):
      raise ValueError('Input is not a Rodriques-Frank vector.\n')
    R = rodrigues
    
    epsilon = 0.0
    if self.lattice == 'cubic':
      return R[0] >= R[1]+epsilon                and R[1] >= R[2]+epsilon    and R[2] >= epsilon

    elif self.lattice == 'hexagonal':
      return R[0] >= math.sqrt(3)*(R[1]-epsilon) and R[1] >= epsilon         and R[2] >= epsilon

    elif self.lattice == 'tetragonal':
      return R[0] >= R[1]-epsilon                and R[1] >= epsilon         and R[2] >= epsilon

    elif self.lattice == 'orthorhombic':
      return R[0] >= epsilon                     and R[1] >= epsilon         and R[2] >= epsilon

    else:
      return True


  def inSST(self,
            vector,
            proper = False,
            color = False):
    """
    Check whether given vector falls into standard stereographic triangle of own symmetry.

    proper considers only vectors with z >= 0, hence uses two neighboring SSTs.
    Return inverse pole figure color if requested.
    Bases are computed from

    basis = {'cubic' :        np.linalg.inv(np.array([[0.,0.,1.],                                   # direction of red
                                                      [1.,0.,1.]/np.sqrt(2.),                       # direction of green
                                                      [1.,1.,1.]/np.sqrt(3.)]).T),                  # direction of blue
             'hexagonal' :    np.linalg.inv(np.array([[0.,0.,1.],                                   # direction of red
                                                      [1.,0.,0.],                                   # direction of green
                                                      [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).T),         # direction of blue
             'tetragonal' :   np.linalg.inv(np.array([[0.,0.,1.],                                   # direction of red
                                                      [1.,0.,0.],                                   # direction of green
                                                      [1.,1.,0.]/np.sqrt(2.)]).T),                  # direction of blue
             'orthorhombic' : np.linalg.inv(np.array([[0.,0.,1.],                                   # direction of red
                                                      [1.,0.,0.],                                   # direction of green
                                                      [0.,1.,0.]]).T),                              # direction of blue
            }
    """
    if self.lattice == 'cubic':
      basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                     [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                     [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
                 'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                     [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                     [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
              }
    elif self.lattice == 'hexagonal':
      basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                     [ 1.            , -np.sqrt(3.)   ,  0. ],
                                     [ 0.            ,  2.            ,  0. ] ]),
               'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                     [-1.            ,  np.sqrt(3.)   ,  0. ],
                                     [ np.sqrt(3.)   , -1.            ,  0. ] ]),
              }
    elif self.lattice == 'tetragonal':
      basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                     [ 1.            , -1.            ,  0. ],
                                     [ 0.            ,  np.sqrt(2.)   ,  0. ] ]),
               'proper':np.array([   [ 0.            ,  0.            ,  1. ],
                                     [-1.            ,  1.            ,  0. ],
                                     [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
              }
    elif self.lattice == 'orthorhombic':
      basis = {'improper':np.array([ [ 0., 0., 1.],
                                     [ 1., 0., 0.],
                                     [ 0., 1., 0.] ]),
                 'proper':np.array([ [ 0., 0., 1.],
                                     [-1., 0., 0.],
                                     [ 0., 1., 0.] ]),
              }
    else:                                                                                           # direct exit for unspecified symmetry
      if color:
        return (True,np.zeros(3,'d'))
      else:
        return True

    v = np.array(vector,dtype=float)
    if proper:                                                                                      # check both improper ...
      theComponents = np.dot(basis['improper'],v)
      inSST = np.all(theComponents >= 0.0)
      if not inSST:                                                                                 # ... and proper SST
        theComponents = np.dot(basis['proper'],v)
        inSST = np.all(theComponents >= 0.0)
    else:      
      v[2] = abs(v[2])                                                                              # z component projects identical 
      theComponents = np.dot(basis['improper'],v)                                                   # for positive and negative values
      inSST = np.all(theComponents >= 0.0)

    if color:                                                                                       # have to return color array
      if inSST:
        rgb = np.power(theComponents/np.linalg.norm(theComponents),0.5)                             # smoothen color ramps
        rgb = np.minimum(np.ones(3,dtype=float),rgb)                                                # limit to maximum intensity
        rgb /= max(rgb)                                                                             # normalize to (HS)V = 1
      else:
        rgb = np.zeros(3,dtype=float)
      return (inSST,rgb)
    else:
      return inSST

# code derived from https://github.com/ezag/pyeuclid
# suggested reading: http://web.mit.edu/2.998/www/QuaternionReport1.pdf


# ******************************************************************************************
class Lattice:
  """
  Lattice system
  
  Currently, this contains only a mapping from Bravais lattice to symmetry
  and orientation relationships. It could include twin and slip systems.
  https://en.wikipedia.org/wiki/Bravais_lattice
  """

  lattices = {
              'triclinic':{'symmetry':None},
              'bct':{'symmetry':'tetragonal'},
              'hex':{'symmetry':'hexagonal'},
              'fcc':{'symmetry':'cubic','c/a':1.0},
              'bcc':{'symmetry':'cubic','c/a':1.0},
             }


  def __init__(self, lattice):
    self.lattice  = lattice
    self.symmetry = Symmetry(self.lattices[lattice]['symmetry'])
    
    
  def __repr__(self):
    """Report basic lattice information"""
    return 'Bravais lattice {} ({} symmetry)'.format(self.lattice,self.symmetry)
    
    
  # Kurdjomov--Sachs orientation relationship for fcc <-> bcc transformation
  # from S. Morito et al. Journal of Alloys and Compounds 577 (2013) 587-S592
  # also see K. Kitahara et al. Acta Materialia 54 (2006) 1279-1288
  KS = {'mapping':{'fcc':0,'bcc':1},
      'planes': np.array([
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[  1,  1, -1],[  0,  1,  1]],
      [[  1,  1, -1],[  0,  1,  1]],
      [[  1,  1, -1],[  0,  1,  1]],
      [[  1,  1, -1],[  0,  1,  1]],
      [[  1,  1, -1],[  0,  1,  1]],
      [[  1,  1, -1],[  0,  1,  1]]],dtype='float'),
      'directions': np.array([     
      [[ -1,  0,  1],[ -1, -1,  1]],
      [[ -1,  0,  1],[ -1,  1, -1]],
      [[  0,  1, -1],[ -1, -1,  1]],
      [[  0,  1, -1],[ -1,  1, -1]],
      [[  1, -1,  0],[ -1, -1,  1]],
      [[  1, -1,  0],[ -1,  1, -1]],
      [[  1,  0, -1],[ -1, -1,  1]],
      [[  1,  0, -1],[ -1,  1, -1]],
      [[ -1, -1,  0],[ -1, -1,  1]],
      [[ -1, -1,  0],[ -1,  1, -1]],
      [[  0,  1,  1],[ -1, -1,  1]],
      [[  0,  1,  1],[ -1,  1, -1]],
      [[  0, -1,  1],[ -1, -1,  1]],
      [[  0, -1,  1],[ -1,  1, -1]],
      [[ -1,  0, -1],[ -1, -1,  1]],
      [[ -1,  0, -1],[ -1,  1, -1]],
      [[  1,  1,  0],[ -1, -1,  1]],
      [[  1,  1,  0],[ -1,  1, -1]],
      [[ -1,  1,  0],[ -1, -1,  1]],
      [[ -1,  1,  0],[ -1,  1, -1]],
      [[  0, -1, -1],[ -1, -1,  1]],
      [[  0, -1, -1],[ -1,  1, -1]],
      [[  1,  0,  1],[ -1, -1,  1]],
      [[  1,  0,  1],[ -1,  1, -1]]],dtype='float')}
      
  # Greninger--Troiano orientation relationship for fcc <-> bcc transformation
  # from Y. He et al. Journal of Applied Crystallography 39 (2006) 72-81
  GT = {'mapping':{'fcc':0,'bcc':1},
      'planes': np.array([
      [[  1,  1,  1],[  1,  0,  1]],
      [[  1,  1,  1],[  1,  1,  0]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[ -1, -1,  1],[ -1,  0,  1]],
      [[ -1, -1,  1],[ -1, -1,  0]],
      [[ -1, -1,  1],[  0, -1,  1]],
      [[ -1,  1,  1],[ -1,  0,  1]],
      [[ -1,  1,  1],[ -1,  1,  0]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  1,  0,  1]],
      [[  1, -1,  1],[  1, -1,  0]],
      [[  1, -1,  1],[  0, -1,  1]],
      [[  1,  1,  1],[  1,  1,  0]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  1,  0,  1]],
      [[ -1, -1,  1],[ -1, -1,  0]],
      [[ -1, -1,  1],[  0, -1,  1]],
      [[ -1, -1,  1],[ -1,  0,  1]],
      [[ -1,  1,  1],[ -1,  1,  0]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[ -1,  0,  1]],
      [[  1, -1,  1],[  1, -1,  0]],
      [[  1, -1,  1],[  0, -1,  1]],
      [[  1, -1,  1],[  1,  0,  1]]],dtype='float'),
      'directions': np.array([
      [[ -5,-12, 17],[-17, -7, 17]],
      [[ 17, -5,-12],[ 17,-17, -7]],
      [[-12, 17, -5],[ -7, 17,-17]],
      [[  5, 12, 17],[ 17,  7, 17]],
      [[-17,  5,-12],[-17, 17, -7]],
      [[ 12,-17, -5],[  7,-17,-17]],
      [[ -5, 12,-17],[-17,  7,-17]],
      [[ 17,  5, 12],[ 17, 17,  7]],
      [[-12,-17,  5],[ -7,-17, 17]],
      [[  5,-12,-17],[ 17, -7,-17]],
      [[-17, -5, 12],[-17,-17,  7]],
      [[ 12, 17,  5],[  7, 17, 17]],
      [[ -5, 17,-12],[-17, 17, -7]],
      [[-12, -5, 17],[ -7,-17, 17]],
      [[ 17,-12, -5],[ 17, -7,-17]],
      [[  5,-17,-12],[ 17,-17, -7]],
      [[ 12,  5, 17],[  7, 17, 17]],
      [[-17, 12, -5],[-17,  7,-17]],
      [[ -5,-17, 12],[-17,-17,  7]],
      [[-12,  5,-17],[ -7, 17,-17]],
      [[ 17, 12,  5],[ 17,  7, 17]],
      [[  5, 17, 12],[ 17, 17,  7]],
      [[ 12, -5,-17],[  7,-17,-17]],
      [[-17,-12,  5],[-17,  7, 17]]],dtype='float')}
 
  # Greninger--Troiano' orientation relationship for fcc <-> bcc transformation
  # from Y. He et al. Journal of Applied Crystallography 39 (2006) 72-81
  GTprime = {'mapping':{'fcc':0,'bcc':1},
      'planes': np.array([
      [[  7, 17, 17],[ 12,  5, 17]],
      [[ 17,  7, 17],[ 17, 12,  5]],
      [[ 17, 17,  7],[  5, 17, 12]],
      [[ -7,-17, 17],[-12, -5, 17]],
      [[-17, -7, 17],[-17,-12,  5]],
      [[-17,-17,  7],[ -5,-17, 12]],
      [[  7,-17,-17],[ 12, -5,-17]],
      [[ 17, -7,-17],[ 17,-12, -5]],
      [[ 17,-17, -7],[  5,-17,-12]],
      [[ -7, 17,-17],[-12,  5,-17]],
      [[-17,  7,-17],[-17, 12, -5]],
      [[-17, 17, -7],[ -5, 17,-12]],
      [[  7, 17, 17],[ 12, 17,  5]],
      [[ 17,  7, 17],[  5, 12, 17]],
      [[ 17, 17,  7],[ 17,  5, 12]],
      [[ -7,-17, 17],[-12,-17,  5]],
      [[-17, -7, 17],[ -5,-12, 17]],
      [[-17,-17,  7],[-17, -5, 12]],
      [[  7,-17,-17],[ 12,-17, -5]],
      [[ 17, -7,-17],[ 5, -12,-17]],
      [[ 17,-17,  7],[ 17, -5,-12]],
      [[ -7, 17,-17],[-12, 17, -5]],
      [[-17,  7,-17],[ -5, 12,-17]],
      [[-17, 17, -7],[-17,  5,-12]]],dtype='float'),
      'directions': np.array([
      [[  0,  1, -1],[  1,  1, -1]],
      [[ -1,  0,  1],[ -1,  1,  1]],
      [[  1, -1,  0],[  1, -1,  1]],
      [[  0, -1, -1],[ -1, -1, -1]],
      [[  1,  0,  1],[  1, -1,  1]],
      [[  1, -1,  0],[  1, -1, -1]],
      [[  0,  1, -1],[ -1,  1, -1]],
      [[  1,  0,  1],[  1,  1,  1]],
      [[ -1, -1,  0],[ -1, -1,  1]],
      [[  0, -1, -1],[  1, -1, -1]],
      [[ -1,  0,  1],[ -1, -1,  1]],
      [[ -1, -1,  0],[ -1, -1, -1]],
      [[  0, -1,  1],[  1, -1,  1]],
      [[  1,  0, -1],[  1,  1, -1]],
      [[ -1,  1,  0],[ -1,  1,  1]],
      [[  0,  1,  1],[ -1,  1,  1]],
      [[ -1,  0, -1],[ -1, -1, -1]],
      [[ -1,  1,  0],[ -1,  1, -1]],
      [[  0, -1,  1],[ -1, -1,  1]],
      [[ -1,  0, -1],[ -1,  1, -1]],
      [[  1,  1,  0],[  1,  1,  1]],
      [[  0,  1,  1],[  1,  1,  1]],
      [[  1,  0, -1],[  1, -1, -1]],
      [[  1,  1,  0],[  1,  1, -1]]],dtype='float')}
    
  # Nishiyama--Wassermann orientation relationship for fcc <-> bcc transformation
  # from H. Kitahara et al. Materials Characterization 54 (2005) 378-386
  NW = {'mapping':{'fcc':0,'bcc':1},
      'planes': np.array([
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[  1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[ -1,  1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[  1, -1,  1],[  0,  1,  1]],
      [[ -1, -1,  1],[  0,  1,  1]],
      [[ -1, -1,  1],[  0,  1,  1]],
      [[ -1, -1,  1],[  0,  1,  1]]],dtype='float'),
      'directions': np.array([
      [[  2, -1, -1],[  0, -1,  1]],
      [[ -1,  2, -1],[  0, -1,  1]],
      [[ -1, -1,  2],[  0, -1,  1]],
      [[ -2, -1, -1],[  0, -1,  1]],
      [[  1,  2, -1],[  0, -1,  1]],
      [[  1, -1,  2],[  0, -1,  1]],
      [[  2,  1, -1],[  0, -1,  1]],
      [[ -1, -2, -1],[  0, -1,  1]],
      [[ -1,  1,  2],[  0, -1,  1]],
      [[ -1,  2,  1],[  0, -1,  1]],
      [[ -1,  2,  1],[  0, -1,  1]],
      [[ -1, -1, -2],[  0, -1,  1]]],dtype='float')}
               
  # Pitsch orientation relationship for fcc <-> bcc transformation              
  # from Y. He et al. Acta Materialia 53 (2005) 1179-1190               
  Pitsch = {'mapping':{'fcc':0,'bcc':1},
      'planes': np.array([
      [[  0,  1,  0],[ -1,  0,  1]],
      [[  0,  0,  1],[  1, -1,  0]],
      [[  1,  0,  0],[  0,  1, -1]],
      [[  1,  0,  0],[  0, -1, -1]],
      [[  0,  1,  0],[ -1,  0, -1]],
      [[  0,  0,  1],[ -1, -1,  0]],
      [[  0,  1,  0],[ -1,  0, -1]],
      [[  0,  0,  1],[ -1, -1,  0]],
      [[  1,  0,  0],[  0, -1, -1]],
      [[  1,  0,  0],[  0, -1,  1]],
      [[  0,  1,  0],[  1,  0, -1]],
      [[  0,  0,  1],[ -1,  1,  0]]],dtype='float'),
      'directions': np.array([
      [[  1,  0,  1],[  1, -1,  1]],
      [[  1,  1,  0],[  1,  1, -1]],
      [[  0,  1,  1],[ -1,  1,  1]],
      [[  0,  1, -1],[ -1,  1, -1]],
      [[ -1,  0,  1],[ -1, -1,  1]],
      [[  1, -1,  0],[  1, -1, -1]],
      [[  1,  0, -1],[  1, -1, -1]],
      [[ -1,  1,  0],[ -1,  1, -1]],
      [[  0, -1,  1],[ -1, -1,  1]],
      [[  0,  1,  1],[ -1,  1,  1]],
      [[  1,  0,  1],[  1, -1,  1]],
      [[  1,  1,  0],[  1,  1, -1]]],dtype='float')}
                              
  # Bain orientation relationship for fcc <-> bcc transformation                        
  # from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
  Bain = {'mapping':{'fcc':0,'bcc':1},
      'planes': np.array([
      [[  1,  0,  0],[  1,  0,  0]],
      [[  0,  1,  0],[  0,  1,  0]],
      [[  0,  0,  1],[  0,  0,  1]]],dtype='float'),
      'directions': np.array([
      [[  0,  1,  0],[  0,  1,  1]],
      [[  0,  0,  1],[  1,  0,  1]],
      [[  1,  0,  0],[  1,  1,  0]]],dtype='float')}
  
  def relationOperations(self,model):
    
    models={'KS':self.KS, 'GT':self.GT, "GT'":self.GTprime, 
            'NW':self.NW, 'Pitsch': self.Pitsch, 'Bain':self.Bain}
    try:
      relationship = models[model]
    except:
      raise KeyError('Orientation relationship "{}" is unknown'.format(model))

    if self.lattice not in relationship['mapping']:
      raise ValueError('Relationship "{}" not supported for lattice "{}"'.format(model,self.lattice))
      
    r = {'lattice':Lattice((set(relationship['mapping'])-{self.lattice}).pop()),                    # target lattice
         'rotations':[] }

    myPlane_id    = relationship['mapping'][self.lattice]
    otherPlane_id = (myPlane_id+1)%2
    myDir_id      = myPlane_id +2
    otherDir_id   = otherPlane_id +2

    for miller in np.hstack((relationship['planes'],relationship['directions'])):
      myPlane     = miller[myPlane_id]/    np.linalg.norm(miller[myPlane_id])
      myDir       = miller[myDir_id]/      np.linalg.norm(miller[myDir_id])
      myMatrix    = np.array([myDir,np.cross(myPlane,myDir),myPlane]).T

      otherPlane  = miller[otherPlane_id]/ np.linalg.norm(miller[otherPlane_id])
      otherDir    = miller[otherDir_id]/   np.linalg.norm(miller[otherDir_id])
      otherMatrix = np.array([otherDir,np.cross(otherPlane,otherDir),otherPlane]).T

      r['rotations'].append(Rotation.fromMatrix(np.dot(otherMatrix,myMatrix.T)))

    return r



class Orientation:
  """
  Crystallographic orientation
  
  A crystallographic orientation contains a rotation and a lattice
  """

  __slots__ = ['rotation','lattice']
  
  def __repr__(self):
    """Report lattice type and orientation"""
    return self.lattice.__repr__()+'\n'+self.rotation.__repr__()

  def __init__(self, rotation, lattice):
  
    if isinstance(lattice, Lattice):
      self.lattice = lattice
    else:
      self.lattice = Lattice(lattice)      # assume string
      
    if isinstance(rotation, Rotation):
      self.rotation = rotation
    else:
      self.rotation = Rotation(rotation)   # assume quaternion
      
  def disorientation(self,
                     other,
                     SST = True,
                     symmetries = False):
    """
    Disorientation between myself and given other orientation.

    Rotation axis falls into SST if SST == True.
    (Currently requires same symmetry for both orientations.
     Look into A. Heinz and P. Neumann 1991 for cases with differing sym.)
    """
    if self.lattice.symmetry != other.lattice.symmetry:
      raise NotImplementedError('disorientation between different symmetry classes not supported yet.')

    mySymEqs    =  self.equivalentOrientations() if SST else self.equivalentOrientations([0])     # take all or only first sym operation
    otherSymEqs = other.equivalentOrientations()
    
    for i,sA in enumerate(mySymEqs):
      aInv = sA.rotation.inversed()
      for j,sB in enumerate(otherSymEqs):
        b = sB.rotation
        r = b*aInv
        for k in range(2):
          r.inverse()
          breaker = self.lattice.symmetry.inFZ(r.asRodrigues(vector=True)) \
                    and (not SST or other.lattice.symmetry.inDisorientationSST(r.asRodrigues(vector=True)))
          if breaker: break
        if breaker: break
      if breaker: break

    return (Orientation(r,self.lattice), i,j, k == 1) if symmetries else r                          # disorientation ...
                                                                                                    # ... own sym, other sym,
                                                                                                    # self-->other: True, self<--other: False


  def inFZ(self):
    return self.lattice.symmetry.inFZ(self.rotation.asRodrigues(vector=True))
  
  def equivalentOrientations(self,members=[]):
    """List of orientations which are symmetrically equivalent"""
    try:
      iter(members)                                                                                 # asking for (even empty) list of members?
    except TypeError:
      return self.__class__(self.lattice.symmetry.symmetryOperations(members)*self.rotation,self.lattice) # no, return rotation object
    else:
      return [self.__class__(q*self.rotation,self.lattice) \
                                    for q in self.lattice.symmetry.symmetryOperations(members)]           # yes, return list of rotations
    
  def relatedOrientations(self,model):
    """List of orientations related by the given orientation relationship"""
    r = self.lattice.relationOperations(model)
    return [self.__class__(self.rotation*o,r['lattice']) for o in r['rotations']]
    
  def reduced(self):
    """Transform orientation to fall into fundamental zone according to symmetry"""
    for me in self.equivalentOrientations():
      if self.lattice.symmetry.inFZ(me.rotation.asRodrigues(vector=True)): break

    return self.__class__(me.rotation,self.lattice)
    
  def inversePole(self,
                  axis,
                  proper = False,
                  SST = True):
    """Axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)"""
    if SST:                                                                                         # pole requested to be within SST
      for i,o in enumerate(self.equivalentOrientations()):                                          # test all symmetric equivalent quaternions
        pole = o.rotation*axis                                                                      # align crystal direction to axis
        if self.lattice.symmetry.inSST(pole,proper): break                                          # found SST version
    else:
      pole = self.rotation*axis                                                                     # align crystal direction to axis

    return (pole,i if SST else 0)
    
    
  def IPFcolor(self,axis):
    """TSL color of inverse pole figure for given axis"""
    color = np.zeros(3,'d')

    for o in self.equivalentOrientations():
      pole = o.rotation*axis                                                                        # align crystal direction to axis
      inSST,color = self.lattice.symmetry.inSST(pole,color=True)
      if inSST: break

    return color                                                                        


  @classmethod
  def fromAverage(cls,
                  orientations,
                  weights = []):
    """Create orientation from average of list of orientations"""
    if not all(isinstance(item, Orientation) for item in orientations):
      raise TypeError("Only instances of Orientation can be averaged.")

    closest = []
    ref = orientations[0]
    for o in orientations:
      closest.append(o.equivalentOrientations(
                    ref.disorientation(o,
                                       SST = False,                                                 # select (o[ther]'s) sym orientation
                                       symmetries = True)[2]).rotation)                             # with lowest misorientation

    return Orientation(Rotation.fromAverage(closest,weights),ref.lattice)


  def average(self,other):
    """Calculate the average rotation"""
    return Orientation.fromAverage([self,other])


####################################################################################################
# Code below available according to the following conditions on https://github.com/MarDiehl/3Drotations
####################################################################################################
# Copyright (c) 2017-2019, Martin Diehl/Max-Planck-Institut für Eisenforschung GmbH
# Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are 
# permitted provided that the following conditions are met:
#
#     - Redistributions of source code must retain the above copyright notice, this list 
#        of conditions and the following disclaimer.
#     - Redistributions in binary form must reproduce the above copyright notice, this 
#        list of conditions and the following disclaimer in the documentation and/or 
#        other materials provided with the distribution.
#     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
#        of its contributors may be used to endorse or promote products derived from 
#        this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
####################################################################################################

def isone(a):
  return np.isclose(a,1.0,atol=1.0e-7,rtol=0.0)
  
def iszero(a):
  return np.isclose(a,0.0,atol=1.0e-12,rtol=0.0)
  
#---------- quaternion ----------

def qu2om(qu):
  """Quaternion to orientation matrix"""
  qq = qu[0]**2-(qu[1]**2 + qu[2]**2 + qu[3]**2)
  om = np.diag(qq + 2.0*np.array([qu[1],qu[2],qu[3]])**2)

  om[1,0] = 2.0*(qu[2]*qu[1]+qu[0]*qu[3])
  om[0,1] = 2.0*(qu[1]*qu[2]-qu[0]*qu[3])
  om[2,1] = 2.0*(qu[3]*qu[2]+qu[0]*qu[1])
  om[1,2] = 2.0*(qu[2]*qu[3]-qu[0]*qu[1])
  om[0,2] = 2.0*(qu[1]*qu[3]+qu[0]*qu[2])
  om[2,0] = 2.0*(qu[3]*qu[1]-qu[0]*qu[2])
  return om if P > 0.0 else om.T
  
  
def qu2eu(qu):
  """Quaternion to Euler angles""" 
  q03 = qu[0]**2+qu[3]**2
  q12 = qu[1]**2+qu[2]**2
  chi = np.sqrt(q03*q12)
  
  if iszero(chi):
    eu = np.array([np.arctan2(-P*2.0*qu[0]*qu[3],qu[0]**2-qu[3]**2), 0.0,   0.0]) if iszero(q12) else \
         np.array([np.arctan2(2.0*qu[1]*qu[2],qu[1]**2-qu[2]**2),         np.pi, 0.0])
  else:
    eu = np.array([np.arctan2((-P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-P*qu[0]*qu[1]-qu[2]*qu[3])*chi ),
                   np.arctan2( 2.0*chi, q03-q12 ), 
                   np.arctan2(( P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-P*qu[0]*qu[1]+qu[2]*qu[3])*chi )])
  
  # reduce Euler angles to definition range, i.e a lower limit of 0.0
  eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
  return eu


def qu2ax(qu):
  """
  Quaternion to axis angle
  
  Modified version of the original formulation, should be numerically more stable
  """
  if iszero(qu[1]**2+qu[2]**2+qu[3]**2):                                                            # set axis to [001] if the angle is 0/360
    ax = [ 0.0, 0.0, 1.0, 0.0 ]
  elif not iszero(qu[0]):
    s = np.sign(qu[0])/np.sqrt(qu[1]**2+qu[2]**2+qu[3]**2)
    omega = 2.0 * np.arccos(np.clip(qu[0],-1.0,1.0))
    ax = [ qu[1]*s, qu[2]*s, qu[3]*s, omega ]
  else:
    ax = [ qu[1], qu[2], qu[3], np.pi]

  return np.array(ax)


def qu2ro(qu):
  """Quaternion to Rodrigues vector"""
  if iszero(qu[0]):
    ro = [qu[1], qu[2], qu[3], np.inf]
  else:
    s = np.linalg.norm([qu[1],qu[2],qu[3]])
    ro = [0.0,0.0,P,0.0] if iszero(s) else \
         [ qu[1]/s,  qu[2]/s,  qu[3]/s, np.tan(np.arccos(np.clip(qu[0],-1.0,1.0)))]                # avoid numerical difficulties
 
  return np.array(ro)


def qu2ho(qu):
  """Quaternion to homochoric"""
  omega = 2.0 * np.arccos(np.clip(qu[0],-1.0,1.0))                                                  # avoid numerical difficulties
  
  if iszero(omega):
    ho = np.array([ 0.0, 0.0, 0.0 ])
  else:
    ho = np.array([qu[1], qu[2], qu[3]])
    f  = 0.75 * ( omega - np.sin(omega) )
    ho = ho/np.linalg.norm(ho) * f**(1./3.)

  return ho


def qu2cu(qu):
  """Quaternion to cubochoric"""
  return ho2cu(qu2ho(qu))


#---------- orientation matrix ----------

def om2qu(om):
  """
  Orientation matrix to quaternion
  
  The original formulation (direct conversion) had (numerical?) issues
  """
  return eu2qu(om2eu(om))


def om2eu(om):
  """Orientation matrix to Euler angles"""
  if abs(om[2,2]) < 1.0:
    zeta = 1.0/np.sqrt(1.0-om[2,2]**2)
    eu = np.array([np.arctan2(om[2,0]*zeta,-om[2,1]*zeta),
                   np.arccos(om[2,2]),
                   np.arctan2(om[0,2]*zeta, om[1,2]*zeta)])
  else:
    eu = np.array([np.arctan2( om[0,1],om[0,0]), np.pi*0.5*(1-om[2,2]),0.0])                        # following the paper, not the reference implementation
  
  # reduce Euler angles to definition range, i.e a lower limit of 0.0
  eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
  return eu


def om2ax(om):
  """Orientation matrix to axis angle"""
  ax=np.empty(4)

  # first get the rotation angle
  t = 0.5*(om.trace() -1.0)
  ax[3] = np.arccos(np.clip(t,-1.0,1.0))
  
  if iszero(ax[3]):
    ax = [ 0.0, 0.0, 1.0, 0.0]
  else:
    w,vr = np.linalg.eig(om)
  # next, find the eigenvalue (1,0j)
    i = np.where(np.isclose(w,1.0+0.0j))[0][0]
    ax[0:3] = np.real(vr[0:3,i])
    diagDelta = np.array([om[1,2]-om[2,1],om[2,0]-om[0,2],om[0,1]-om[1,0]])
    ax[0:3] = np.where(iszero(diagDelta), ax[0:3],np.abs(ax[0:3])*np.sign(-P*diagDelta))
  
  return np.array(ax)
  

def om2ro(om):
  """Orientation matrix to Rodriques vector"""
  return eu2ro(om2eu(om))
  

def om2ho(om):
  """Orientation matrix to homochoric"""
  return ax2ho(om2ax(om))


def om2cu(om):
  """Orientation matrix to cubochoric"""
  return ho2cu(om2ho(om))
  
  
#---------- Euler angles ----------

def eu2qu(eu):
  """Euler angles to quaternion"""
  ee = 0.5*eu
  cPhi = np.cos(ee[1])
  sPhi = np.sin(ee[1])
  qu = np.array([     cPhi*np.cos(ee[0]+ee[2]),
                   -P*sPhi*np.cos(ee[0]-ee[2]),
                   -P*sPhi*np.sin(ee[0]-ee[2]),
                   -P*cPhi*np.sin(ee[0]+ee[2]) ])
  if qu[0] < 0.0: qu*=-1
  return qu


def eu2om(eu):
  """Euler angles to orientation matrix"""
  c = np.cos(eu)
  s = np.sin(eu)

  om = np.array([[+c[0]*c[2]-s[0]*s[2]*c[1], +s[0]*c[2]+c[0]*s[2]*c[1], +s[2]*s[1]],
                 [-c[0]*s[2]-s[0]*c[2]*c[1], -s[0]*s[2]+c[0]*c[2]*c[1], +c[2]*s[1]],
                 [+s[0]*s[1],                -c[0]*s[1],                +c[1]     ]])

  om[np.where(iszero(om))] = 0.0
  return om


def eu2ax(eu):
  """Euler angles to axis angle""" 
  t = np.tan(eu[1]*0.5)
  sigma = 0.5*(eu[0]+eu[2])
  delta = 0.5*(eu[0]-eu[2])
  tau   = np.linalg.norm([t,np.sin(sigma)])
  alpha = np.pi if iszero(np.cos(sigma)) else \
          2.0*np.arctan(tau/np.cos(sigma))
  
  if iszero(alpha):
    ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
  else:
    ax = -P/tau * np.array([ t*np.cos(delta), t*np.sin(delta), np.sin(sigma) ])                     # passive axis-angle pair so a minus sign in front
    ax = np.append(ax,alpha) 
    if alpha < 0.0: ax *= -1.0                                                                      # ensure alpha is positive
  
  return ax


def eu2ro(eu):
  """Euler angles to Rodrigues vector"""
  ro = eu2ax(eu)                                                                                    # convert to axis angle representation
  if ro[3] >= np.pi:                                                                                # Differs from original implementation. check convention 5
    ro[3] = np.inf
  elif iszero(ro[3]):
    ro = np.array([ 0.0, 0.0, P, 0.0 ])
  else:
    ro[3] = np.tan(ro[3]*0.5)
  
  return ro


def eu2ho(eu):
  """Euler angles to homochoric"""
  return ax2ho(eu2ax(eu))


def eu2cu(eu):
  """Euler angles to cubochoric"""
  return ho2cu(eu2ho(eu))


#---------- axis angle ----------

def ax2qu(ax):
  """Axis angle to quaternion"""
  if iszero(ax[3]):
    qu = np.array([ 1.0, 0.0, 0.0, 0.0 ])
  else:
    c = np.cos(ax[3]*0.5)
    s = np.sin(ax[3]*0.5)
    qu = np.array([ c, ax[0]*s, ax[1]*s, ax[2]*s ])
  
  return qu


def ax2om(ax):
  """Axis angle to orientation matrix"""
  c = np.cos(ax[3])
  s = np.sin(ax[3])
  omc = 1.0-c
  om=np.diag(ax[0:3]**2*omc + c)

  for idx in [[0,1,2],[1,2,0],[2,0,1]]:
    q = omc*ax[idx[0]] * ax[idx[1]]
    om[idx[0],idx[1]] = q + s*ax[idx[2]]
    om[idx[1],idx[0]] = q - s*ax[idx[2]]

  return om if P < 0.0 else om.T

  
def ax2eu(ax):
  """Orientation matrix to Euler angles"""
  return om2eu(ax2om(ax))


def ax2ro(ax):
  """Axis angle to Rodrigues vector""" 
  if iszero(ax[3]):
    ro = [ 0.0, 0.0, P, 0.0 ]
  else:
    ro = [ax[0], ax[1], ax[2]]
    # 180 degree case
    ro += [np.inf] if np.isclose(ax[3],np.pi,atol=1.0e-15,rtol=0.0) else \
          [np.tan(ax[3]*0.5)]

  return np.array(ro)


def ax2ho(ax):
  """Axis angle to homochoric""" 
  f = (0.75 * ( ax[3] - np.sin(ax[3]) ))**(1.0/3.0)
  ho = ax[0:3] * f
  return ho


def ax2cu(ax):
  """Axis angle to cubochoric"""
  return ho2cu(ax2ho(ax))


#---------- Rodrigues--Frank ----------

def ro2qu(ro):
  """Rodrigues vector to quaternion"""
  return ax2qu(ro2ax(ro))


def ro2om(ro):
 """Rodgrigues vector to orientation matrix"""
 return ax2om(ro2ax(ro))


def ro2eu(ro):
  """Rodrigues vector to orientation matrix"""
  return om2eu(ro2om(ro))
 
 
def ro2ax(ro):
  """Rodrigues vector to axis angle"""
  ta = ro[3]
  
  if iszero(ta):
    ax = [ 0.0, 0.0, 1.0, 0.0 ]
  elif not np.isfinite(ta):
    ax = [ ro[0], ro[1], ro[2], np.pi ]
  else:
    angle = 2.0*np.arctan(ta)
    ta = 1.0/np.linalg.norm(ro[0:3])
    ax = [ ro[0]/ta, ro[1]/ta, ro[2]/ta, angle ]

  return np.array(ax)


def ro2ho(ro):
  """Rodrigues vector to homochoric""" 
  if iszero(np.sum(ro[0:3]**2.0)):
    ho = [ 0.0, 0.0, 0.0 ]
  else:
    f = 2.0*np.arctan(ro[3]) -np.sin(2.0*np.arctan(ro[3])) if np.isfinite(ro[3]) else np.pi
    ho = ro[0:3] * (0.75*f)**(1.0/3.0)

  return np.array(ho)


def ro2cu(ro):
  """Rodrigues vector to cubochoric"""
  return ho2cu(ro2ho(ro))


#---------- homochoric ----------

def ho2qu(ho):
  """Homochoric to quaternion"""
  return ax2qu(ho2ax(ho))


def ho2om(ho):
  """Homochoric to orientation matrix"""
  return ax2om(ho2ax(ho))


def ho2eu(ho):
  """Homochoric to Euler angles"""
  return ax2eu(ho2ax(ho))
  

def ho2ax(ho):
  """Homochoric to axis angle"""
  tfit = np.array([+1.0000000000018852,      -0.5000000002194847,
                   -0.024999992127593126,    -0.003928701544781374,
                   -0.0008152701535450438,   -0.0002009500426119712,
                   -0.00002397986776071756,  -0.00008202868926605841,
                   +0.00012448715042090092,  -0.0001749114214822577,
                   +0.0001703481934140054,   -0.00012062065004116828,
                   +0.000059719705868660826, -0.00001980756723965647,
                   +0.000003953714684212874, -0.00000036555001439719544])
  # normalize h and store the magnitude
  hmag_squared = np.sum(ho**2.)
  if iszero(hmag_squared):
    ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
  else:
    hm = hmag_squared
  
  # convert the magnitude to the rotation angle
    s = tfit[0] + tfit[1] * hmag_squared
    for i in range(2,16):
      hm *= hmag_squared
      s  += tfit[i] * hm
    ax = np.append(ho/np.sqrt(hmag_squared),2.0*np.arccos(np.clip(s,-1.0,1.0)))
  return ax


def ho2ro(ho):
  """Axis angle to Rodriques vector"""
  return ax2ro(ho2ax(ho))


def ho2cu(ho):
  """Homochoric to cubochoric"""
  return  Lambert.BallToCube(ho)


#---------- cubochoric ----------

def cu2qu(cu):
  """Cubochoric to quaternion"""
  return ho2qu(cu2ho(cu))


def cu2om(cu):
  """Cubochoric to orientation matrix"""
  return ho2om(cu2ho(cu))


def cu2eu(cu):
  """Cubochoric to Euler angles"""
  return ho2eu(cu2ho(cu))


def cu2ax(cu):
  """Cubochoric to axis angle"""
  return ho2ax(cu2ho(cu))


def cu2ro(cu):
  """Cubochoric to Rodrigues vector"""
  return ho2ro(cu2ho(cu))


def cu2ho(cu):
  """Cubochoric to homochoric"""
  return  Lambert.CubeToBall(cu)
