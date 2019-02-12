# -*- coding: UTF-8 no BOM -*-

import math,os
import numpy as np
from . import Lambert

P = -1

####################################################################################################
class Quaternion2:
    u"""
    Quaternion with basic operations
    
    q is the real part, p = (x, y, z) are the imaginary parts.
    Defintion of multiplication depends on variable P, P ∉ {-1,1}.
    """

    def __init__(self,
                 q    = 0.0,
                 p    = np.zeros(3,dtype=float)):
      """Initializes to identity unless specified"""
      self.q = q
      self.p = np.array(p)


    def __copy__(self):
      """Copy"""
      return self.__class__(q=self.q,
                            p=self.p.copy())

    copy = __copy__


    def __iter__(self):
      """Components"""
      return iter(self.asList())

    def asArray(self):
      """As numpy array"""
      return np.array((self.q,self.p[0],self.p[1],self.p[2]))
 
    def asList(self):
      return [self.q]+list(self.p)


    def __repr__(self):
      """Readable string"""
      return 'Quaternion(real={q:+.6f}, imag=<{p[0]:+.6f}, {p[1]:+.6f}, {p[2]:+.6f}>)'.format(q=self.q,p=self.p)
      

    def __add__(self, other):
      """Addition"""
      if isinstance(other, Quaternion2):
        return self.__class__(q=self.q + other.q,
                              p=self.p + other.p)
      else:
        return NotImplemented

    def __iadd__(self, other):
      """In-place addition"""
      if isinstance(other, Quaternion2):
        self.q += other.q
        self.p += other.p
        return self
      else:
        return NotImplemented
      
    def __pos__(self):                                                                      
      """Unary positive operator"""                                                                       
      return self
      
      
    def __sub__(self, other):
      """Subtraction"""
      if isinstance(other, Quaternion2):
        return self.__class__(q=self.q - other.q,
                              p=self.p - other.p)
      else:
        return NotImplemented

    def __isub__(self, other):
      """In-place subtraction"""
      if isinstance(other, Quaternion2):
        self.q -= other.q
        self.p -= other.p
        return self
      else:
        return NotImplemented
        
    def __neg__(self):                                                                      
      """Unary positive operator"""  
      self.q *= -1.0
      self.p *= -1.0                                                                     
      return self
      
      
    def __mul__(self, other):
      """Multiplication with quaternion or scalar"""
      if isinstance(other, Quaternion2):
        return self.__class__(q=self.q*other.q - np.dot(self.p,other.p),
                              p=self.q*other.p + other.q*self.p + P * np.cross(self.p,other.p))
      elif isinstance(other, (int, float)):
        return self.__class__(q=self.q*other,
                              p=self.p*other)
      else:
        return NotImplemented

    def __imul__(self, other):
      """In-place multiplication with quaternion or scalar"""
      if isinstance(other, Quaternion2):
        self.q = self.q*other.q - np.dot(self.p,other.p)
        self.p = self.q*other.p + other.q*self.p + P * np.cross(self.p,other.p)
        return self
      elif isinstance(other, (int, float)):
        self *= other
        return self
      else:
        return NotImplemented
      

    def __truediv__(self, other):
      """Divsion with quaternion or scalar"""
      if isinstance(other, Quaternion2):
        s = other.conjugate()/abs(other)**2.
        return self.__class__(q=self.q * s,
                              p=self.p * s)
      elif isinstance(other, (int, float)):
        self.q /= other
        self.p /= other
        return self
      else:
        return NotImplemented

    def __itruediv__(self, other):
      """In-place divsion with quaternion or scalar"""
      if isinstance(other, Quaternion2):
        s = other.conjugate()/abs(other)**2.
        self *= s
        return self
      elif isinstance(other, (int, float)):
        self.q /= other
        return self
      else:
        return NotImplemented


    def __pow__(self, exponent):
      """Power"""
      if isinstance(exponent, (int, float)):
        omega = np.acos(self.q)
        return self.__class__(q=         np.cos(exponent*omega),
                              p=self.p * np.sin(exponent*omega)/np.sin(omega))
      else:
        return NotImplemented

    def __ipow__(self, exponent):
      """In-place power"""
      if isinstance(exponent, (int, float)):
        omega = np.acos(self.q)
        self.q  = np.cos(exponent*omega)
        self.p *= np.sin(exponent*omega)/np.sin(omega)
      else:
        return NotImplemented
      

    def __abs__(self):
      """Norm"""
      return math.sqrt(self.q ** 2 + np.dot(self.p,self.p))

    magnitude = __abs__


    def __eq__(self,other):
      """Equal (sufficiently close) to each other"""
      return np.isclose(( self-other).magnitude(),0.0) \
          or np.isclose((-self-other).magnitude(),0.0)

    def __ne__(self,other):
      """Not equal (sufficiently close) to each other"""
      return not self.__eq__(other)


    def normalize(self):
      d = self.magnitude()
      if d > 0.0:
          self.q /= d
          self.p /= d
      return self
      
    def normalized(self):
      return self.copy().normalize()


    def conjugate(self):
      self.p = -self.p
      return self
      
    def conjugated(self):
      return self.copy().conjugate()


    def homomorph(self):
      if self.q < 0.0:
        self.q = -self.q
        self.p = -self.p
      return self
      
    def homomorphed(self):
      return self.copy().homomorph()
      


####################################################################################################
class Rotation:
    u"""
    Orientation stored as unit quaternion.
    
    All methods and naming conventions based on Rowenhorst_etal2015
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
      """Initializes to identity unless specified"""
      self.quaternion = Quaternion2(q=quaternion[0],p=quaternion[1:4])
      self.quaternion.homomorph() # ToDo: Needed?

    def __repr__(self):
      """Value in selected representation"""
      return '\n'.join([
            'Quaternion: {}'.format(self.quaternion),
            'Matrix:\n{}'.format( '\n'.join(['\t'.join(list(map(str,self.asMatrix()[i,:]))) for i in range(3)]) ),
            'Bunge Eulers / deg: {}'.format('\t'.join(list(map(str,self.asEulers(degrees=True)))) ),
              ])


    def asQuaternion(self):
      return self.quaternion.asArray()
      
    def asEulers(self,
                 degrees = False):
      return np.degrees(qu2eu(self.quaternion.asArray())) if degrees else qu2eu(self.quaternion.asArray())
    
    def asAngleAxis(self,
                    degrees = False):
                      
      ax = qu2ax(self.quaternion.asArray())
      if degrees: ax[0] = np.degrees(ax[0])
      
      return ax
      
    def asMatrix(self):
      return qu2om(self.quaternion.asArray())

    def asRodrigues(self):
      return qu2ro(self.quaternion.asArray())
  
    def asHomochoric(self):
     return qu2ho(self.quaternion.asArray())
         
    def asCubochoric(self):
      return qu2cu(self.quaternion.asArray())
      
      
  
    @classmethod
    def fromQuaternion(cls,
                      quaternion,
                      P = -1):

      qu = quaternion
      if P > 0: qu[1:3] *= -1                                                                       # convert from P=1 to P=-1
      if qu[0] < 0.0:
        raise ValueError('Quaternion has negative first component.\n{}'.format(qu[0]))
      if not np.isclose(np.linalg.norm(qu), 1.0):
        raise ValueError('Quaternion is not of unit length.\n{} {} {} {}'.format(*qu))

      return cls(qu)
      
    @classmethod
    def fromEulers(cls,
                   eulers,
                   degrees = False):
                     
      eu = np.radians(eulers) if degrees else eulers
      if np.any(eu < 0.0) or np.any(eu > 2.0*np.pi) or eu[1] > np.pi:
        raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π].\n{} {} {}.'.format(*eu))
        
      return cls(eu2qu(eu))
      
    @classmethod
    def fromAngleAxis(cls,
                      angleAxis,
                      degrees = False,
                      P = -1):
                        
      ax = angleAxis
      if P > 0: ax[1:3] *= -1                                                                       # convert from P=1 to P=-1
      if degrees: ax[0] = np.degrees(ax[0])
      if ax[0] < 0.0 or ax[0] > np.pi:
        raise ValueError('Axis angle rotation angle outside of [0..π].\n'.format(ax[0]))
      if not np.isclose(np.linalg.norm(ax[1:3]), 1.0):
        raise ValueError('Axis angle rotation axis is not of unit length.\n{} {} {}'.format(*ax[1:3]))
      
      return cls(ax2qu(ax))
      
    @classmethod
    def fromMatrix(cls,
                   matrix):
      
      om = matrix 
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
                      P = -1):
        
      ro = rodrigues
      if P > 0: ro[1:3] *= -1                                                                       # convert from P=1 to P=-1
      if not np.isclose(np.linalg.norm(ro[1:3]), 1.0):
        raise ValueError('Rodrigues rotation axis is not of unit length.\n{} {} {}'.format(*ro[1:3]))
      if ro[0] < 0.0:
        raise ValueError('Rodriques rotation angle not positive.\n'.format(ro[0]))
        
      return cls(ro2qu(ro))


    def __mul__(self, other):
      """
      Multiplication
      
      Rotation: Details needed (active/passive), more cases (3,3), (3,3,3,3) need to be considered
      """
      if isinstance(other, Rotation):
        return self.__class__((self.quaternion * other.quaternion).asArray())
      elif isinstance(other, np.ndarray):
        if other.shape == (3,):
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
        elif other.shape == (3,3,):
          raise NotImplementedError
        elif other.shape == (3,3,3,3):
          raise NotImplementedError
        else:
          return NotImplemented
      else:
        return NotImplemented

# ******************************************************************************************
class Quaternion:
    u"""
    Orientation represented as unit quaternion.
    
    All methods and naming conventions based on Rowenhorst_etal2015
    Convention 1: coordinate frames are right-handed
    Convention 2: a rotation angle ω is taken to be positive for a counterclockwise rotation
                  when viewing from the end point of the rotation axis towards the origin
    Convention 3: rotations will be interpreted in the passive sense
    Convention 4: Euler angle triplets are implemented using the Bunge convention,
                  with the angular ranges as [0, 2π],[0, π],[0, 2π]
    Convention 5: the rotation angle ω is limited to the interval [0, π]
    Convention 6: P = -1 (as default)
    
    w is the real part, (x, y, z) are the imaginary parts.

    Vector "a" (defined in coordinate system "A") is passively rotated
               resulting in new coordinates "b" when expressed in system "B".
    b = Q * a
    b = np.dot(Q.asMatrix(),a)
    """

    def __init__(self,
                 quat = None,
                 q    =  1.0,
                 p    = np.zeros(3,dtype=float)):
      """Initializes to identity unless specified"""
      self.q = quat[0] if quat is not None else q
      self.p = np.array(quat[1:4]) if quat is not None else p
      self.homomorph()

    def __iter__(self):
      """Components"""
      return iter(self.asList())

    def __copy__(self):
      """Copy"""
      return self.__class__(q=self.q,p=self.p.copy())

    copy = __copy__

    def __repr__(self):
      """Readable string"""
      return 'Quaternion(real={q:+.6f}, imag=<{p[0]:+.6f}, {p[1]:+.6f}, {p[2]:+.6f}>)'.format(q=self.q,p=self.p)

    def __pow__(self, exponent):
      """Power"""
      omega = math.acos(self.q)
      return self.__class__(q=         math.cos(exponent*omega),
                            p=self.p * math.sin(exponent*omega)/math.sin(omega))

    def __ipow__(self, exponent):
      """In-place power"""
      omega = math.acos(self.q)
      self.q  = math.cos(exponent*omega)
      self.p *= math.sin(exponent*omega)/math.sin(omega)
      return self

    def __mul__(self, other):
      """Multiplication"""
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      try:                                                           # quaternion
          return self.__class__(q=self.q*other.q - np.dot(self.p,other.p),
                                p=self.q*other.p + other.q*self.p + P * np.cross(self.p,other.p))
      except: pass
      try:                                                           # vector (perform passive rotation)
          ( x, y, z) = self.p
          (Vx,Vy,Vz) = other[0:3]
          A = self.q*self.q - np.dot(self.p,self.p)
          B = 2.0 * (x*Vx + y*Vy + z*Vz)
          C = 2.0 * P*self.q

          return np.array([
            A*Vx + B*x + C*(y*Vz - z*Vy),
            A*Vy + B*y + C*(z*Vx - x*Vz),
            A*Vz + B*z + C*(x*Vy - y*Vx),
            ])
      except: pass
      try:                                                           # scalar
          return self.__class__(q=self.q*other,
                                p=self.p*other)
      except:
          return self.copy()

    def __imul__(self, other):
      """In-place multiplication"""
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      try:                                                           # Quaternion
          self.q = self.q*other.q - np.dot(self.p,other.p)
          self.p = self.q*other.p + other.q*self.p + P * np.cross(self.p,other.p)
      except: pass
      return self

    def __div__(self, other):
      """Division"""
      if isinstance(other, (int,float)):
        return self.__class__(q=self.q / other,
                              p=self.p / other)
      else:
          return NotImplemented

    def __idiv__(self, other):
      """In-place division"""
      if isinstance(other, (int,float)):
          self.q /= other
          self.p /= other
      return self

    def __add__(self, other):
      """Addition"""
      if isinstance(other, Quaternion):
        return self.__class__(q=self.q + other.q,
                              p=self.p + other.p)
      else:
          return NotImplemented

    def __iadd__(self, other):
      """In-place addition"""
      if isinstance(other, Quaternion):
          self.q += other.q
          self.p += other.p
      return self

    def __sub__(self, other):
      """Subtraction"""
      if isinstance(other, Quaternion):
          return self.__class__(q=self.q - other.q,
                                p=self.p - other.p)
      else:
          return NotImplemented

    def __isub__(self, other):
      """In-place subtraction"""
      if isinstance(other, Quaternion):
          self.q -= other.q
          self.p -= other.p
      return self

    def __neg__(self):
      """Additive inverse"""
      self.q = -self.q
      self.p = -self.p
      return self

    def __abs__(self):
      """Norm"""
      return math.sqrt(self.q ** 2 + np.dot(self.p,self.p))

    magnitude = __abs__

    def __eq__(self,other):
      """Equal (sufficiently close) to each other"""
      return np.isclose(( self-other).magnitude(),0.0) \
          or np.isclose((-self-other).magnitude(),0.0)

    def __ne__(self,other):
      """Not equal (sufficiently close) to each other"""
      return not self.__eq__(other)

    def __cmp__(self,other):
      """Linear ordering"""
      return (1 if np.linalg.norm(self.asRodrigues()) > np.linalg.norm(other.asRodrigues()) else 0) \
           - (1 if np.linalg.norm(self.asRodrigues()) < np.linalg.norm(other.asRodrigues()) else 0)

    def magnitude_squared(self):
      return self.q ** 2 + np.dot(self.p,self.p)

    def normalize(self):
      d = self.magnitude()
      if d > 0.0:
          self.q /= d
          self.p /= d
      return self

    def conjugate(self):
      self.p = -self.p
      return self

    def homomorph(self):
      if self.q < 0.0:
        self.q = -self.q
        self.p = -self.p
      return self

    def normalized(self):
      return self.copy().normalize()

    def conjugated(self):
      return self.copy().conjugate()

    def homomorphed(self):
      return self.copy().homomorph()

    def asList(self):
      return [self.q]+list(self.p)

    def asM(self):                                                   # to find Averaging Quaternions (see F. Landis Markley et al.)
      return np.outer(self.asList(),self.asList())
      
    def asMatrix(self):
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      qbarhalf = 0.5*(self.q**2 - np.dot(self.p,self.p))
      return 2.0*np.array(
        [[        qbarhalf + self.p[0]**2          ,
           self.p[0]*self.p[1] -P* self.q*self.p[2],
           self.p[0]*self.p[2] +P* self.q*self.p[1] ],
         [ self.p[0]*self.p[1] +P* self.q*self.p[2],
                  qbarhalf + self.p[1]**2          ,
           self.p[1]*self.p[2] -P* self.q*self.p[0] ],
         [ self.p[0]*self.p[2] -P* self.q*self.p[1],
           self.p[1]*self.p[2] +P* self.q*self.p[0],
                  qbarhalf + self.p[2]**2           ],
        ])

    def asAngleAxis(self,
                    degrees = False,
                    flat = False):

      angle = 2.0*math.acos(self.q)

      if np.isclose(angle,0.0):
        angle = 0.0
        axis  = np.array([0.0,0.0,1.0])
      elif np.isclose(self.q,0.0):
        angle = math.pi
        axis  = self.p
      else:
        axis  = np.sign(self.q)*self.p/np.linalg.norm(self.p)

      angle = np.degrees(angle) if degrees else angle

      return np.hstack((angle,axis)) if flat else (angle,axis)

    def asRodrigues(self):
      return np.inf*np.ones(3) if np.isclose(self.q,0.0) else self.p/self.q

    def asEulers(self,
                 degrees = False):
      """Orientation as Bunge-Euler angles."""
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P   = -1.0
      q03 = self.q**2 + self.p[2]**2
      q12 = self.p[0]**2 + self.p[1]**2
      chi = np.sqrt(q03*q12)
      
      if np.isclose(chi,0.0) and np.isclose(q12,0.0):
        eulers = np.array([math.atan2(-2*P*self.q*self.p[2],self.q**2-self.p[2]**2),0,0])
      elif np.isclose(chi,0.0) and np.isclose(q03,0.0):
        eulers = np.array([math.atan2( 2  *self.p[0]*self.p[1],self.p[0]**2-self.p[1]**2),np.pi,0])
      else:
        eulers = np.array([math.atan2((self.p[0]*self.p[2]-P*self.q*self.p[1])/chi,(-P*self.q*self.p[0]-self.p[1]*self.p[2])/chi),
                           math.atan2(2*chi,q03-q12),
                           math.atan2((P*self.q*self.p[1]+self.p[0]*self.p[2])/chi,( self.p[1]*self.p[2]-P*self.q*self.p[0])/chi),
                          ])

      eulers %= 2.0*math.pi                                          # enforce positive angles
      return np.degrees(eulers) if degrees else eulers


#    # Static constructors
    @classmethod
    def fromIdentity(cls):
      return cls()


    @classmethod
    def fromRandom(cls,randomSeed = None):
      import binascii
      if randomSeed is None:
        randomSeed = int(binascii.hexlify(os.urandom(4)),16)
      np.random.seed(randomSeed)
      r = np.random.random(3)
      A = math.sqrt(max(0.0,r[2]))
      B = math.sqrt(max(0.0,1.0-r[2]))
      w = math.cos(2.0*math.pi*r[0])*A
      x = math.sin(2.0*math.pi*r[1])*B
      y = math.cos(2.0*math.pi*r[1])*B
      z = math.sin(2.0*math.pi*r[0])*A
      return cls(quat=[w,x,y,z])


    @classmethod
    def fromRodrigues(cls, rodrigues):
      if not isinstance(rodrigues, np.ndarray): rodrigues = np.array(rodrigues)
      norm = np.linalg.norm(rodrigues)
      halfangle = math.atan(norm)
      s = math.sin(halfangle)
      c = math.cos(halfangle)
      return cls(q=c,p=s*rodrigues/norm)


    @classmethod
    def fromAngleAxis(cls,
                      angle,
                      axis,
                      degrees = False):
      if not isinstance(axis, np.ndarray): axis = np.array(axis,dtype=float)
      axis = axis.astype(float)/np.linalg.norm(axis)
      angle = np.radians(angle) if degrees else angle
      s = math.sin(0.5 * angle)
      c = math.cos(0.5 * angle)
      return cls(q=c,p=axis*s)


    @classmethod
    def fromEulers(cls,
                   eulers,
                   degrees = False):
      if not isinstance(eulers, np.ndarray): eulers = np.array(eulers,dtype=float)
      eulers = np.radians(eulers) if degrees else eulers

      sigma = 0.5*(eulers[0]+eulers[2])
      delta = 0.5*(eulers[0]-eulers[2])
      c = np.cos(0.5*eulers[1])
      s = np.sin(0.5*eulers[1])

      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      w =       c * np.cos(sigma) 
      x =  -P * s * np.cos(delta) 
      y =  -P * s * np.sin(delta)
      z =  -P * c * np.sin(sigma) 
      return cls(quat=[w,x,y,z])


# Modified Method to calculate Quaternion from Orientation Matrix,
# Source: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

    @classmethod
    def fromMatrix(cls, m):
      if m.shape != (3,3) and np.prod(m.shape) == 9:
        m = m.reshape(3,3)

      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      w =   0.5*math.sqrt(max(0.0,1.0+m[0,0]+m[1,1]+m[2,2]))
      x = P*0.5*math.sqrt(max(0.0,1.0+m[0,0]-m[1,1]-m[2,2]))
      y = P*0.5*math.sqrt(max(0.0,1.0-m[0,0]+m[1,1]-m[2,2]))
      z = P*0.5*math.sqrt(max(0.0,1.0-m[0,0]-m[1,1]+m[2,2]))

      x *= -1 if m[2,1] < m[1,2] else 1
      y *= -1 if m[0,2] < m[2,0] else 1
      z *= -1 if m[1,0] < m[0,1] else 1

      return cls(quat=np.array([w,x,y,z])/math.sqrt(w**2 + x**2 + y**2 + z**2))


    @classmethod
    def new_interpolate(cls, q1, q2, t):
        """
        Interpolation
 
        See http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872_2007014421.pdf
        for (another?) way to interpolate quaternions.
        """
        assert isinstance(q1, Quaternion) and isinstance(q2, Quaternion)
        Q = cls()

        costheta = q1.q*q2.q + np.dot(q1.p,q2.p)
        if costheta < 0.:
            costheta = -costheta
            q1 = q1.conjugated()
        elif costheta > 1.:
            costheta = 1.

        theta = math.acos(costheta)
        if abs(theta) < 0.01:
            Q.q = q2.q
            Q.p = q2.p
            return Q

        sintheta = math.sqrt(1.0 - costheta * costheta)
        if abs(sintheta) < 0.01:
            Q.q = (q1.q + q2.q) * 0.5
            Q.p = (q1.p + q2.p) * 0.5
            return Q

        ratio1 = math.sin((1.0 - t) * theta) / sintheta
        ratio2 = math.sin(       t  * theta) / sintheta

        Q.q = q1.q * ratio1 + q2.q * ratio2
        Q.p = q1.p * ratio1 + q2.p * ratio2
        return Q


# ******************************************************************************************
class Symmetry:

  lattices = [None,'orthorhombic','tetragonal','hexagonal','cubic',]

  def __init__(self, symmetry = None):
    """Lattice with given symmetry, defaults to None"""
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

  def symmetryQuats(self,who = []):
    """List of symmetry operations as quaternions."""
    if self.lattice == 'cubic':
      symQuats =  [
                    [ 1.0,              0.0,              0.0,              0.0              ],
                    [ 0.0,              1.0,              0.0,              0.0              ],
                    [ 0.0,              0.0,              1.0,              0.0              ],
                    [ 0.0,              0.0,              0.0,              1.0              ],
                    [ 0.0,              0.0,              0.5*math.sqrt(2), 0.5*math.sqrt(2) ],
                    [ 0.0,              0.0,              0.5*math.sqrt(2),-0.5*math.sqrt(2) ],
                    [ 0.0,              0.5*math.sqrt(2), 0.0,              0.5*math.sqrt(2) ],
                    [ 0.0,              0.5*math.sqrt(2), 0.0,             -0.5*math.sqrt(2) ],
                    [ 0.0,              0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0              ],
                    [ 0.0,             -0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0              ],
                    [ 0.5,              0.5,              0.5,              0.5              ],
                    [-0.5,              0.5,              0.5,              0.5              ],
                    [-0.5,              0.5,              0.5,             -0.5              ],
                    [-0.5,              0.5,             -0.5,              0.5              ],
                    [-0.5,             -0.5,              0.5,              0.5              ],
                    [-0.5,             -0.5,              0.5,             -0.5              ],
                    [-0.5,             -0.5,             -0.5,              0.5              ],
                    [-0.5,              0.5,             -0.5,             -0.5              ],
                    [-0.5*math.sqrt(2), 0.0,              0.0,              0.5*math.sqrt(2) ],
                    [ 0.5*math.sqrt(2), 0.0,              0.0,              0.5*math.sqrt(2) ],
                    [-0.5*math.sqrt(2), 0.0,              0.5*math.sqrt(2), 0.0              ],
                    [-0.5*math.sqrt(2), 0.0,             -0.5*math.sqrt(2), 0.0              ],
                    [-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0,              0.0              ],
                    [-0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0,              0.0              ],
                  ]
    elif self.lattice == 'hexagonal':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [-0.5*math.sqrt(3), 0.0, 0.0,-0.5 ],
                    [ 0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
                    [ 0.0,0.0,0.0,1.0 ],
                    [-0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
                    [-0.5*math.sqrt(3), 0.0, 0.0, 0.5 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,-0.5*math.sqrt(3), 0.5, 0.0 ],
                    [ 0.0, 0.5,-0.5*math.sqrt(3), 0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,-0.5,-0.5*math.sqrt(3), 0.0 ],
                    [ 0.0, 0.5*math.sqrt(3), 0.5, 0.0 ],
                  ]
    elif self.lattice == 'tetragonal':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,0.0,0.0,1.0 ],
                    [ 0.0, 0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0 ],
                    [ 0.0,-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0 ],
                    [ 0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
                    [-0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
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

    return list(map(Quaternion,
                    np.array(symQuats)[np.atleast_1d(np.array(who)) if who != [] else range(len(symQuats))]))
    
    
  def equivalentQuaternions(self,
                            quaternion,
                            who = []):
    """List of symmetrically equivalent quaternions based on own symmetry."""
    return [q*quaternion for q in self.symmetryQuats(who)]


  def inFZ(self,R):
    """Check whether given Rodrigues vector falls into fundamental zone of own symmetry."""
    if isinstance(R, Quaternion): R = R.asRodrigues()                                               # translate accidentally passed quaternion
# fundamental zone in Rodrigues space is point symmetric around origin
    R = abs(R)                                                                                      
    if self.lattice == 'cubic':
      return     math.sqrt(2.0)-1.0 >= R[0] \
             and math.sqrt(2.0)-1.0 >= R[1] \
             and math.sqrt(2.0)-1.0 >= R[2] \
             and 1.0 >= R[0] + R[1] + R[2]
    elif self.lattice == 'hexagonal':
      return     1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2] \
             and 2.0 >= math.sqrt(3)*R[0] + R[1] \
             and 2.0 >= math.sqrt(3)*R[1] + R[0] \
             and 2.0 >= math.sqrt(3) + R[2]
    elif self.lattice == 'tetragonal':
      return     1.0 >= R[0] and 1.0 >= R[1] \
             and math.sqrt(2.0) >= R[0] + R[1] \
             and math.sqrt(2.0) >= R[2] + 1.0
    elif self.lattice == 'orthorhombic':
      return     1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2]
    else:
      return True


  def inDisorientationSST(self,R):
    """
    Check whether given Rodrigues vector (of misorientation) falls into standard stereographic triangle of own symmetry.

    Determination of disorientations follow the work of A. Heinz and P. Neumann:
    Representation of Orientation and Disorientation Data for Cubic, Hexagonal, Tetragonal and Orthorhombic Crystals
    Acta Cryst. (1991). A47, 780-789
    """
    if isinstance(R, Quaternion): R = R.asRodrigues()       # translate accidentially passed quaternion

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
class Orientation:

  __slots__ = ['quaternion','symmetry']

  def __init__(self,
               quaternion = Quaternion.fromIdentity(),
               Rodrigues  = None,
               angleAxis  = None,
               matrix     = None,
               Eulers     = None,
               random     = False,                                                                  # integer to have a fixed seed or True for real random
               symmetry   = None,
               degrees    = False,
              ):
    if random:                                                                                      # produce random orientation
      if isinstance(random, bool ):
        self.quaternion = Quaternion.fromRandom()
      else:
        self.quaternion = Quaternion.fromRandom(randomSeed=random)
    elif isinstance(Eulers, np.ndarray) and Eulers.shape == (3,):                                   # based on given Euler angles
      self.quaternion = Quaternion.fromEulers(Eulers,degrees=degrees)
    elif isinstance(matrix, np.ndarray) :                                                           # based on given rotation matrix
      self.quaternion = Quaternion.fromMatrix(matrix)
    elif isinstance(angleAxis, np.ndarray) and angleAxis.shape == (4,):                             # based on given angle and rotation axis
      self.quaternion = Quaternion.fromAngleAxis(angleAxis[0],angleAxis[1:4],degrees=degrees)
    elif isinstance(Rodrigues, np.ndarray) and Rodrigues.shape == (3,):                             # based on given Rodrigues vector
      self.quaternion = Quaternion.fromRodrigues(Rodrigues)
    elif isinstance(quaternion, Quaternion):                                                        # based on given quaternion
      self.quaternion = quaternion.homomorphed()
    elif (isinstance(quaternion, np.ndarray) and quaternion.shape == (4,)) or \
         (isinstance(quaternion, list)       and len(quaternion)  ==  4  ):                         # based on given quaternion-like array
      self.quaternion = Quaternion(quat=quaternion).homomorphed()

    self.symmetry = Symmetry(symmetry)

  def __copy__(self):
    """Copy"""
    return self.__class__(quaternion=self.quaternion,symmetry=self.symmetry.lattice)

  copy = __copy__


  def __repr__(self):
    """Value as all implemented representations"""
    return '\n'.join([
              'Symmetry: {}'.format(self.symmetry),
              'Quaternion: {}'.format(self.quaternion),
              'Matrix:\n{}'.format( '\n'.join(['\t'.join(list(map(str,self.asMatrix()[i,:]))) for i in range(3)]) ),
              'Bunge Eulers / deg: {}'.format('\t'.join(list(map(str,self.asEulers(degrees=True)))) ),
              ])

  def asQuaternion(self):
    return self.quaternion.asList()

  def asEulers(self,
               degrees = False,
              ):
    return self.quaternion.asEulers(degrees)
  eulers = property(asEulers)

  def asRodrigues(self):
    return self.quaternion.asRodrigues()
  rodrigues = property(asRodrigues)

  def asAngleAxis(self,
                  degrees = False,
                  flat = False):
    return self.quaternion.asAngleAxis(degrees,flat)
  angleAxis = property(asAngleAxis)

  def asMatrix(self):
    return self.quaternion.asMatrix()
  matrix = property(asMatrix)

  def inFZ(self):
    return self.symmetry.inFZ(self.quaternion.asRodrigues())
  infz = property(inFZ)

  def equivalentQuaternions(self,
                            who = []):
    return self.symmetry.equivalentQuaternions(self.quaternion,who)

  def equivalentOrientations(self,
                             who = []):
    return [Orientation(quaternion = q, symmetry = self.symmetry.lattice) for q in self.equivalentQuaternions(who)]

  def reduced(self):
    """Transform orientation to fall into fundamental zone according to symmetry"""
    for me in self.symmetry.equivalentQuaternions(self.quaternion):
      if self.symmetry.inFZ(me.asRodrigues()): break

    return Orientation(quaternion=me,symmetry=self.symmetry.lattice)


  def disorientation(self,
                     other,
                     SST = True):
    """
    Disorientation between myself and given other orientation.

    Rotation axis falls into SST if SST == True.
    (Currently requires same symmetry for both orientations.
     Look into A. Heinz and P. Neumann 1991 for cases with differing sym.)
    """
    if self.symmetry != other.symmetry: raise TypeError('disorientation between different symmetry classes not supported yet.')

    misQ = other.quaternion*self.quaternion.conjugated()
    mySymQs    =  self.symmetry.symmetryQuats() if SST else self.symmetry.symmetryQuats()[:1]       # take all or only first sym operation
    otherSymQs = other.symmetry.symmetryQuats()
    
    for i,sA in enumerate(mySymQs):
      for j,sB in enumerate(otherSymQs):
        theQ = sB*misQ*sA.conjugated()
        for k in range(2):
          theQ.conjugate()
          breaker = self.symmetry.inFZ(theQ) \
                    and (not SST or other.symmetry.inDisorientationSST(theQ))
          if breaker: break
        if breaker: break
      if breaker: break

# disorientation, own sym, other sym, self-->other: True, self<--other: False
    return (Orientation(quaternion = theQ,symmetry = self.symmetry.lattice),
            i,j, k == 1)                                                                             


  def inversePole(self,
                  axis,
                  proper = False,
                  SST = True):
    """Axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)"""
    if SST:                                                                                         # pole requested to be within SST
      for i,q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):                   # test all symmetric equivalent quaternions
        pole = q*axis                                                                               # align crystal direction to axis
        if self.symmetry.inSST(pole,proper): break                                                  # found SST version
    else:
      pole = self.quaternion*axis                                                                   # align crystal direction to axis

    return (pole,i if SST else 0)

  def IPFcolor(self,axis):
    """TSL color of inverse pole figure for given axis"""
    color = np.zeros(3,'d')

    for q in self.symmetry.equivalentQuaternions(self.quaternion):
      pole = q*axis                                                                                 # align crystal direction to axis
      inSST,color = self.symmetry.inSST(pole,color=True)
      if inSST: break

    return color

  @classmethod
  def average(cls,
              orientations,
              multiplicity = []):
    """
    Average orientation

    ref: F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
         Averaging Quaternions,
         Journal of Guidance, Control, and Dynamics, Vol. 30, No. 4 (2007), pp. 1193-1197.
         doi: 10.2514/1.28949
    usage:
         a = Orientation(Eulers=np.radians([10, 10, 0]), symmetry='hexagonal')
         b = Orientation(Eulers=np.radians([20, 0, 0]),  symmetry='hexagonal')
         avg = Orientation.average([a,b])
    """
    if not all(isinstance(item, Orientation) for item in orientations):
      raise TypeError("Only instances of Orientation can be averaged.")

    N = len(orientations)
    if multiplicity == [] or not multiplicity:
      multiplicity = np.ones(N,dtype='i')

    reference = orientations[0]                                                                     # take first as reference
    for i,(o,n) in enumerate(zip(orientations,multiplicity)):
      closest = o.equivalentOrientations(reference.disorientation(o,SST = False)[2])[0]             # select sym orientation with lowest misorientation
      M = closest.quaternion.asM() * n if i == 0 else M + closest.quaternion.asM() * n              # noqa add (multiples) of this orientation to average noqa
    eig, vec = np.linalg.eig(M/N)

    return Orientation(quaternion = Quaternion(quat = np.real(vec.T[eig.argmax()])),
                       symmetry = reference.symmetry.lattice)


  def related(self,
              relationModel,
              direction,
              targetSymmetry = 'cubic'):
    """
    Orientation relationship

    positive number: fcc --> bcc
    negative number: bcc --> fcc
    """
    if relationModel not in ['KS','GT','GTdash','NW','Pitsch','Bain']:  return None
    if int(direction) == 0:  return None

    # KS from S. Morito et al./Journal of Alloys and Compounds 5775 (2013) S587-S592
    # for KS rotation matrices also check K. Kitahara et al./Acta Materialia 54 (2006) 1279-1288
    # GT from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
    # GT' from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
    # NW from H. Kitahara et al./Materials Characterization 54 (2005) 378-386
    # Pitsch from Y. He et al./Acta Materialia 53 (2005) 1179-1190
    # Bain from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81

    variant  = int(abs(direction))-1
    (me,other)  = (0,1) if direction > 0 else (1,0)

    planes = {'KS': \
                    np.array([[[  1,  1,  1],[  0,  1,  1]],
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
                              [[  1,  1, -1],[  0,  1,  1]]]),
              'GT': \
                    np.array([[[  1,  1,  1],[  1,  0,  1]],
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
                              [[  1, -1,  1],[  1,  0,  1]]]),
              'GTdash': \
                    np.array([[[  7, 17, 17],[ 12,  5, 17]],
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
                              [[-17, 17, -7],[-17,  5,-12]]]),
              'NW': \
                    np.array([[[  1,  1,  1],[  0,  1,  1]],
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
                              [[ -1, -1,  1],[  0,  1,  1]]]),
              'Pitsch': \
                    np.array([[[  0,  1,  0],[ -1,  0,  1]],
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
                              [[  0,  0,  1],[ -1,  1,  0]]]),
              'Bain': \
                    np.array([[[  1,  0,  0],[  1,  0,  0]],
                              [[  0,  1,  0],[  0,  1,  0]],
                              [[  0,  0,  1],[  0,  0,  1]]]),
              }

    normals = {'KS': \
                    np.array([[[ -1,  0,  1],[ -1, -1,  1]],
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
                              [[  1,  0,  1],[ -1,  1, -1]]]),
              'GT': \
                    np.array([[[ -5,-12, 17],[-17, -7, 17]],
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
                              [[-17,-12,  5],[-17,  7, 17]]]),
              'GTdash': \
                    np.array([[[  0,  1, -1],[  1,  1, -1]],
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
                              [[  1,  1,  0],[  1,  1, -1]]]),
              'NW': \
                    np.array([[[  2, -1, -1],[  0, -1,  1]],
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
                              [[ -1, -1, -2],[  0, -1,  1]]]),
              'Pitsch': \
                    np.array([[[  1,  0,  1],[  1, -1,  1]],
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
                              [[  1,  1,  0],[  1,  1, -1]]]),
              'Bain': \
                    np.array([[[  0,  1,  0],[  0,  1,  1]],
                              [[  0,  0,  1],[  1,  0,  1]],
                              [[  1,  0,  0],[  1,  1,  0]]]),
              }
    myPlane   = [float(i) for i in planes[relationModel][variant,me]]                               # map(float, planes[...]) does not work in python 3
    myPlane  /= np.linalg.norm(myPlane)
    myNormal  = [float(i) for i in normals[relationModel][variant,me]]                              # map(float, planes[...]) does not work in python 3
    myNormal /= np.linalg.norm(myNormal)
    myMatrix  = np.array([myNormal,np.cross(myPlane,myNormal),myPlane]).T

    otherPlane   = [float(i) for i in planes[relationModel][variant,other]]                         # map(float, planes[...]) does not work in python 3
    otherPlane  /= np.linalg.norm(otherPlane)
    otherNormal  = [float(i) for i in normals[relationModel][variant,other]]                        # map(float, planes[...]) does not work in python 3
    otherNormal /= np.linalg.norm(otherNormal)
    otherMatrix  = np.array([otherNormal,np.cross(otherPlane,otherNormal),otherPlane]).T

    rot=np.dot(otherMatrix,myMatrix.T)

    return Orientation(matrix=np.dot(rot,self.asMatrix()),symmetry=targetSymmetry)

####################################################################################################
# Code below available according to the followin conditions on https://github.com/MarDiehl/3Drotations
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
  return np.isclose(a,1.0,atol=1.0e-15,rtol=0.0)
  
def iszero(a):
  return np.isclose(a,0.0,atol=1.0e-300,rtol=0.0)
  
  
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


def eu2qu(eu):
  """Euler angles to quaternion"""
  ee = 0.5*eu
  cPhi = np.cos(ee[1])
  sPhi = np.sin(ee[1])
  qu = np.array([     cPhi*np.cos(ee[0]+ee[2]),
                   -P*sPhi*np.cos(ee[0]-ee[2]),
                   -P*sPhi*np.sin(ee[0]-ee[2]),
                   -P*cPhi*np.sin(ee[0]+ee[2]) ])
  #if qu[0] < 0.0: qu.homomorph()                                      !ToDo: Check with original
  return qu


def om2eu(om):
  """Euler angles to orientation matrix"""
  if isone(om[2,2]**2):
    eu = np.array([np.arctan2( om[0,1],om[0,0]), np.pi*0.5*(1-om[2,2]),0.0])                        # following the paper, not the reference implementation
  else:
    zeta = 1.0/np.sqrt(1.0-om[2,2]**2)
    eu = np.array([np.arctan2(om[2,0]*zeta,-om[2,1]*zeta),
                   np.arccos(om[2,2]),
                   np.arctan2(om[0,2]*zeta, om[1,2]*zeta)])
  
  # reduce Euler angles to definition range, i.e a lower limit of 0.0
  eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
  return eu


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


def qu2eu(qu):
  """Quaternion to Euler angles""" 
  q03 = qu[0]**2+qu[3]**2
  q12 = qu[1]**2+qu[2]**2
  chi = np.sqrt(q03*q12)
  
  if iszero(chi):
    eu = np.array([np.arctan2(-P*2.0*qu[0]*qu[3],qu[0]**2-qu[3]**2), 0.0,   0.0]) if iszero(q12) else \
         np.array([np.arctan2(2.0*qu[1]*qu[2],qu[1]**2-qu[2]**2),         np.pi, 0.0])
  else:
    #chiInv = 1.0/chi         ToDo: needed for what?
    eu = np.array([np.arctan2((-P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-P*qu[0]*qu[1]-qu[2]*qu[3])*chi ),
                   np.arctan2( 2.0*chi, q03-q12 ), 
                   np.arctan2(( P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-P*qu[0]*qu[1]+qu[2]*qu[3])*chi )])
  
  # reduce Euler angles to definition range, i.e a lower limit of 0.0
  eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
  return eu


def ax2ho(ax):
  """Axis angle to homochoric""" 
  f = (0.75 * ( ax[3] - np.sin(ax[3]) ))**(1.0/3.0)
  ho = ax[0:3] * f
  return ho


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
    ax = np.append(ho/np.sqrt(hmag_squared),2.0*np.arccos(s))                                       # ToDo: Check sanity check in reference implementation 

  return ax


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


def ax2qu(ax):
  """Axis angle to quaternion"""
  if iszero(ax[3]):
    qu = np.array([ 1.0, 0.0, 0.0, 0.0 ])
  else:
    c = np.cos(ax[3]*0.5)
    s = np.sin(ax[3]*0.5)
    qu = np.array([ c, ax[0]*s, ax[1]*s, ax[2]*s ])
  
  return qu


def ro2ho(ro):
  """Rodrigues vector to homochoric""" 
  if iszero(np.sum(ro[0:3]**2.0)):
    ho = [ 0.0, 0.0, 0.0 ]
  else:
    f = 2.0*np.arctan(ro[3]) -np.sin(2.0*np.arctan(ro[3])) if np.isfinite(ro[3]) else np.pi
    ho = ro[0:3] * (0.75*f)**(1.0/3.0)

  return np.array(ho)


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


def om2qu(om):
  """Orientation matrix to quaternion"""
  s = [+om[0,0] +om[1,1] +om[2,2] +1.0, 
       +om[0,0] -om[1,1] -om[2,2] +1.0,
       -om[0,0] +om[1,1] -om[2,2] +1.0, 
       -om[0,0] -om[1,1] +om[2,2] +1.0]
  s = np.maximum(np.zeros(4),s)
  qu = np.sqrt(s)*0.5*np.array([1.0,P,P,P])
  # verify the signs (q0 always positive)
  #ToDo: Here I donot understand the original shortcut from paper to implementation

  qu /= np.linalg.norm(qu)
  if any(isone(abs(qu))): qu[np.where(np.logical_not(isone(qu)))] = 0.0
  if om[2,1] < om[1,2]: qu[1] *= -1.0
  if om[0,2] < om[2,0]: qu[2] *= -1.0
  if om[1,0] < om[0,1]: qu[3] *= -1.0
  if any(om2ax(om)[0:3]*qu[1:4] < 0.0): print(om2ax(om),qu) # something is wrong here
  return qu


def qu2ax(qu):
  """Quaternion to axis angle"""
  omega = 2.0 * np.arccos(qu[0])
  if iszero(omega):                                                                                 #  return axis as [001]  if the angle is zero
    ax = [ 0.0, 0.0, 1.0, 0.0 ]
  elif not iszero(qu[0]):
    s = np.sign(qu[0])/np.sqrt(qu[1]**2+qu[2]**2+qu[3]**2)
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
         [ qu[1]/s,  qu[2]/s,  qu[3]/s, np.tan(np.arccos(qu[0]))]
 
  return np.array(ro)


def qu2ho(qu):
  """Quaternion to homochoric"""
  omega = 2.0 * np.arccos(qu[0])
  
  if iszero(omega):
    ho = np.array([ 0.0, 0.0, 0.0 ])
  else:
    ho = np.array([qu[1], qu[2], qu[3]])
    f  = 0.75 * ( omega - np.sin(omega) )
    ho = ho/np.linalg.norm(ho) * f**(1./3.)

  return ho


def ho2cu(ho):
  """Homochoric to cubochoric"""
  return  Lambert.BallToCube(ho)


def cu2ho(cu):
  """Cubochoric to homochoric"""
  return  Lambert.CubeToBall(cu)


def ro2eu(ro):
  """Rodrigues vector to orientation matrix"""
  return om2eu(ro2om(ro))


def eu2ho(eu):
  """Euler angles to homochoric"""
  return ax2ho(eu2ax(eu))


def om2ro(om):
  """Orientation matrix to Rodriques vector"""
  return eu2ro(om2eu(om))


def om2ho(om):
  """Orientation matrix to homochoric"""
  return ax2ho(om2ax(om))


def ax2eu(ax):
  """Orientation matrix to Euler angles"""
  return om2eu(ax2om(ax))


def ro2om(ro):
 """Rodgrigues vector to orientation matrix"""
 return ax2om(ro2ax(ro))


def ro2qu(ro):
  """Rodrigues vector to quaternion"""
  return ax2qu(ro2ax(ro))


def ho2eu(ho):
  """Homochoric to Euler angles"""
  return ax2eu(ho2ax(ho))


def ho2om(ho):
  """Homochoric to orientation matrix"""
  return ax2om(ho2ax(ho))


def ho2ro(ho):
  """Axis angle to Rodriques vector"""
  return ax2ro(ho2ax(ho))


def ho2qu(ho):
  """Homochoric to quaternion"""
  return ax2qu(ho2ax(ho))


def eu2cu(eu):
  """Euler angles to cubochoric"""
  return ho2cu(eu2ho(eu))


def om2cu(om):
  """Orientation matrix to cubochoric"""
  return ho2cu(om2ho(om))


def ax2cu(ax):
  """Axis angle to cubochoric"""
  return ho2cu(ax2ho(ax))


def ro2cu(ro):
  """Rodrigues vector to cubochoric"""
  return ho2cu(ro2ho(ro))


def qu2cu(qu):
  """Quaternion to cubochoric"""
  return ho2cu(qu2ho(qu))


def cu2eu(cu):
  """Cubochoric to Euler angles"""
  return ho2eu(cu2ho(cu))


def cu2om(cu):
  """Cubochoric to orientation matrix"""
  return ho2om(cu2ho(cu))


def cu2ax(cu):
  """Cubochoric to axis angle"""
  return ho2ax(cu2ho(cu))


def cu2ro(cu):
  """Cubochoric to Rodrigues vector"""
  return ho2ro(cu2ho(cu))


def cu2qu(cu):
  """Cubochoric to quaternion"""
  return ho2qu(cu2ho(cu))
