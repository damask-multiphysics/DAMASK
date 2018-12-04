# -*- coding: UTF-8 no BOM -*-

###################################################
# NOTE: everything here needs to be a np array #
###################################################

import math,os
import numpy as np

# ******************************************************************************************
class Rodrigues:

    def __init__(self, vector = np.zeros(3)):
      self.vector = vector

    def asQuaternion(self):
      norm = np.linalg.norm(self.vector)
      halfAngle = np.arctan(norm)
      return Quaternion(np.cos(halfAngle),np.sin(halfAngle)*self.vector/norm)

    def asAngleAxis(self):
      norm = np.linalg.norm(self.vector)
      halfAngle = np.arctan(norm)
      return (2.0*halfAngle,self.vector/norm)



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

    w is the real part, (x, y, z) are the imaginary parts.

    Vector "a" (defined in coordinate system "A") is passively rotated
               resulting in new coordinates "b" when expressed in system "B".
    b = Q * a
    b = np.dot(Q.asMatrix(),a)
    """

    def __init__(self,
                 quatArray = [1.0,0.0,0.0,0.0]):
      """Initializes to identity unless specified"""
      (self.w,
       self.x,
       self.y,
       self.z ) = quatArray
      self.homomorph()

    def __iter__(self):
      """Components"""
      return iter([self.w,self.x,self.y,self.z])

    def __copy__(self):
      """Copy"""
      Q = Quaternion([self.w,self.x,self.y,self.z])
      return Q

    copy = __copy__

    def __repr__(self):
      """Readable string"""
      return 'Quaternion(real=%+.6f, imag=<%+.6f, %+.6f, %+.6f>)' % (self.w, self.x, self.y, self.z)

    def __pow__(self, exponent):
      """Power"""
      omega = math.acos(self.w)
      vRescale = math.sin(exponent*omega)/math.sin(omega)
      Q = Quaternion()
      Q.w = math.cos(exponent*omega)
      Q.x = self.x * vRescale
      Q.y = self.y * vRescale
      Q.z = self.z * vRescale
      return Q

    def __ipow__(self, exponent):
      """In-place power"""
      omega    = math.acos(self.w)
      vRescale = math.sin(exponent*omega)/math.sin(omega)
      self.w  = np.cos(exponent*omega)
      self.x *= vRescale
      self.y *= vRescale
      self.z *= vRescale
      return self

    def __mul__(self, other):
      """Multiplication"""
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      try:                                                          # quaternion
          Aw = self.w
          Ax = self.x
          Ay = self.y
          Az = self.z
          Bw = other.w
          Bx = other.x
          By = other.y
          Bz = other.z
          Q = Quaternion()
          Q.w = - Ax * Bx - Ay * By - Az * Bz + Aw * Bw
          Q.x = + Ax * Bw + Aw * Bx + P * (Ay * Bz - Az * By)
          Q.y = + Ay * Bw + Aw * By + P * (Az * Bx - Ax * Bz)
          Q.z = + Az * Bw + Aw * Bz + P * (Ax * By - Ay * Bx)
          return Q
      except: pass
      try:                                                         # vector (perform active rotation, i.e. q*v*q.conjugated)
          w = self.w
          x = self.x
          y = self.y
          z = self.z
          Vx = other[0]
          Vy = other[1]
          Vz = other[2]

          A = w**2 - x**2 - y**2 - z**2
          B = 2.0*(x*Vx +  y*Vy + z*Vz)
                    
          return np.array([
            A*Vx + B*x + 2*P*w * (y*Vz - z*Vy),
            A*Vy + B*y + 2*P*w * (z*Vx - x*Vz),
            A*Vz + B*z + 2*P*w * (x*Vy - y*Vx),
            ])
      except: pass
      try:                                                        # scalar
          Q = self.copy()
          Q.w *= other
          Q.x *= other
          Q.y *= other
          Q.z *= other
          return Q
      except:
          return self.copy()

    def __imul__(self, other):
      """In-place multiplication"""
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      try:                                                        # Quaternion
          Aw = self.w
          Ax = self.x
          Ay = self.y
          Az = self.z
          Bw = other.w
          Bx = other.x
          By = other.y
          Bz = other.z
          self.w = - Ax * Bx - Ay * By - Az * Bz + Aw * Bw
          self.x = + Ax * Bw + Aw * Bx + P * (Ay * Bz - Az * By)
          self.y = + Ay * Bw + Aw * By + P * (Az * Bx - Ax * Bz)
          self.z = + Az * Bw + Aw * Bz + P * (Ax * By - Ay * Bx)
      except: pass
      return self

    def __div__(self, other):
      """Division"""
      if isinstance(other, (int,float)):
        w = self.w / other
        x = self.x / other
        y = self.y / other
        z = self.z / other
        return self.__class__([w,x,y,z])
      else:
          return NotImplemented

    def __idiv__(self, other):
      """In-place division"""
      if isinstance(other, (int,float)):
          self.w /= other
          self.x /= other
          self.y /= other
          self.z /= other
      return self

    def __add__(self, other):
      """Addition"""
      if isinstance(other, Quaternion):
        w = self.w + other.w
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return self.__class__([w,x,y,z])
      else:
          return NotImplemented

    def __iadd__(self, other):
      """In-place addition"""
      if isinstance(other, Quaternion):
          self.w += other.w
          self.x += other.x
          self.y += other.y
          self.z += other.z
      return self

    def __sub__(self, other):
      """Subtraction"""
      if isinstance(other, Quaternion):
          Q = self.copy()
          Q.w -= other.w
          Q.x -= other.x
          Q.y -= other.y
          Q.z -= other.z
          return Q
      else:
          return self.copy()

    def __isub__(self, other):
      """In-place subtraction"""
      if isinstance(other, Quaternion):
          self.w -= other.w
          self.x -= other.x
          self.y -= other.y
          self.z -= other.z
      return self

    def __neg__(self):
      """Additive inverse"""
      self.w = -self.w
      self.x = -self.x
      self.y = -self.y
      self.z = -self.z
      return self

    def __abs__(self):
      """Norm"""
      return math.sqrt(self.w ** 2 + \
                       self.x ** 2 + \
                       self.y ** 2 + \
                       self.z ** 2)

    magnitude = __abs__

    def __eq__(self,other):
      """Equal at e-8 precision"""
      return (abs(self.w-other.w) < 1e-8 and \
              abs(self.x-other.x) < 1e-8 and \
              abs(self.y-other.y) < 1e-8 and \
              abs(self.z-other.z) < 1e-8) \
              or \
             (abs(-self.w-other.w) < 1e-8 and \
              abs(-self.x-other.x) < 1e-8 and \
              abs(-self.y-other.y) < 1e-8 and \
              abs(-self.z-other.z) < 1e-8)

    def __ne__(self,other):
      """Not equal at e-8 precision"""
      return not self.__eq__(self,other)

    def __cmp__(self,other):
      """Linear ordering"""
      return (self.Rodrigues()>other.Rodrigues()) - (self.Rodrigues()<other.Rodrigues())

    def magnitude_squared(self):
      return self.w ** 2 + \
             self.x ** 2 + \
             self.y ** 2 + \
             self.z ** 2

    def identity(self):
      self.w = 1.
      self.x = 0.
      self.y = 0.
      self.z = 0.
      return self

    def normalize(self):
      d = self.magnitude()
      if d > 0.0:
          self /= d
      return self

    def conjugate(self):
      self.x = -self.x
      self.y = -self.y
      self.z = -self.z
      return self

    def inverse(self):
      d = self.magnitude()
      if d > 0.0:
        self.conjugate()
        self /= d
      return self

    def homomorph(self):
      if self.w < 0.0:
        self.w = -self.w
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z
      return self

    def normalized(self):
      return self.copy().normalize()

    def conjugated(self):
      return self.copy().conjugate()

    def inversed(self):
      return self.copy().inverse()

    def homomorphed(self):
      return self.copy().homomorph()

    def asList(self):
      return [i for i in self]

    def asM(self):                                                # to find Averaging Quaternions (see F. Landis Markley et al.)
      return np.outer([i for i in self],[i for i in self])
      
    def asMatrix(self):
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      qbarhalf = 0.5*(self.w**2 - self.x**2 - self.y**2 - self.z**2)
      return 2.0*np.array(
        [[      qbarhalf + self.x**2      , self.x*self.y -P* self.w*self.z, self.x*self.z +P* self.w*self.y],
         [ self.x*self.y +P* self.w*self.z,      qbarhalf + self.y**2      , self.y*self.z -P* self.w*self.x],
         [ self.x*self.z -P* self.w*self.y, self.y*self.z +P* self.w*self.x,      qbarhalf + self.z**2      ],
        ])

    def asAngleAxis(self,
                    degrees = False):
      if self.w > 1:
          self.normalize()

      s = math.sqrt(1. - self.w**2)
      x = 2*self.w**2 - 1.
      y = 2*self.w * s

      angle = math.atan2(y,x)
      if angle < 0.0:
        angle *= -1.
        s     *= -1.

      return (np.degrees(angle) if degrees else angle,
              np.array([1.0, 0.0, 0.0] if np.abs(angle) < 1e-6 else [self.x / s, self.y / s, self.z / s]))

    def asRodrigues(self):
      return np.inf*np.ones(3) if self.w == 0.0 else np.array([self.x, self.y, self.z])/self.w

    def asEulers(self,
                 degrees = False):
      """Orientation as Bunge-Euler angles."""
      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P   = -1.0
      q03 = self.w**2 + self.z**2
      q12 = self.x**2 + self.y**2
      chi = np.sqrt(q03*q12)
      
      if abs(chi) < 1e-10 and abs(q12) < 1e-10:
        eulers = np.array([math.atan2(-2*P*self.w*self.z,self.w**2-self.z**2),0,0])
      elif abs(chi) < 1e-10 and abs(q03) < 1e-10:
        eulers = np.array([math.atan2( 2  *self.x*self.y,self.x**2-self.y**2),np.pi,0])
      else:
        eulers = np.array([math.atan2((self.x*self.z-P*self.w*self.y)/chi,(-P*self.w*self.x-self.y*self.z)/chi),
                           math.atan2(2*chi,q03-q12),
                           math.atan2((P*self.w*self.y+self.x*self.z)/chi,( self.y*self.z-P*self.w*self.x)/chi),
                          ])

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
      w = math.cos(2.0*math.pi*r[0])*math.sqrt(r[2])
      x = math.sin(2.0*math.pi*r[1])*math.sqrt(1.0-r[2])
      y = math.cos(2.0*math.pi*r[1])*math.sqrt(1.0-r[2])
      z = math.sin(2.0*math.pi*r[0])*math.sqrt(r[2])
      return cls([w,x,y,z])


    @classmethod
    def fromRodrigues(cls, rodrigues):
      if not isinstance(rodrigues, np.ndarray): rodrigues = np.array(rodrigues)
      halfangle = math.atan(np.linalg.norm(rodrigues))
      c = math.cos(halfangle)
      w = c
      x,y,z = rodrigues/c
      return cls([w,x,y,z])


    @classmethod
    def fromAngleAxis(cls,
                      angle,
                      axis,
                      degrees = False):
      if not isinstance(axis, np.ndarray): axis = np.array(axis,dtype='d')
      axis = axis.astype(float)/np.linalg.norm(axis)
      angle = np.radians(angle) if degrees else angle
      s = math.sin(0.5 * angle)
      w = math.cos(0.5 * angle)
      x = axis[0] * s
      y = axis[1] * s
      z = axis[2] * s
      return cls([w,x,y,z])


    @classmethod
    def fromEulers(cls,
                   eulers,
                   degrees = False):
      if not isinstance(eulers, np.ndarray): eulers = np.array(eulers,dtype='d')
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
      return cls([w,x,y,z])


# Modified Method to calculate Quaternion from Orientation Matrix,
# Source: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

    @classmethod
    def fromMatrix(cls, m):
      if m.shape != (3,3) and np.prod(m.shape) == 9:
        m = m.reshape(3,3)

      # Rowenhorst_etal2015 MSMSE: value of P is selected as -1
      P = -1.0
      w =  0.5*math.sqrt(1.+m[0,0]+m[1,1]+m[2,2])
      x = P*0.5*math.sqrt(1.+m[0,0]-m[1,1]-m[2,2])
      y = P*0.5*math.sqrt(1.-m[0,0]+m[1,1]-m[2,2])
      z = P*0.5*math.sqrt(1.-m[0,0]-m[1,1]+m[2,2])

      x *= -1 if m[2,1] < m[1,2] else 1
      y *= -1 if m[0,2] < m[2,0] else 1
      z *= -1 if m[1,0] < m[0,1] else 1

      return cls(np.array([w,x,y,z])/math.sqrt(w**2 + x**2 + y**2 + z**2))


    @classmethod
    def new_interpolate(cls, q1, q2, t):
        """
        Interpolation
 
        See http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872_2007014421.pdf
        for (another?) way to interpolate quaternions.
        """
        assert isinstance(q1, Quaternion) and isinstance(q2, Quaternion)
        Q = cls()

        costheta = q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z
        if costheta < 0.:
            costheta = -costheta
            q1 = q1.conjugated()
        elif costheta > 1:
            costheta = 1

        theta = math.acos(costheta)
        if abs(theta) < 0.01:
            Q.w = q2.w
            Q.x = q2.x
            Q.y = q2.y
            Q.z = q2.z
            return Q

        sintheta = math.sqrt(1.0 - costheta * costheta)
        if abs(sintheta) < 0.01:
            Q.w = (q1.w + q2.w) * 0.5
            Q.x = (q1.x + q2.x) * 0.5
            Q.y = (q1.y + q2.y) * 0.5
            Q.z = (q1.z + q2.z) * 0.5
            return Q

        ratio1 = math.sin((1 - t) * theta) / sintheta
        ratio2 = math.sin(t * theta) / sintheta

        Q.w = q1.w * ratio1 + q2.w * ratio2
        Q.x = q1.x * ratio1 + q2.x * ratio2
        Q.y = q1.y * ratio1 + q2.y * ratio2
        Q.z = q1.z * ratio1 + q2.z * ratio2
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
    """Readbable string"""
    return '%s' % (self.lattice)


  def __eq__(self, other):
    """Equal"""
    return self.lattice == other.lattice

  def __neq__(self, other):
    """Not equal"""
    return not self.__eq__(other)

  def __cmp__(self,other):
    """Linear ordering"""
    myOrder = Symmetry.lattices.index(self.lattice)
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
    if isinstance(R, Quaternion): R = R.asRodrigues()                                               # translate accidentially passed quaternion
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
    """
#     basis = {'cubic' :        np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,1.]/np.sqrt(2.),                     # direction of green
#                                                       [1.,1.,1.]/np.sqrt(3.)]).transpose()),      # direction of blue
#              'hexagonal' :    np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,0.],                                 # direction of green
#                                                       [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).transpose()),      # direction of blue
#              'tetragonal' :   np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,0.],                                 # direction of green
#                                                       [1.,1.,0.]/np.sqrt(2.)]).transpose()),      # direction of blue
#              'orthorhombic' : np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,0.],                                 # direction of green
#                                                       [0.,1.,0.]]).transpose()),                  # direction of blue
#             }

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

    v = np.array(vector,dtype = float)
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
        rgb = np.minimum(np.ones(3,'d'),rgb)                                                        # limit to maximum intensity
        rgb /= max(rgb)                                                                             # normalize to (HS)V = 1
      else:
        rgb = np.zeros(3,'d')
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
    elif isinstance(quaternion, np.ndarray) and quaternion.shape == (4,):                           # based on given quaternion-like array
      self.quaternion = Quaternion(quaternion).homomorphed()

    self.symmetry = Symmetry(symmetry)

  def __copy__(self):
    """Copy"""
    return self.__class__(quaternion=self.quaternion,symmetry=self.symmetry.lattice)

  copy = __copy__


  def __repr__(self):
    """Value as all implemented representations"""
    return 'Symmetry: %s\n' % (self.symmetry) + \
           'Quaternion: %s\n' % (self.quaternion) + \
           'Matrix:\n%s\n' % ( '\n'.join(['\t'.join(map(str,self.asMatrix()[i,:])) for i in range(3)]) ) + \
           'Bunge Eulers / deg: %s' % ('\t'.join(map(str,self.asEulers(degrees=True))) )

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
                  degrees = False):
    return self.quaternion.asAngleAxis(degrees)
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

    return Orientation(quaternion = Quaternion(quatArray = np.real(vec.T[eig.argmax()])),
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
    # for KS rotation matrices also check  K. Kitahara et al./Acta Materialia 54 (2006) 1279-1288
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
