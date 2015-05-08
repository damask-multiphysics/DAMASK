# -*- coding: UTF-8 no BOM -*-

###################################################
# NOTE: everything here needs to be a np array #
###################################################

import math,random
import numpy as np

# ******************************************************************************************
class Rodrigues:
# ******************************************************************************************

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
# ******************************************************************************************
    # All methods and naming conventions based off
    # http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions

    # w is the real part, (x, y, z) are the imaginary parts

    def __init__(self, quatArray=[1.0,0.0,0.0,0.0]):
      self.w, \
      self.x, \
      self.y, \
      self.z = quatArray
      self = self.homomorph()

    def __iter__(self):
      return iter([self.w,self.x,self.y,self.z])

    def __copy__(self):
      Q = Quaternion([self.w,self.x,self.y,self.z])
      return Q

    copy = __copy__

    def __repr__(self):
      return 'Quaternion(real=%+.4f, imag=<%+.4f, %+.4f, %+.4f>)' % \
          (self.w, self.x, self.y, self.z)

    def __pow__(self, exponent):
      omega = math.acos(self.w)
      vRescale = math.sin(exponent*omega)/math.sin(omega)
      Q = Quaternion()
      Q.x = self.x * vRescale
      Q.y = self.y * vRescale
      Q.z = self.z * vRescale
      Q.w = math.cos(exponent*omega)
      return Q

    def __ipow__(self, exponent):
      omega = math.acos(self.w)
      vRescale = math.sin(exponent*omega)/math.sin(omega)
      self.x *= vRescale
      self.y *= vRescale
      self.z *= vRescale
      self.w = np.cos(exponent*omega)
      return self

    def __mul__(self, other):
      try:                                                          # quaternion
          Ax = self.x
          Ay = self.y
          Az = self.z
          Aw = self.w
          Bx = other.x
          By = other.y
          Bz = other.z
          Bw = other.w
          Q = Quaternion()
          Q.x = + Ax * Bw + Ay * Bz - Az * By + Aw * Bx
          Q.y = - Ax * Bz + Ay * Bw + Az * Bx + Aw * By
          Q.z = + Ax * By - Ay * Bx + Az * Bw + Aw * Bz
          Q.w = - Ax * Bx - Ay * By - Az * Bz + Aw * Bw
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

          return np.array([\
             w * w * Vx + 2 * y * w * Vz - 2 * z * w * Vy + \
             x * x * Vx + 2 * y * x * Vy + 2 * z * x * Vz - \
             z * z * Vx - y * y * Vx,
             2 * x * y * Vx + y * y * Vy + 2 * z * y * Vz + \
             2 * w * z * Vx - z * z * Vy + w * w * Vy - \
             2 * x * w * Vz - x * x * Vy,
             2 * x * z * Vx + 2 * y * z * Vy + \
             z * z * Vz - 2 * w * y * Vx - y * y * Vz + \
             2 * w * x * Vy - x * x * Vz + w * w * Vz ])
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
      try:                                                        # Quaternion
          Ax = self.x
          Ay = self.y
          Az = self.z
          Aw = self.w
          Bx = other.x
          By = other.y
          Bz = other.z
          Bw = other.w
          self.x =  Ax * Bw + Ay * Bz - Az * By + Aw * Bx
          self.y = -Ax * Bz + Ay * Bw + Az * Bx + Aw * By
          self.z =  Ax * By - Ay * Bx + Az * Bw + Aw * Bz
          self.w = -Ax * Bx - Ay * By - Az * Bz + Aw * Bw
      except: pass
      return self

    def __div__(self, other):
      if isinstance(other, (int,float,long)):
        w = self.w / other
        x = self.x / other
        y = self.y / other
        z = self.z / other
        return self.__class__([w,x,y,z])
      else:
          return NotImplemented

    def __idiv__(self, other):
      if isinstance(other, (int,float,long)):
          self.w /= other
          self.x /= other
          self.y /= other
          self.z /= other
      return self

    def __add__(self, other):
      if isinstance(other, Quaternion):
        w = self.w + other.w
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return self.__class__([w,x,y,z])
      else:
          return NotImplemented

    def __iadd__(self, other):
      if isinstance(other, Quaternion):
          self.w += other.w
          self.x += other.x
          self.y += other.y
          self.z += other.z
      return self

    def __sub__(self, other):
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
      if isinstance(other, Quaternion):
          self.w -= other.w
          self.x -= other.x
          self.y -= other.y
          self.z -= other.z
      return self

    def __neg__(self):
      self.w = -self.w
      self.x = -self.x
      self.y = -self.y
      self.z = -self.z
      return self

    def __abs__(self):
      return math.sqrt(self.w ** 2 + \
                       self.x ** 2 + \
                       self.y ** 2 + \
                       self.z ** 2)

    magnitude = __abs__

    def __eq__(self,other):
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
      return not __eq__(self,other)

    def __cmp__(self,other):
      return cmp(self.Rodrigues(),other.Rodrigues())

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

    def rotateBy_angleaxis(self, angle, axis):
      self *= Quaternion.fromAngleAxis(angle, axis)
      return self

    def rotateBy_Eulers(self, heading, attitude, bank):
      self *= Quaternion.fromEulers(eulers, type)
      return self

    def rotateBy_matrix(self, m):
      self *= Quaternion.fromMatrix(m)
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
      return np.array([[1.0-2.0*(self.y*self.y+self.z*self.z),     2.0*(self.x*self.y-self.z*self.w),     2.0*(self.x*self.z+self.y*self.w)],
                       [    2.0*(self.x*self.y+self.z*self.w), 1.0-2.0*(self.x*self.x+self.z*self.z),     2.0*(self.y*self.z-self.x*self.w)],
                       [    2.0*(self.x*self.z-self.y*self.w),     2.0*(self.x*self.w+self.y*self.z), 1.0-2.0*(self.x*self.x+self.y*self.y)]])

    def asAngleAxis(self):
      if self.w > 1:
          self.normalize()

      s = math.sqrt(1. - self.w**2)
      x = 2*self.w**2 - 1.
      y = 2*self.w * s

      angle = math.atan2(y,x)

      if angle < 1e-3:
          return angle, np.array([1.0, 0.0, 0.0])
      else:
          return angle, np.array([self.x / s, self.y / s, self.z / s])

    def asRodrigues(self):
      if self.w != 0.0:
          return np.array([self.x, self.y, self.z])/self.w
      else:
          return np.array([float('inf')]*3)

    def asEulers(self,type='bunge'):
      '''
      conversion taken from:
      Melcher, A.; Unser, A.; Reichhardt, M.; Nestler, B.; Pötschke, M.; Selzer, M.
      Conversion of EBSD data by a quaternion based algorithm to be used for grain structure simulations
      Technische Mechanik 30 (2010) pp 401--413
      '''
      angles = [0.0,0.0,0.0]

      if type.lower() == 'bunge' or type.lower() == 'zxz':
        if   abs(self.x - self.y) < 1e-8:
          x = self.w**2 - self.z**2
          y = 2.*self.w*self.z
          angles[0] = math.atan2(y,x)
        elif abs(self.w - self.z) < 1e-8:
          x = self.x**2 - self.y**2
          y = 2.*self.x*self.y
          angles[0] = math.atan2(y,x)
          angles[1] = math.pi
        else:
          chi = math.sqrt((self.w**2 + self.z**2)*(self.x**2 + self.y**2))

          x = (self.w * self.x - self.y * self.z)/2./chi
          y = (self.w * self.y + self.x * self.z)/2./chi
          angles[0] = math.atan2(y,x)

          x = self.w**2 + self.z**2 - (self.x**2 + self.y**2)
          y = 2.*chi
          angles[1] = math.atan2(y,x)

          x = (self.w * self.x + self.y * self.z)/2./chi
          y = (self.z * self.x - self.y * self.w)/2./chi
          angles[2] = math.atan2(y,x)

#        if angles[0] < 0.0:
#          angles[0] += 2*math.pi
#        if angles[1] < 0.0:
#          angles[1] += math.pi
#          angles[2] *= -1
#        if angles[2] < 0.0:
#          angles[2] += 2*math.pi

      return angles


#    # Static constructors
    @classmethod
    def fromIdentity(cls):
        return cls()


    @classmethod
    def fromRandom(cls):
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        w = math.cos(2.0*math.pi*r1)*math.sqrt(r3)
        x = math.sin(2.0*math.pi*r2)*math.sqrt(1.0-r3)
        y = math.cos(2.0*math.pi*r2)*math.sqrt(1.0-r3)
        z = math.sin(2.0*math.pi*r1)*math.sqrt(r3)
        return cls([w,x,y,z])


    @classmethod
    def fromRodrigues(cls, rodrigues):
        if not isinstance(rodrigues, np.ndarray): rodrigues = np.array(rodrigues)
        halfangle = math.atan(np.linalg.norm(rodrigues))
        c = math.cos(halfangle)
        w = c
        x,y,z = c*rodrigues
        return cls([w,x,y,z])


    @classmethod
    def fromAngleAxis(cls, angle, axis):
        if not isinstance(axis, np.ndarray): axis = np.array(axis)
        axis /= np.linalg.norm(axis)
        s = math.sin(angle / 2.0)
        w = math.cos(angle / 2.0)
        x = axis[0] * s
        y = axis[1] * s
        z = axis[2] * s
        return cls([w,x,y,z])


    @classmethod
    def fromEulers(cls, eulers, type = 'Bunge'):
        c1 = math.cos(eulers[0] / 2.0)
        s1 = math.sin(eulers[0] / 2.0)
        c2 = math.cos(eulers[1] / 2.0)
        s2 = math.sin(eulers[1] / 2.0)
        c3 = math.cos(eulers[2] / 2.0)
        s3 = math.sin(eulers[2] / 2.0)

        if type.lower() == 'bunge' or type.lower() == 'zxz':
          w =   c1 * c2 * c3 - s1 * c2 * s3
          x =   c1 * s2 * c3 + s1 * s2 * s3
          y = - c1 * s2 * s3 + s1 * s2 * c3
          z =   c1 * c2 * s3 + s1 * c2 * c3
        else:
#          print 'unknown Euler convention'
          w = c1 * c2 * c3 - s1 * s2 * s3
          x = s1 * s2 * c3 + c1 * c2 * s3
          y = s1 * c2 * c3 + c1 * s2 * s3
          z = c1 * s2 * c3 - s1 * c2 * s3
        return cls([w,x,y,z])


    @classmethod
    def fromMatrix(cls, m):
      if m[0,0] + m[1,1] + m[2,2] > 0.00000001:
        t = m[0,0] + m[1,1] + m[2,2] + 1.0
        s = 0.5/math.sqrt(t)

        return cls(
          [ s*t,
            (m[1,2] - m[2,1])*s,
            (m[2,0] - m[0,2])*s,
            (m[0,1] - m[1,0])*s,
          ])

      elif m[0,0] > m[1,1] and m[0,0] > m[2,2]:
        t = m[0,0] - m[1,1] - m[2,2] + 1.0
        s = 0.5/math.sqrt(t)

        return cls(
          [ (m[1,2] - m[2,1])*s,
            s*t,
            (m[0,1] + m[1,0])*s,
            (m[2,0] + m[0,2])*s,
          ])

      elif m[1,1] > m[2,2]:
        t = -m[0,0] + m[1,1] - m[2,2] + 1.0
        s = 0.5/math.sqrt(t)

        return cls(
          [ (m[2,0] - m[0,2])*s,
            (m[0,1] + m[1,0])*s,
            s*t,
            (m[1,2] + m[2,1])*s,
          ])

      else:
        t = -m[0,0] - m[1,1] + m[2,2] + 1.0
        s = 0.5/math.sqrt(t)

        return cls(
          [ (m[0,1] - m[1,0])*s,
            (m[2,0] + m[0,2])*s,
            (m[1,2] + m[2,1])*s,
            s*t,
          ])


    @classmethod
    def new_interpolate(cls, q1, q2, t):
# see http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872_2007014421.pdf for (another?) way to interpolate quaternions

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
# ******************************************************************************************

  lattices = [None,'orthorhombic','tetragonal','hexagonal','cubic',]

  def __init__(self, symmetry = None):
    if isinstance(symmetry, basestring) and symmetry.lower() in Symmetry.lattices:
      self.lattice = symmetry.lower()
    else:
      self.lattice = None


  def __copy__(self):
    return self.__class__(self.lattice)

  copy = __copy__


  def __repr__(self):
    return '%s' % (self.lattice)


  def __eq__(self, other):
    return self.lattice == other.lattice


  def __neq__(self, other):
    return not self.__eq__(other)

  def __cmp__(self,other):
    return cmp(Symmetry.lattices.index(self.lattice),Symmetry.lattices.index(other.lattice))

  def equivalentQuaternions(self,quaternion):
    '''
    List of symmetrically equivalent quaternions based on own symmetry.
    '''
    if self.lattice == 'cubic':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,0.0,0.0,1.0 ],
                    [ 0.0, 0.0, 0.5*math.sqrt(2), 0.5*math.sqrt(2) ],
                    [ 0.0, 0.0, 0.5*math.sqrt(2),-0.5*math.sqrt(2) ],
                    [ 0.0, 0.5*math.sqrt(2), 0.0, 0.5*math.sqrt(2) ],
                    [ 0.0, 0.5*math.sqrt(2), 0.0,-0.5*math.sqrt(2) ],
                    [ 0.0, 0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0 ],
                    [ 0.0,-0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0 ],
                    [ 0.5, 0.5, 0.5, 0.5 ],
                    [-0.5, 0.5, 0.5, 0.5 ],
                    [-0.5, 0.5, 0.5,-0.5 ],
                    [-0.5, 0.5,-0.5, 0.5 ],
                    [-0.5,-0.5, 0.5, 0.5 ],
                    [-0.5,-0.5, 0.5,-0.5 ],
                    [-0.5,-0.5,-0.5, 0.5 ],
                    [-0.5, 0.5,-0.5,-0.5 ],
                    [-0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
                    [ 0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
                    [-0.5*math.sqrt(2), 0.0, 0.5*math.sqrt(2), 0.0 ],
                    [-0.5*math.sqrt(2), 0.0,-0.5*math.sqrt(2), 0.0 ],
                    [-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0, 0.0 ],
                    [-0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0, 0.0 ],
                  ]
    elif self.lattice == 'hexagonal':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,0.0,0.0,1.0 ],
                    [-0.5*math.sqrt(3), 0.0, 0.0, 0.5 ],
                    [-0.5*math.sqrt(3), 0.0, 0.0,-0.5 ],
                    [ 0.0, 0.5*math.sqrt(3), 0.5, 0.0 ],
                    [ 0.0,-0.5*math.sqrt(3), 0.5, 0.0 ],
                    [ 0.0, 0.5,-0.5*math.sqrt(3), 0.0 ],
                    [ 0.0,-0.5,-0.5*math.sqrt(3), 0.0 ],
                    [ 0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
                    [-0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
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

    return [quaternion*Quaternion(q) for q in symQuats]


  def inFZ(self,R):
    '''
    Check whether given Rodrigues vector falls into fundamental zone of own symmetry.
    '''
    if isinstance(R, Quaternion): R = R.asRodrigues()       # translate accidentially passed quaternion
    R = abs(R)                                              # fundamental zone in Rodrigues space is point symmetric around origin
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
    '''
    Check whether given Rodrigues vector (of misorientation) falls into standard stereographic triangle of own symmetry.
    Determination of disorientations follow the work of A. Heinz and P. Neumann:
    Representation of Orientation and Disorientation Data for Cubic, Hexagonal, Tetragonal and Orthorhombic Crystals
    Acta Cryst. (1991). A47, 780-789
    '''
    if isinstance(R, Quaternion): R = R.asRodrigues()       # translate accidentially passed quaternion

    epsilon = 0.0

    if self.lattice == 'cubic':
      return R[0] >= R[1]+epsilon                and R[1] >= R[2]+epsilon    and R[2] >= epsilon and self.inFZ(R)

    elif self.lattice == 'hexagonal':
      return R[0] >= math.sqrt(3)*(R[1]+epsilon) and R[1] >= epsilon         and R[2] >= epsilon and self.inFZ(R)

    elif self.lattice == 'tetragonal':
      return R[0] >= R[1]+epsilon                and R[1] >= epsilon         and R[2] >= epsilon and self.inFZ(R)

    elif self.lattice == 'orthorhombic':
      return R[0] >= epsilon                     and R[1] >= epsilon         and R[2] >= epsilon and self.inFZ(R)

    else:
      return True


  def inSST(self,vector,color = False):
    '''
    Check whether given vector falls into standard stereographic triangle of own symmetry.
    Return inverse pole figure color if requested.
    '''
#     basis = {'cubic' :        np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,1.]/np.sqrt(2.),                     # direction of green
#                                                             [1.,1.,1.]/np.sqrt(3.)]).transpose()),      # direction of blue
#              'hexagonal' :    np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,0.],                                    # direction of green
#                                                 [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).transpose()),      # direction of blue
#              'tetragonal' :   np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,0.],                                    # direction of green
#                                                             [1.,1.,0.]/np.sqrt(2.)]).transpose()),      # direction of blue
#              'orthorhombic' : np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,0.],                                    # direction of green
#                                                             [0.,1.,0.]]).transpose()),                     # direction of blue
#             }
    if self.lattice == 'cubic':
      basis = np.array([ [-1.            ,  0.            ,  1. ],
                         [ np.sqrt(2.), -np.sqrt(2.),        0. ],
                         [ 0.            ,  np.sqrt(3.),     0. ] ])
    elif self.lattice == 'hexagonal':
      basis = np.array([ [ 0.            ,  0.            ,  1. ],
                         [ 1.            , -np.sqrt(3.),     0. ],
                         [ 0.            ,  2.            ,  0. ] ])
    elif self.lattice == 'tetragonal':
      basis = np.array([ [ 0.            ,  0.            ,  1. ],
                         [ 1.            , -1.            ,  0. ],
                         [ 0.            ,  np.sqrt(2.),     0. ] ])
    elif self.lattice == 'orthorhombic':
      basis = np.array([ [ 0., 0., 1.],
                         [ 1., 0., 0.],
                         [ 0., 1., 0.] ])
    else:
      basis = None

    if basis == None:
      theComponents = -np.ones(3,'d')
    else:
      theComponents = np.dot(basis,np.array([vector[0],vector[1],abs(vector[2])]))

    inSST = np.all(theComponents >= 0.0)

    if color:                                                                      # have to return color array
      if inSST:
        rgb = np.power(theComponents/np.linalg.norm(theComponents),0.5)      # smoothen color ramps
        rgb = np.minimum(np.ones(3,'d'),rgb)                                 # limit to maximum intensity
        rgb /= max(rgb)                                                            # normalize to (HS)V = 1
      else:
        rgb = np.zeros(3,'d')
      return (inSST,rgb)
    else:
      return inSST

# code derived from http://pyeuclid.googlecode.com/svn/trunk/euclid.py
# suggested reading: http://web.mit.edu/2.998/www/QuaternionReport1.pdf



# ******************************************************************************************
class Orientation:
# ******************************************************************************************

  __slots__ = ['quaternion','symmetry']

  def __init__(self,
               quaternion = Quaternion.fromIdentity(),
               Rodrigues  = None,
               angleAxis  = None,
               matrix     = None,
               Eulers     = None,
               random     = False,
               symmetry   = None,
              ):
    if random:                                                                         # produce random orientation
      self.quaternion = Quaternion.fromRandom()
    elif isinstance(Eulers, np.ndarray) and Eulers.shape == (3,):                      # based on given Euler angles
      self.quaternion = Quaternion.fromEulers(Eulers,'bunge')
    elif isinstance(matrix, np.ndarray) and matrix.shape == (3,3):                     # based on given rotation matrix
      self.quaternion = Quaternion.fromMatrix(matrix)
    elif isinstance(angleAxis, np.ndarray) and angleAxis.shape == (4,):                # based on given angle and rotation axis
      self.quaternion = Quaternion.fromAngleAxis(angleAxis[0],angleAxis[1:4])
    elif isinstance(Rodrigues, np.ndarray) and Rodrigues.shape == (3,):                # based on given Rodrigues vector
      self.quaternion = Quaternion.fromRodrigues(Rodrigues)
    elif isinstance(quaternion, Quaternion):                                           # based on given quaternion
      self.quaternion = quaternion.homomorphed()
    elif isinstance(quaternion, np.ndarray) and quaternion.shape == (4,):              # based on given quaternion
      self.quaternion = Quaternion(quaternion).homomorphed()

    self.symmetry = Symmetry(symmetry)

  def __copy__(self):
    return self.__class__(quaternion=self.quaternion,symmetry=self.symmetry.lattice)

  copy = __copy__


  def __repr__(self):
    return 'Symmetry: %s\n' % (self.symmetry) + \
           'Quaternion: %s\n' % (self.quaternion) + \
           'Matrix:\n%s\n' % ( '\n'.join(['\t'.join(map(str,self.asMatrix()[i,:])) for i in range(3)]) ) + \
           'Bunge Eulers / deg: %s' % ('\t'.join(map(lambda x:str(np.degrees(x)),self.asEulers('bunge'))) )

  def asQuaternion(self):
    return self.quaternion.asList()


  def asEulers(self,type='bunge'):
    return self.quaternion.asEulers(type)


  def asRodrigues(self):
    return self.quaternion.asRodrigues()


  def asMatrix(self):
    return self.quaternion.asMatrix()


  def reduced(self):
    '''
    Transform orientation to fall into fundamental zone according to own (or given) symmetry
    '''

    for me in self.symmetry.equivalentQuaternions(self.quaternion):
      if self.symmetry.inFZ(me.asRodrigues()): break

    return Orientation(quaternion=me,symmetry=self.symmetry.lattice)


  def disorientation(self,other):
    '''
    Disorientation between myself and given other orientation
    (either reduced according to my own symmetry or given one)
    '''

    lowerSymmetry = min(self.symmetry,other.symmetry)
    breaker = False

    for me in self.symmetry.equivalentQuaternions(self.quaternion):
      me.conjugate()
      for they in other.symmetry.equivalentQuaternions(other.quaternion):
#        theQ = me * they
        theQ = they * me
#        if theQ.x < 0.0 or theQ.y < 0.0 or theQ.z < 0.0: theQ.conjugate()                   # speed up scanning since minimum angle is usually found for positive x,y,z
        breaker = lowerSymmetry.inDisorientationSST(theQ.asRodrigues()) #\
#               or lowerSymmetry.inDisorientationSST(theQ.conjugated().asRodrigues())
        if breaker: break
      if breaker: break

    return Orientation(quaternion=theQ,symmetry=self.symmetry.lattice) #, me.conjugated(), they


  def inversePole(self,axis,SST = True):
    '''
    axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)
    '''

    for i,q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):                 # test all symmetric equivalent orientations
      if SST:                                                                                   # pole requested to be within SST
          pole = q.conjugated()*axis                                                            # align crystal direction to axis
          if self.symmetry.inSST(pole): break
      else:
        pole = q.conjugated()*axis                                                              # align crystal direction to axis

    return pole

  def IPFcolor(self,axis):
    '''
    TSL color of inverse pole figure for given axis
    '''

    color = np.zeros(3,'d')

    for i,q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):
      pole = q.conjugated()*axis                               # align crystal direction to axis
      inSST,color = self.symmetry.inSST(pole,color=True)
      if inSST: break

    return color

  @classmethod
  def getAverageOrientation(cls, orientationList):
    """RETURN THE AVERAGE ORIENTATION
    ref: F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
         Averaging Quaternions,
         Journal of Guidance, Control, and Dynamics, Vol. 30, No. 4 (2007), pp. 1193-1197.
         doi: 10.2514/1.28949
    usage:
         a = Orientation(Eulers=np.radians([10, 10, 0]), symmetry='hexagonal')
         b = Orientation(Eulers=np.radians([20, 0, 0]), symmetry='hexagonal')
         avg = Orientation.getAverageOrientation([a,b])"""
    if not all(isinstance(item, Orientation) for item in orientationList):
        raise TypeError("Only instances of Orientation can be averaged.")
    n = len(orientationList)
    tmp_m = orientationList.pop(0).quaternion.asM()
    for tmp_o in orientationList:
      tmp_m += tmp_o.quaternion.asM()
    eig, vec = np.linalg.eig(tmp_m/n)
    return Orientation( quaternion=Quaternion(quatArray=vec.T[eig.argmax()]) )


  def related(self, relationModel, direction, targetSymmetry):
  
    if relationModel not in ['KS','GT',"GT'",'NW','Bain']:  return None

    variant  = int(abs(direction))
    me    = 0 if direction > 0 else 1
    other = 1 if direction > 0 else 0

    planes = {'KS': \
                    np.array([[[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]]]),
              'GT': \
                    np.array([[[  1,  1,  1],[  1,  1,  0]],\
                              [[  1,  1,  1],[  1,  0,  1]],\
                              [[ -1, -1,  1],[ -1, -1,  0]],\
                              [[ -1, -1,  1],[ -1,  0,  1]],\
                              [[ -1,  1,  1],[ -1,  1,  0]],\
                              [[ -1,  1,  1],[ -1,  0,  1]],\
                              [[  1, -1,  1],[  1, -1,  0]],\
                              [[  1, -1,  1],[  1,  0,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  1,  1,  0]],\
                              [[ -1, -1,  1],[  0, -1,  1]],\
                              [[ -1, -1,  1],[ -1, -1,  0]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[ -1,  1,  0]],\
                              [[  1, -1,  1],[  0, -1,  1]],\
                              [[  1, -1,  1],[  1, -1,  0]],\
                              [[  1,  1,  1],[  1,  0,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[ -1, -1,  1],[ -1,  0,  1]],\
                              [[ -1, -1,  1],[  0, -1,  1]],\
                              [[ -1,  1,  1],[ -1,  0,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  1,  0,  1]],\
                              [[  1, -1,  1],[  0, -1,  1]]]),
              "GT'": \
                    np.array([[[ 17,  7, 17],[ 17, 12,  5]],\
                              [[-17,  7,-17],[-17, 12, -5]],\
                              [[-17, -7, 17],[-17,-12,  5]],\
                              [[ 17, -7,-17],[ 17,-12, -5]],\
                              [[ 17, 17,  7],[ 17,  5, 12]],\
                              [[-17,-17,  7],[-17, -5, 12]],\
                              [[ 17,-17, -7],[ 17, -5,-12]],\
                              [[-17, 17, -7],[-17,  5,-12]],\
                              [[ 17, 17,  7],[  5, 17, 12]],\
                              [[-17,-17,  7],[ -5,-17, 12]],\
                              [[ 17,-17, -7],[  5,-17,-12]],\
                              [[-17, 17, -7],[ -5, 17,-12]],\
                              [[  7, 17, 17],[ 12, 17,  5]],\
                              [[  7,-17,-17],[ 12,-17, -5]],\
                              [[ -7,-17, 17],[-12,-17,  5]],\
                              [[ -7, 17,-17],[-12, 17, -5]],\
                              [[  7, 17, 17],[ 12,  5, 17]],\
                              [[  7,-17,-17],[ 12, -5,-17]],\
                              [[ -7,-17, 17],[-12, -5, 17]],\
                              [[ -7, 17,-17],[-12,  5,-17]],\
                              [[ 17,  7, 17],[  5, 12, 17]],\
                              [[-17,  7,-17],[ -5, 12,-17]],\
                              [[-17, -7, 17],[ -5,-12, 17]],\
                              [[ 17, -7,-17],[  5,-12,-17]]]),
              'NW': \
                    np.array([[[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[  1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[ -1,  1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1, -1,  1],[  0,  1,  1]],\
                              [[  1,  1, -1],[  0,  1,  1]],\
                              [[  1,  1, -1],[  0,  1,  1]],\
                              [[  1,  1, -1],[  0,  1,  1]]]),
              'Bain': \
                    np.array([[[  1,  0,  0],[  1,  0,  0]],\
                              [[  0,  1,  0],[  0,  1,  0]],\
                              [[  0,  0,  1],[  0,  0,  1]]]),
              }

    normals = {'KS': \
                    np.array([[[ -1,  0,  1],[ -1, -1,  1]],\
                              [[ -1,  0,  1],[ -1,  1, -1]],\
                              [[  0,  1, -1],[ -1, -1,  1]],\
                              [[  0,  1, -1],[ -1,  1, -1]],\
                              [[  1, -1,  0],[ -1, -1,  1]],\
                              [[  1, -1,  0],[ -1,  1, -1]],\
                              [[ -1,  0,  1],[ -1, -1,  1]],\
                              [[ -1,  0,  1],[ -1,  1, -1]],\
                              [[  0,  1, -1],[ -1, -1,  1]],\
                              [[  0,  1, -1],[ -1,  1, -1]],\
                              [[  1, -1,  0],[ -1, -1,  1]],\
                              [[  1, -1,  0],[ -1,  1, -1]],\
                              [[ -1,  0,  1],[ -1, -1,  1]],\
                              [[ -1,  0,  1],[ -1,  1, -1]],\
                              [[  0,  1, -1],[ -1, -1,  1]],\
                              [[  0,  1, -1],[ -1,  1, -1]],\
                              [[  1, -1,  0],[ -1, -1,  1]],\
                              [[  1, -1,  0],[ -1,  1, -1]],\
                              [[ -1,  0,  1],[ -1, -1,  1]],\
                              [[ -1,  0,  1],[ -1,  1, -1]],\
                              [[  0,  1, -1],[ -1, -1,  1]],\
                              [[  0,  1, -1],[ -1,  1, -1]],\
                              [[  1, -1,  0],[ -1, -1,  1]],\
                              [[  1, -1,  0],[ -1,  1, -1]]]),
              'GT': \
                    np.array([[[ 17, -5,-12],[ 17,-17, -7]],\
                              [[ 17,-12, -5],[ 17, -7,-17]],\
                              [[-17,  5,-12],[-17, 17, -7]],\
                              [[-17, 12, -5],[-17,  7,-17]],\
                              [[ 17,  5, 12],[ 17, 17,  7]],\
                              [[ 17, 12,  5],[ 17,  7, 17]],\
                              [[-17, -5, 12],[-17,-17,  7]],\
                              [[-17,-12,  5],[-17, -7, 17]],\
                              [[-12, 17, -5],[ -7, 17,-17]],\
                              [[ -5, 17,-12],[-17, 17, -7]],\
                              [[ 12,-17, -5],[  7,-17,-17]],\
                              [[  5,-17,-12],[ 17,-17, -7]],\
                              [[-12,-17,  5],[ -7,-17, 17]],\
                              [[ -5,-17, 12],[-17,-17,  7]],\
                              [[ 12, 17,  5],[  7, 17, 17]],\
                              [[  5, 17, 12],[ 17, 17,  7]],\
                              [[ -5,-12, 17],[-17, -7, 17]],\
                              [[-12, -5, 17],[ -7,-17, 17]],\
                              [[  5, 12, 17],[ 17,  7, 17]],\
                              [[ 12,  5, 17],[  7, 17, 17]],\
                              [[ -5, 12,-17],[-17,  7,-17]],\
                              [[-12,  5,-17],[ -7, 17,-17]],\
                              [[  5,-12,-17],[ 17, -7,-17]],\
                              [[ 12, -5,-17],[  7,-17,-17]]]),
              "GT'": \
                    np.array([[[ -1,  0,  1],[ -1,  1,  1]],\
                              [[ -1,  0,  1],[ -1, -1,  1]],\
                              [[  1,  0,  1],[  1, -1,  1]],\
                              [[  1,  0,  1],[  1,  1,  1]],\
                              [[ -1,  1,  0],[ -1,  1,  1]],\
                              [[ -1,  1,  0],[ -1,  1, -1]],\
                              [[  1,  1,  0],[  1,  1,  1]],\
                              [[  1,  1,  0],[  1,  1, -1]],\
                              [[  1, -1,  0],[  1, -1,  1]],\
                              [[  1, -1,  0],[  1, -1, -1]],\
                              [[ -1, -1,  0],[ -1, -1,  1]],\
                              [[ -1, -1,  0],[ -1, -1, -1]],\
                              [[  0, -1,  1],[  1, -1,  1]],\
                              [[  0, -1,  1],[ -1, -1,  1]],\
                              [[  0,  1,  1],[ -1,  1,  1]],\
                              [[  0,  1,  1],[  1,  1,  1]],\
                              [[  0,  1, -1],[  1,  1, -1]],\
                              [[  0,  1, -1],[ -1,  1, -1]],\
                              [[  0, -1, -1],[ -1, -1, -1]],\
                              [[  0, -1, -1],[  1, -1, -1]],\
                              [[  1,  0, -1],[  1,  1, -1]],\
                              [[  1,  0, -1],[  1, -1, -1]],\
                              [[ -1,  0, -1],[ -1, -1, -1]],\
                              [[ -1,  0, -1],[ -1,  1, -1]]]),
              'NW': \
                    np.array([[[  1, -1,  0],[  1,  0,  0]],\
                              [[  1,  0, -1],[  1,  0,  0]],\
                              [[  0, -1,  1],[  1,  0,  0]],\
                              [[  1,  1,  0],[  1,  0,  0]],\
                              [[  0,  1, -1],[  1,  0,  0]],\
                              [[  1,  0,  1],[  1,  0,  0]],\
                              [[  1,  1,  0],[  1,  0,  0]],\
                              [[  0,  1,  1],[  1,  0,  0]],\
                              [[ -1,  0,  1],[  1,  0,  0]],\
                              [[  1,  0,  1],[  1,  0,  0]],\
                              [[  1, -1,  0],[  1,  0,  0]],\
                              [[  0,  1,  1],[  1,  0,  0]]]),
              'Bain': \
                    np.array([[[  0,  1,  0],[  0,  1,  1]],
                              [[  0,  0,  1],[  1,  0,  1]],
                              [[  1,  0,  0],[  1,  1,  0]]])]
              }
    myMatrix = np.array([[planes [relationModel][variant,me]],\
                         [normals[relationModel][variant,me]],\
                         [np.cross(normals[relationModel][variant,me],planes[relationModel][variant,me])]])
    otherMatrix = np.array([[planes [relationModel][variant,other]],\
                            [normals[relationModel][variant,other]],\
                            [np.cross(normals[relationModel][variant,other],planes[relationModel][variant,other])]])
  
  def getRotation(self,variant):
    fccN=np.array([1.,1.,1.])
    fccN=fccN/np.linalg.norm(fccN)
    fccN=self.orientation.asMatrix().dot(fccN)

    fccD=np.array([17.,-5.,-12.])
    fccD=fccD/np.linalg.norm(fccD)
    fccD=self.orientation.asMatrix().dot(fccD)


    bccN=np.array([1.,1.,0.])
    bccN=bccN/np.linalg.norm(bccN)
    bccN=self.orientation.asMatrix().dot(bccN)

    bccD=np.array([17.,-17.,-7.])
    bccD=bccD/np.linalg.norm(bccD)
    bccD=self.orientation.asMatrix().dot(bccD)

    B = np.array(np.outer(bccN,fccN.T)+np.outer(bccD,fccD.T))*0.5
    U,S,VT = np.linalg.svd(B)
    M=np.diag([1,1,np.linalg.det(U)*np.linalg.det(VT)])
    R=(U.dot(M)).dot(VT)
    return Orientation(matrix=R)
    
     

