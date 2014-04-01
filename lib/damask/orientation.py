# -*- coding: UTF-8 no BOM -*-

import numpy,math,random

# ******************************************************************************************
class Rodrigues:
# ******************************************************************************************

    def __init__(self, vector = numpy.zeros(3)):
      self.vector = vector
    
    def asQuaternion(self):
      norm = numpy.linalg.norm(self.vector)
      halfAngle = numpy.arctan(norm)
      return Quaternion(numpy.cos(halfAngle),numpy.sin(halfAngle)*self.vector/norm)

    def asAngleAxis(self):
      norm = numpy.linalg.norm(self.vector)
      halfAngle = numpy.arctan(norm)
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

          return numpy.array([\
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
      try:                                                          # Quaternion
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

    def asMatrix(self):
      return numpy.array([[1.0-2.0*(self.y*self.y+self.z*self.z),     2.0*(self.x*self.y-self.z*self.w),     2.0*(self.x*self.z+self.y*self.w)],
                          [    2.0*(self.x*self.y+self.z*self.w), 1.0-2.0*(self.x*self.x+self.z*self.z),     2.0*(self.y*self.z-self.x*self.w)],
                          [    2.0*(self.x*self.z-self.y*self.w),     2.0*(self.x*self.w+self.y*self.z), 1.0-2.0*(self.x*self.x+self.y*self.y)]])
    
    def asAngleAxis(self):
      if self.w > 1:
          self.normalize()
      angle = 2 * math.acos(self.w)
      s = math.sqrt(1 - self.w ** 2)
      if s < 0.001:
          return angle, nunmpy.array([1.0, 0.0, 0.0])
      else:
          return angle, numpy.array([self.x / s, self.y / s, self.z / s])

    def asRodrigues(self):
      if self.w != 0.0:
          return numpy.array([self.x, self.y, self.z])/self.w
      else:
          return numpy.array([float('inf')]*3)

    def asEulers(self,type='bunge'):
      angles = [0.0,0.0,0.0]
      if type.lower() == 'bunge' or type.lower() == 'zxz':
        angles[0] =  math.atan2( self.x*self.z+self.y*self.w,
                                -self.y*self.z+self.x*self.w)
#          angles[1] =  math.acos(-self.x**2-self.y**2+self.z**2+self.w**2)
        angles[1] =  math.acos(1.0 - 2*(self.x**2+self.y**2))
        angles[2] =  math.atan2( self.x*self.z-self.y*self.w,
                                +self.y*self.z+self.x*self.w)

        if angles[0] < 0.0:
          angles[0] += 2*math.pi
        if angles[1] < 0.0:
          angles[1] += math.pi
          angles[2] *= -1
        if angles[2] < 0.0:
          angles[2] += 2*math.pi
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
        if not isinstance(rodrigues, numpy.ndarray): rodrigues = numpy.array(rodrigues)
        halfangle = math.atan(numpy.linalg.norm(rodrigues))
        c = math.cos(halfangle)
        w = c
        x,y,z = c*rodrigues
        return cls([w,x,y,z])


    @classmethod
    def fromAngleAxis(cls, angle, axis):
        if not isinstance(axis, numpy.ndarray): axis = numpy.array(axis)
        axis /= numpy.linalg.norm(axis)
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
          print 'unknown Euler convention'
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
#     basis = {'cubic' :        numpy.linalg.inv(numpy.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,1.]/numpy.sqrt(2.),                     # direction of green
#                                                             [1.,1.,1.]/numpy.sqrt(3.)]).transpose()),      # direction of blue
#              'hexagonal' :    numpy.linalg.inv(numpy.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,0.],                                    # direction of green
#                                                 [numpy.sqrt(3.),1.,0.]/numpy.sqrt(4.)]).transpose()),      # direction of blue
#              'tetragonal' :   numpy.linalg.inv(numpy.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,0.],                                    # direction of green
#                                                             [1.,1.,0.]/numpy.sqrt(2.)]).transpose()),      # direction of blue
#              'orthorhombic' : numpy.linalg.inv(numpy.array([[0.,0.,1.],                                    # direction of red
#                                                             [1.,0.,0.],                                    # direction of green
#                                                             [0.,1.,0.]]).transpose()),                     # direction of blue
#             }
    if self.lattice == 'cubic':
      basis = numpy.array([ [-1.            ,  0.            ,  1. ],
                            [ numpy.sqrt(2.), -numpy.sqrt(2.),  0. ],
                            [ 0.            ,  numpy.sqrt(3.),  0. ] ])
    elif self.lattice == 'hexagonal':
      basis = numpy.array([ [ 0.            ,  0.            ,  1. ],
                            [ 1.            , -numpy.sqrt(3.),  0. ],
                            [ 0.            ,  2.            ,  0. ] ])
    elif self.lattice == 'tetragonal':
      basis = numpy.array([ [ 0.            ,  0.            ,  1. ],
                            [ 1.            , -1.            ,  0. ],
                            [ 0.            ,  numpy.sqrt(2.),  0. ] ])
    elif self.lattice == 'orthorhombic':
      basis = numpy.array([ [ 0., 0., 1.],
                            [ 1., 0., 0.],
                            [ 0., 1., 0.] ])
    else:
      basis = None

    if basis == None:
      theComponents = -numpy.ones(3,'d')
    else:
      theComponents = numpy.dot(basis,numpy.array([vector[0],vector[1],abs(vector[2])]))
      
    inSST = numpy.all(theComponents >= 0.0)

    if color:                                                                      # have to return color array
      if inSST:
        rgb = numpy.power(theComponents/numpy.linalg.norm(theComponents),0.5)      # smoothen color ramps
        rgb = numpy.minimum(numpy.ones(3,'d'),rgb)                                 # limit to maximum intensity
        rgb /= max(rgb)                                                            # normalize to (HS)V = 1
      else:
        rgb = numpy.zeros(3,'d')
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
    elif isinstance(Eulers, numpy.ndarray) and Eulers.shape == (3,):                   # based on given Euler angles
      self.quaternion = Quaternion.fromEulers(Eulers,'bunge')
    elif isinstance(matrix, numpy.ndarray) and matrix.shape == (3,3):                  # based on given rotation matrix
      self.quaternion = Quaternion.fromMatrix(matrix)
    elif isinstance(angleAxis, numpy.ndarray) and angleAxis.shape == (4,):             # based on given angle and rotation axis
      self.quaternion = Quaternion.fromAngleAxis(angleAxis[0],angleAxis[1:4])
    elif isinstance(Rodrigues, numpy.ndarray) and Rodrigues.shape == (3,):             # based on given Rodrigues vector
      self.quaternion = Quaternion.fromRodrigues(Rodrigues)
    elif isinstance(quaternion, Quaternion):                                           # based on given quaternion
      self.quaternion = quaternion.homomorphed()
    elif isinstance(quaternion, numpy.ndarray) and quaternion.shape == (4,):           # based on given quaternion
      self.quaternion = Quaternion(quaternion).homomorphed()

    self.symmetry = Symmetry(symmetry)

  def __copy__(self):
    return self.__class__(quaternion=self.quaternion,symmetry=self.symmetry.lattice)

  copy = __copy__


  def __repr__(self):
    return 'Symmetry: %s\n' % (self.symmetry) + \
           'Quaternion: %s\n' % (self.quaternion) + \
           'Matrix:\n%s\n' % ( '\n'.join(['\t'.join(map(str,self.asMatrix()[i,:])) for i in range(3)]) ) + \
           'Bunge Eulers: %s' % ('\t'.join(map(lambda x:str(numpy.degrees(x)),self.asEulers('Bunge'))) )

  def asQuaternion(self):
    return self.quaternion.asList()


  def asEulers(self,type):
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
        theQ = (me * they).homomorph()
        if theQ.x < 0.0 or theQ.y < 0.0 or theQ.z < 0.0: theQ.conjugate()                   # speed up scanning since minimum angle is usually found for positive x,y,z
        if lowerSymmetry.inDisorientationSST(theQ.asRodrigues()):
          breaker = True
          break
      if breaker: break
      
    return Orientation(quaternion=theQ,symmetry=self.symmetry.lattice)


  def IPFcolor(self,axis):
    '''
    TSL color of inverse pole figure for given axis
    '''

    color = numpy.zeros(3,'d')

    for i,q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):
      pole = q.conjugated()*axis                               # align crystal direction to axis
      inSST,color = self.symmetry.inSST(pole,color=True)
      if inSST: break

    return color
