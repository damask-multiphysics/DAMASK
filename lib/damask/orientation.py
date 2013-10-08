# code derived from http://pyeuclid.googlecode.com/svn/trunk/euclid.py

import math,random

# ******************************************************************************************
class Vector3:
# ******************************************************************************************
    __slots__ = ['x', 'y', 'z']

    def __init__(self, array):
        assert hasattr(array, '__len__') and len(array) == 3
        self.x = array[0]
        self.y = array[1]
        self.z = array[2]

    def __copy__(self):
        return self.__class__(self.x, self.y, self.z)

    copy = __copy__

    def __repr__(self):
        return 'Vector3(%.4f, %.4f, %.4f)' % (self.x,
                                              self.y,
                                              self.z)

    def __eq__(self, other):
        if isinstance(other, Vector3):
            return self.x == other.x and \
                   self.y == other.y and \
                   self.z == other.z
        else:
            assert hasattr(other, '__len__') and len(other) == 3
            return self.x == other[0] and \
                   self.y == other[1] and \
                   self.z == other[2]

    def __neq__(self, other):
        return not self.__eq__(other)

    def __cmp__(self,other):
        return cmp((self.z,self.y,self.x),(other.z,other.y,other.x))

    def __nonzero__(self):
        return self.x != 0 or self.y != 0 or self.z != 0

    def __len__(self):
        return 3

    def __getitem__(self, key):
        return (self.x, self.y, self.z)[key]

    def __setitem__(self, key, value):
        l = [self.x, self.y, self.z]
        l[key] = value
        self.x, self.y, self.z = l

    def __iter__(self):
        return iter((self.x, self.y, self.z))

    def __getattr__(self, name):
        try:
            return tuple([(self.x, self.y, self.z)['xyz'.index(c)] \
                          for c in name])
        except ValueError:
            raise AttributeError, name

    def __add__(self, other):
        if isinstance(other, Vector3):
            return Vector3([self.x + other.x,
                            self.y + other.y,
                            self.z + other.z])
        else:
            assert hasattr(other, '__len__') and len(other) == 3
            return Vector3([self.x + other[0],
                            self.y + other[1],
                            self.z + other[2]])
    __radd__ = __add__

    def __iadd__(self, other):
        if isinstance(other, Vector3):
            self.x += other.x
            self.y += other.y
            self.z += other.z
        else:
            assert hasattr(other, '__len__') and len(other) == 3
            self.x += other[0]
            self.y += other[1]
            self.z += other[2]
        return self

    def __sub__(self, other):
        if isinstance(other, Vector3):
            return Vector3([self.x - other.x,
                            self.y - other.y,
                            self.z - other.z])
        else:
            assert hasattr(other, '__len__') and len(other) == 3
            return Vector3([self.x - other[0],
                            self.y - other[1],
                            self.z - other[2]])

   
    def __rsub__(self, other):
        if isinstance(other, Vector3):
            return Vector3([other.x - self.x,
                            other.y - self.y,
                            other.z - self.z])
        else:
            assert hasattr(other, '__len__') and len(other) == 3
            return Vector3([other.x - self[0],
                            other.y - self[1],
                            other.z - self[2]])

    def __mul__(self, other):
        if isinstance(other, Vector3):
            # TODO component-wise mul/div in-place and on Vector2; docs.
            return Vector3([self.x * other.x,
                            self.y * other.y,
                            self.z * other.z])
        else: 
            assert type(other) in (int, long, float)
            return Vector3([self.x * other,
                            self.y * other,
                            self.z * other])

    __rmul__ = __mul__

    def __imul__(self, other):
        assert type(other) in (int, long, float)
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def __div__(self, other):
        assert type(other) in (int, long, float)
        self.x /= other
        self.y /= other
        self.z /= other
        return self


    def __rdiv__(self, other):
        assert type(other) in (int, long, float)
        return Vector3([operator.div(other, self.x),
                        operator.div(other, self.y),
                        operator.div(other, self.z)])

    def __floordiv__(self, other):
        assert type(other) in (int, long, float)
        return Vector3([operator.floordiv(self.x, other),
                        operator.floordiv(self.y, other),
                        operator.floordiv(self.z, other)])


    def __rfloordiv__(self, other):
        assert type(other) in (int, long, float)
        return Vector3([operator.floordiv(other, self.x),
                        operator.floordiv(other, self.y),
                        operator.floordiv(other, self.z)])

    def __truediv__(self, other):
        assert type(other) in (int, long, float)
        return Vector3([operator.truediv(self.x, other),
                        operator.truediv(self.y, other),
                        operator.truediv(self.z, other)])


    def __rtruediv__(self, other):
        assert type(other) in (int, long, float)
        return Vector3([operator.truediv(other, self.x),
                        operator.truediv(other, self.y),
                        operator.truediv(other, self.z)])
    
    def __neg__(self):
        return Vector3([-self.x,
                        -self.y,
                        -self.z])

    __pos__ = __copy__
    
    def __abs__(self):
        return math.sqrt(self.x ** 2 + \
                         self.y ** 2 + \
                         self.z ** 2)

    magnitude = __abs__

    def magnitude_squared(self):
        return self.x ** 2 + \
               self.y ** 2 + \
               self.z ** 2

    def normalize(self):
        d = self.magnitude()
        if d:
            self.x /= d
            self.y /= d
            self.z /= d
        return self

    def normalized(self):
        d = self.magnitude()
        if d:
            return Vector3([self.x / d, 
                            self.y / d, 
                            self.z / d])
        return self.copy()

    def dot(self, other):
      if isinstance(other, Vector3):
        return self.x * other.x + \
               self.y * other.y + \
               self.z * other.z
      else:
        assert hasattr(other, '__len__') and len(other) == 3
        return self.x * other[0] + \
               self.y * other[1] + \
               self.z * other[2]

    def cross(self, other):
      if isinstance(other, Vector3):
        return Vector3([  self.y * other.z - self.z * other.y,
                        - self.x * other.z + self.z * other.x,
                          self.x * other.y - self.y * other.x])
      else:
        assert hasattr(other, '__len__') and len(other) == 3
        return Vector3([  self.y * other[2] - self.z * other[1],
                        - self.x * other[2] + self.z * other[0],
                          self.x * other[1] - self.y * other[0]])

    def reflect(self, normal):
        # assume normal is normalized
      if isinstance(other, Vector3):
        d = 2 * (self.x * normal.x + self.y * normal.y + self.z * normal.z)
        return Vector3([self.x - d * normal.x,
                        self.y - d * normal.y,
                        self.z - d * normal.z])
      else:
        assert hasattr(other, '__len__') and len(other) == 3
        d = 2 * (self.x * normal[0] + self.y * normal[1] + self.z * normal[2])
        return Vector3([self.x - d * normal[0],
                        self.y - d * normal[1],
                        self.z - d * normal[2]])

    def inSST(self,symmetry):
      '''
      Determination of disorientations follow the work of A. Heinz and P. Neumann:
      Representation of Orientation and Disorientation Data for
      Cubic, Hexagonal, Tetragonal and Orthorhombic Crystals
      Acta Cryst. (1991). A47, 780-789
      '''
      return {'cubic':       lambda R: R.inFZ('cubic')       and R.x >= R.y and R.y >= R.z and R.z >= 0.0,
              'hexagonal':   lambda R: R.inFZ('hexagonal')   and R.x >= math.sqrt(3)*R.y and R.y >= 0.0 and R.z >= 0.0,
              'tetragonal':  lambda R: R.inFZ('tetragonal')  and R.x >= R.y and R.y >= 0.0 and R.z >= 0.0,
              'orthorombic': lambda R: R.inFZ('orthorombic') and R.x >= 0.0 and R.y >= 0.0 and R.z >= 0.0,
             }.get(symmetry.lower(),False)(self)

    def inFZ(self,symmetry):
        return {'cubic':       lambda R:     math.sqrt(2.0)-1.0 >= abs(R.x) \
                                         and math.sqrt(2.0)-1.0 >= abs(R.y) \
                                         and math.sqrt(2.0)-1.0 >= abs(R.z) \
                                         and 1.0 >= abs(R.x) + abs(R.y) + abs(R.z),
                'hexagonal':   lambda R:     1.0 >= abs(R.x) and 1.0 >= abs(R.y) and 1.0 >= abs(R.z) \
                                         and 2.0 >= math.sqrt(3)*abs(R.x) + abs(R.y) \
                                         and 2.0 >= math.sqrt(3)*abs(R.y) + abs(R.x) \
                                         and 2.0 >= math.sqrt(3) + abs(R.z),
                'tetragonal':  lambda R:     1.0 >= abs(R.x) and 1.0 >= abs(R.y) \
                                         and math.sqrt(2.0) >= abs(R.x) + abs(R.y) \
                                         and math.sqrt(2.0) >= abs(R.z) + 1.0,
                'orthorombic': lambda R:     1.0 >= abs(R.x) and 1.0 >= abs(R.y) and 1.0 >= abs(R.z),
               }.get(symmetry.lower(),False)(self)

# ******************************************************************************************
class Quaternion:
# ******************************************************************************************
    # All methods and naming conventions based off 
    # http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions

    # w is the real part, (x, y, z) are the imaginary parts
    __slots__ = ['w', 'x', 'y', 'z']

    def __init__(self, quatArray=[1.0,0.0,0.0,0.0]):
        self.w, \
        self.x, \
        self.y, \
        self.z = quatArray

    def __iter__(self):
        return iter([self.w,self.x,self.y,self.z])
        
    def __copy__(self):
        Q = Quaternion([self.w,self.x,self.y,self.z])
        return Q

    copy = __copy__

    def __repr__(self):
        return 'Quaternion(real=%.4f, imag=<%.4f, %.4f, %.4f>)' % \
            (self.w, self.x, self.y, self.z)

    def __mul__(self, other):
        if isinstance(other, Quaternion):
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
        elif isinstance(other, Vector3):    # active rotation of vector (i.e. within fixed coordinate system)
            w = self.w
            x = self.x
            y = self.y
            z = self.z
            Vx = other.x
            Vy = other.y
            Vz = other.z
                
            return other.__class__(\
               w * w * Vx + 2 * y * w * Vz - 2 * z * w * Vy + \
               x * x * Vx + 2 * y * x * Vy + 2 * z * x * Vz - \
               z * z * Vx - y * y * Vx,
               2 * x * y * Vx + y * y * Vy + 2 * z * y * Vz + \
               2 * w * z * Vx - z * z * Vy + w * w * Vy - \
               2 * x * w * Vz - x * x * Vy,
               2 * x * z * Vx + 2 * y * z * Vy + \
               z * z * Vz - 2 * w * y * Vx - y * y * Vz + \
               2 * w * x * Vy - x * x * Vz + w * w * Vz)
        elif isinstance(other, (int,float,long)):
            Q = self.copy()
            Q.w *= other
            Q.x *= other
            Q.y *= other
            Q.z *= other
            return Q
        else:
            return self.copy()

    def __imul__(self, other):
        if isinstance(other, Quaternion):
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
        return self

    def __div__(self, other):
        if isinstance(other, (int,float,long)):
            Q = self.copy()
            Q.w /= other
            Q.x /= other
            Q.y /= other
            Q.z /= other
            return Q
        else:
            return self.copy()

    def __idiv__(self, other):
        if isinstance(other, (int,float,long)):
            self.w /= other
            self.x /= other
            self.y /= other
            self.z /= other
        return self

    def __div__(self, other):
        if isinstance(other, (int,float,long)):
            Q = self.copy()
            Q.w /= other
            Q.x /= other
            Q.y /= other
            Q.z /= other
            return Q
        else:
            return self.copy()

    def __idiv__(self, other):
        if isinstance(other, (int,float,long)):
            self.w /= other
            self.x /= other
            self.y /= other
            self.z /= other
        return self

    def __add__(self, other):
        if isinstance(other, Quaternion):
            Q = self.copy()
            Q.w += other.w
            Q.x += other.x
            Q.y += other.y
            Q.z += other.z
            return Q
        else:
            return self.copy()

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
        self.w = 1
        self.x = 0
        self.y = 0
        self.z = 0
        return self

    def rotateBy_angleaxis(self, angle, axis):
        self *= Quaternion.from_angleaxis(angle, axis)
        return self

    def rotateBy_Eulers(self, heading, attitude, bank):
        self *= Quaternion.from_Eulers(eulers, type)
        return self

    def rotateBy_matrix(self, m):
        self *= Quaternion.from_matrix(m)
        return self

    def normalize(self):
        d = self.magnitude()
        if d != 0:
            self /= d
        return self

    def conjugate(self):
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z
        return self

    def inverse(self):
        d = self.magnitude()
        if d != 0:
          self.conjugate()
          self /= d
        return self

    def reduce(self, symmetry='cubic'):
        for q in self.symmetricEquivalents(symmetry):
          if q.inFZ(symmetry): break
        self.w = q.w
        self.x = q.x
        self.y = q.y
        self.z = q.z
        return self

    def normalized(self):
        return self.copy().normalize()

    def conjugated(self):
        return self.copy().conjugate()

    def inversed(self):
        return self.copy().inverse()

    def reduced(self):
        return self.copy().reduce()

    def angleaxis(self):
        if self.w > 1:
            self.normalize()
        angle = 2 * math.acos(self.w)
        s = math.sqrt(1 - self.w ** 2)
        if s < 0.001:
            return angle, Vector3([1, 0, 0])
        else:
            return angle, Vector3([self.x / s, self.y / s, self.z / s])

    def Rodrigues(self):
        if self.w != 0.0:
            return Vector3([self.x / self.w, self.y / self.w, self.z / self.w])
        else:
            return Vector3([float('inf')]*3)

    def Eulers(self,type='bunge'):
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

    def inSST(self,symmetry):
      return self.Rodrigues().inSST(symmetry)

    def inFZ(self,symmetry):
      return self.Rodrigues().inFZ(symmetry)

    def symmetricEquivalents(self,symmetry):
      symQuats = {'cubic':
                      [
                        [ 1.0,0.0,0.0,0.0 ],
                        [ 0.0,-1.0,0.0,0.0 ],
                        [ 0.0,0.0,-1.0,0.0 ],
                        [ 0.0,0.0,0.0,1.0 ],
                        [ 0.0,0.0,0.5*math.sqrt(2), 0.5*math.sqrt(2) ],
                        [ 0.0,0.0,0.5*math.sqrt(2),-0.5*math.sqrt(2) ],
                        [ 0.0, 0.5*math.sqrt(2),-0.5*math.sqrt(2),0.0 ],
                        [ 0.0,-0.5*math.sqrt(2),-0.5*math.sqrt(2),0.0 ],
                        [ 0.0, 0.5*math.sqrt(2), 0.0, 0.5*math.sqrt(2) ],
                        [ 0.0, 0.5*math.sqrt(2), 0.0,-0.5*math.sqrt(2) ],
                        [  0.5, 0.5, 0.5, 0.5 ],
                        [ -0.5, 0.5, 0.5, 0.5 ],
                        [ -0.5, 0.5, 0.5,-0.5 ],
                        [ -0.5, 0.5,-0.5, 0.5 ],
                        [ -0.5,-0.5, 0.5, 0.5 ],
                        [ -0.5,-0.5, 0.5,-0.5 ],
                        [ -0.5,-0.5,-0.5, 0.5 ],
                        [ -0.5, 0.5,-0.5,-0.5 ],
                        [ -0.5*math.sqrt(2),0.0,0.0,0.5*math.sqrt(2) ],
                        [  0.5*math.sqrt(2),0.0,0.0,0.5*math.sqrt(2) ],
                        [ -0.5*math.sqrt(2),0.0,-0.5*math.sqrt(2),0.0 ],
                        [ -0.5*math.sqrt(2),0.0, 0.5*math.sqrt(2),0.0 ],
                        [ -0.5*math.sqrt(2), 0.5*math.sqrt(2),0.0,0.0 ],
                        [ -0.5*math.sqrt(2),-0.5*math.sqrt(2),0.0,0.0 ],
                      ],
                  'hexagonal':
                      [
                        [ 1.0,0.0,0.0,0.0],
                        [ 0.0,1.0,0.0,0.0],
                        [ 0.0,0.0,1.0,0.0],
                        [ 0.0,0.0,0.0,1.0],
                        [ 0.0, 0.5,-0.5*math.sqrt(3),0.0],
                        [ 0.0,-0.5,-0.5*math.sqrt(3),0.0],
                        [ 0.5,0.0,0.0,0.5*math.sqrt(3)],
                        [-0.5,0.0,0.0,0.5*math.sqrt(3)],
                        [-0.5*math.sqrt(3),0.0,0.0, 0.5],
                        [-0.5*math.sqrt(3),0.0,0.0,-0.5],
                        [0.0, 0.5*math.sqrt(3), 0.5,0.0],
                        [0.0,-0.5*math.sqrt(3), 0.5,0.0],
                      ],
                  'tetragonal':
                      [
                        [ 0.0,0.0,1.0,0.0 ],
                        [ 0.0, 0.5*math.sqrt(2),0.5*math.sqrt(2),0.0 ],
                        [ 0.0,-0.5*math.sqrt(2),0.5*math.sqrt(2),0.0 ],
                        [ 0.0,0.0,0.0,1.0 ],
                        [ 0.0,1.0,0.0,0.0 ],
                        [  0.5*math.sqrt(2),0.0,0.0,0.5*math.sqrt(2) ],
                        [ -0.5*math.sqrt(2),0.0,0.0,0.5*math.sqrt(2) ],
                        [ 1.0,0.0,0.0,0.0 ],
                      ],
                  'orthorombic':
                      [
                        [ 0.0,0.0,1.0,0.0 ],
                        [ 0.0,0.0,0.0,1.0 ],
                        [ 0.0,1.0,0.0,0.0 ],
                        [ 1.0,0.0,0.0,0.0 ],
                      ],
               }
      equivalents = [self*Quaternion(q) for q in symQuats[symmetry.lower()]]
      return equivalents


    def disorientation(self,other,symmetry):

      QinSST = Quaternion()
      breaker = False
      for me in self.symmetricEquivalents(symmetry):
        me.inverse()
        for they in other.symmetricEquivalents(symmetry):
          theQ = me * they
          if theQ.w < 0.0: theQ.w = -theQ.w
          if theQ.inSST(symmetry):
            QinSST = theQ.copy()
            breaker = True
            break
        if breaker: break
        
      return QinSST

#    # Static constructors
    def new_identity(cls):
        return cls()
    new_identity = classmethod(new_identity)

    def from_random(cls):
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        Q = cls()
        Q.w = math.cos(2.0*math.pi*r1)*math.sqrt(r3)
        Q.x = math.sin(2.0*math.pi*r2)*math.sqrt(1.0-r3)
        Q.y = math.cos(2.0*math.pi*r2)*math.sqrt(1.0-r3)
        Q.z = math.sin(2.0*math.pi*r1)*math.sqrt(r3)
        return Q
    from_random = classmethod(from_random)

    def from_Rodrigues(cls, rodrigues):
        if not isinstance(rodrigues, Vector3): rodrigues = Vector3(rodrigues)
        halfangle = math.atan(rodrigues.magnitude())
        c = math.cos(halfangle)
        Q = cls()
        Q.w = c
        Q.x = rodrigues.x * c
        Q.y = rodrigues.y * c
        Q.z = rodrigues.z * c
        return Q
    from_Rodrigues = classmethod(from_Rodrigues)

    def from_angleaxis(cls, angle, axis):
        assert(isinstance(axis, Vector3))
        axis = axis.normalized()
        s = math.sin(angle / 2)
        Q = cls()
        Q.w = math.cos(angle / 2)
        Q.x = axis.x * s
        Q.y = axis.y * s
        Q.z = axis.z * s
        return Q
    from_angleaxis = classmethod(from_angleaxis)

    def from_Eulers(cls, eulers, type = 'Bunge'):
        c1 = math.cos(eulers[0] / 2.0)
        s1 = math.sin(eulers[0] / 2.0)
        c2 = math.cos(eulers[1] / 2.0)
        s2 = math.sin(eulers[1] / 2.0)
        c3 = math.cos(eulers[2] / 2.0)
        s3 = math.sin(eulers[2] / 2.0)
        Q = cls()

        if type.lower() == 'bunge' or type.lower() == 'zxz':
          Q.w =   c1 * c2 * c3 - s1 * c2 * s3
          Q.x =   c1 * s2 * c3 + s1 * s2 * s3
          Q.y = - c1 * s2 * s3 + s1 * s2 * c3
          Q.z =   c1 * c2 * s3 + s1 * c2 * c3
        else:
          print 'unknown Euler convention'
          Q.w = c1 * c2 * c3 - s1 * s2 * s3
          Q.x = s1 * s2 * c3 + c1 * c2 * s3
          Q.y = s1 * c2 * c3 + c1 * s2 * s3
          Q.z = c1 * s2 * c3 - s1 * c2 * s3
        return Q
    from_Eulers = classmethod(from_Eulers)
    
    def from_matrix(cls, m):
      if m[0*4 + 0] + m[1*4 + 1] + m[2*4 + 2] > 0.00000001:
        t = m[0*4 + 0] + m[1*4 + 1] + m[2*4 + 2] + 1.0
        s = 0.5/math.sqrt(t)
        
        return cls(
          s*t,
          (m[1*4 + 2] - m[2*4 + 1])*s,
          (m[2*4 + 0] - m[0*4 + 2])*s,
          (m[0*4 + 1] - m[1*4 + 0])*s
          )
        
      elif m[0*4 + 0] > m[1*4 + 1] and m[0*4 + 0] > m[2*4 + 2]:
        t = m[0*4 + 0] - m[1*4 + 1] - m[2*4 + 2] + 1.0
        s = 0.5/math.sqrt(t)
        
        return cls(
          (m[1*4 + 2] - m[2*4 + 1])*s,
          s*t,
          (m[0*4 + 1] + m[1*4 + 0])*s,
          (m[2*4 + 0] + m[0*4 + 2])*s
          )
        
      elif m[1*4 + 1] > m[2*4 + 2]:
        t = -m[0*4 + 0] + m[1*4 + 1] - m[2*4 + 2] + 1.0
        s = 0.5/math.sqrt(t)
        
        return cls(
          (m[2*4 + 0] - m[0*4 + 2])*s,
          (m[0*4 + 1] + m[1*4 + 0])*s,
          s*t,
          (m[1*4 + 2] + m[2*4 + 1])*s
          )
        
      else:
        t = -m[0*4 + 0] - m[1*4 + 1] + m[2*4 + 2] + 1.0
        s = 0.5/math.sqrt(t)
        
        return cls(
          (m[0*4 + 1] - m[1*4 + 0])*s,
          (m[2*4 + 0] + m[0*4 + 2])*s,
          (m[1*4 + 2] + m[2*4 + 1])*s,
          s*t
          )
    from_matrix = classmethod(from_matrix)
    
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
    new_interpolate = classmethod(new_interpolate)
