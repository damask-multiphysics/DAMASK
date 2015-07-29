#!/usr/bin/env python
# encoding: utf-8
# filename: corientation.pyx

#     __ __ __________  ____  __ ____  ______  ____
#    / //_// ____/ __ \/ __ \/ //_/ / / / __ \/ __ \
#   / ,<  / __/ / / / / / / / ,< / / / / / / / / / /
#  / /| |/ /___/ /_/ / /_/ / /| / /_/ / /_/ / /_/ /
# /_/ |_/_____/_____/\____/_/ |_\____/_____/\____/


######################################################
# This is a Cython implementation of original DAMASK #
# orientation class, mainly for speed improvement.   #
######################################################

import math, random, os
import numpy as np
#cimport numpy as np


##
# This Rodrigues class is odd, not sure if it will function
# properly or not
cdef class Rodrigues:
    """Rodrigues representation of orientation """
    cdef public double[3] r

    def __init__(self, vector):
        if isinstance(vector, Rodrigues):
            self.r[0] = vector.r[0]
            self.r[1] = vector.r[1]
            self.r[2] = vector.r[2]
        else:
            self.r[0] = vector[0]
            self.r[1] = vector[1]
            self.r[2] = vector[2]

    def asQuaternion(self):
        cdef double norm, halfAngle
        cdef double[4] q

        norm = np.linalg.norm(self.vector)
        halfAngle = np.arctan(norm)
        q[0] = np.cos(halfAngle)
        tmp = np.sin(halfAngle)*self.vector/norm
        q[1],q[2],q[3] = tmp[0],tmp[1],tmp[2]

        return Quaternion(q)

    def asAngleAxis(self):
        cdef double norm, halfAngle

        norm = np.linalg.norm(self.vector)
        halfAngle = np.arctan(norm)

        return (2.0*halfAngle,self.vector/norm)


##
# The Quaternion class do the heavy lifting of orientation
# calculation
cdef class Quaternion:
    """ Quaternion representation of orientation """
    # All methods and naming conventions based off
    # http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions
    cdef public double w,x,y,z

    def __init__(self, data):
        """ copy constructor friendly """
        cdef double[4] q

        if isinstance(data, Quaternion):
            q[0] = data.w
            q[1] = data.x
            q[2] = data.y
            q[3] = data.z
        else:
            q[0] = data[0]
            q[1] = data[1]
            q[2] = data[2]
            q[3] = data[3]

        self.Quaternion(q)

    cdef Quaternion(self, double* quatArray):
        """
        @para:  quatArray = <w, x, y, z>
        w is the real part, (x, y, z) are the imaginary parts
        """
        if quatArray[0] < 0:
            self.w = -quatArray[0]
            self.x = -quatArray[1]
            self.y = -quatArray[2]
            self.z = -quatArray[3]
        else:
            self.w = quatArray[0]
            self.x = quatArray[1]
            self.y = quatArray[2]
            self.z = quatArray[3]

    def __copy__(self):
        cdef double[4] q = [self.w,self.x,self.y,self.z]
        return Quaternion(q)

    copy = __copy__

    def __iter__(self):
      return iter([self.w,self.x,self.y,self.z])

    def __repr__(self):
        return 'Quaternion(real={:.4f},imag=<{:.4f},{:.4f}, {:.4f}>)'.format(self.w,
                                                                             self.x,
                                                                             self.y,
                                                                             self.z)

    def __pow__(self, exponent, modulo):
        # declare local var for speed gain
        cdef double omega, vRescale
        cdef double[4] q

        omega = math.acos(self.w)
        vRescale = math.sin(exponent*omega)/math.sin(omega)

        q[0] = math.cos(exponent*omega)
        q[1] = self.x*vRescale
        q[2] = self.y*vRescale
        q[3] = self.z*vRescale
        return Quaternion(q)

    def __ipow__(self, exponent):
        self = self.__pow__(self, exponent, 1.0)
        return self

    def __mul__(self, other):
        # declare local var for speed gain
        cdef double Aw,Ax,Ay,Az,Bw,Bx,By,Bz
        cdef double w,x,y,z,Vx,Vy,Vz
        cdef double[4] q

        # quaternion * quaternion
        try:
            Aw = self.w
            Ax = self.x
            Ay = self.y
            Az = self.z
            Bw = other.w
            Bx = other.x
            By = other.y
            Bz = other.z
            q[0] = - Ax * Bx - Ay * By - Az * Bz + Aw * Bw
            q[1] = + Ax * Bw + Ay * Bz - Az * By + Aw * Bx
            q[2] = - Ax * Bz + Ay * Bw + Az * Bx + Aw * By
            q[3] = + Ax * By - Ay * Bx + Az * Bw + Aw * Bz
            return Quaternion(q)
        except:
            pass
        # vector (perform active rotation, i.e. q*v*q.conjugated)
        try:
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
        except:
            pass
        # quaternion * scalar
        try:
            Q = self.copy()
            Q.w *= other
            Q.x *= other
            Q.y *= other
            Q.z *= other
            return Q
        except:
            return self.copy()

    def __imul__(self, other):
        if isinstance(other, Quaternion):
            self = self.__mul__(other)
            return self
        else:
            return NotImplemented

    def __div__(self, other):
        cdef double[4] q

        if isinstance(other, (int,float,long)):
            q[0] = self.w / other
            q[1] = self.x / other
            q[2] = self.y / other
            q[3] = self.z / other
            return Quaternion(q)
        else:
            NotImplemented

    def __idiv__(self, other):
        self = self.__div__(other)
        return self

    def __add__(self, other):
        cdef double[4] q

        if isinstance(other, Quaternion):
            q[0] = self.w + other.w
            q[1] = self.x + other.x
            q[2] = self.y + other.y
            q[3] = self.z + other.z
            return self.__class__(q)
        else:
            return NotImplemented

    def __iadd__(self, other):
        self = self.__add__(other)
        return self

    def __sub__(self, other):
        cdef double[4] q

        if isinstance(other, Quaternion):
            q[0] = self.w - other.w
            q[1] = self.x - other.x
            q[2] = self.y - other.y
            q[3] = self.z - other.z
            return self.__class__(q)
        else:
            return NotImplemented

    def __isub__(self, other):
        self = self.__sub__(other)
        return self

    def __neg__(self):
        cdef double[4] q

        q[0] = -self.w
        q[1] = -self.x
        q[2] = -self.y
        q[3] = -self.z

        return self.__class__(q)

    def __abs__(self):
        cdef double tmp

        tmp = self.w**2 + self.x**2 + self.y**2 + self.z**2
        tmp = math.sqrt(tmp)
        return tmp

    magnitude = __abs__

    def __richcmp__(Quaternion self, Quaternion other, int op):
        cdef bint tmp

        tmp = (abs(self.w-other.w) < 1e-8 and \
               abs(self.x-other.x) < 1e-8 and \
               abs(self.y-other.y) < 1e-8 and \
               abs(self.z-other.z) < 1e-8) \
               or \
               (abs(-self.w-other.w) < 1e-8 and \
               abs(-self.x-other.x) < 1e-8 and \
               abs(-self.y-other.y) < 1e-8 and \
               abs(-self.z-other.z) < 1e-8)
        if op == 2:    #__eq__
            return tmp
        elif op ==3:   #__ne__
            return not tmp
        else:
            return NotImplemented

    def __cmp__(self,other):
        # not sure if this actually works or not
        return cmp(self.Rodrigues(),other.Rodrigues())

    def magnitude_squared(self):
        cdef double tmp

        tmp = self.w**2 + self.x**2 + self.y**2 + self.z**2
        return tmp

    def identity(self):
        self.w = 1.0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        return self

    def rotateBy_angleaxis(self, angle, axis):
        self *= Quaternion.fromAngleAxis(angle, axis)
        return self

    def rotateBy_Eulers(self, eulers):
        self *= Quaternion.fromEulers(eulers, type)
        return self

    def rotateBy_matrix(self, m):
        self *= Quaternion.fromMatrix(m)
        return self

    def normalize(self):
        cdef double d

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
        cdef double d

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

    # return a copy of me
    def normalized(self):
        cdef Quaternion q

        q = Quaternion(self.normalize())
        return q

    def conjugated(self):
        cdef Quaternion q

        q = Quaternion(self.conjugate())
        return q

    def asList(self):
        cdef double[4] q = [self.w, self.x, self.y, self.z]

        return list(q)

    def asM(self):                                                # to find Averaging Quaternions (see F. Landis Markley et al.)
        return np.outer([i for i in self],[i for i in self])

    def asMatrix(self):
        return np.array([[1.0-2.0*(self.y*self.y+self.z*self.z),     2.0*(self.x*self.y-self.z*self.w),     2.0*(self.x*self.z+self.y*self.w)],
                         [    2.0*(self.x*self.y+self.z*self.w), 1.0-2.0*(self.x*self.x+self.z*self.z),     2.0*(self.y*self.z-self.x*self.w)],
                         [    2.0*(self.x*self.z-self.y*self.w),     2.0*(self.x*self.w+self.y*self.z), 1.0-2.0*(self.x*self.x+self.y*self.y)]])

    def asAngleAxis(self):
        cdef double s,x,y

        if self.w > 1:
            self.normalize()

        s = math.sqrt(1. - self.w**2)
        x = 2*self.w**2 - 1.
        y = 2*self.w * s

        angle = math.atan2(y,x)

        return angle, np.array([1.0, 0.0, 0.0] if angle < 1e-3 else [self.x/s, self.y/s, self.z/s])

    def asRodrigues(self):
        if self.w != 0.0:
            return np.array([self.x, self.y, self.z])/self.w
        else:
            return np.array([float('inf')]*3)

    def asEulers(self,type='bunge',degrees=False):
        """conversion taken from:
            Melcher, A.; Unser, A.; Reichhardt, M.; Nestler, B.; Pötschke, M.; Selzer, M.
            Conversion of EBSD data by a quaternion based algorithm to be used for grain structure simulations
            Technische Mechanik 30 (2010) pp 401--413
        """
        cdef double x,y

        angles = [0.0,0.0,0.0]

        if type.lower() == 'bunge' or type.lower() == 'zxz':
            if abs(self.x - self.y) < 1e-8:
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

        return np.degrees(angles) if degrees else angles

    @staticmethod
    def fromIdentity():
        cdef double[4] q = [1.0, 0.0, 0.0, 0.0]

        return Quaternion(q)

    @staticmethod
    def fromRandom(randomSeed=None):
        cdef double r1,r2,r3
        cdef double[4] q

        if randomSeed == None:
            randomSeed = int(os.urandom(4).encode('hex'), 16)
        random.seed(randomSeed)

        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        q[0] = math.cos(2.0*math.pi*r1)*math.sqrt(r3)
        q[1] = math.sin(2.0*math.pi*r2)*math.sqrt(1.0-r3)
        q[2] = math.cos(2.0*math.pi*r2)*math.sqrt(1.0-r3)
        q[3] = math.sin(2.0*math.pi*r1)*math.sqrt(r3)
        return Quaternion(q)

    @staticmethod
    def fromRodrigues(cls, rodrigues):
        if not isinstance(rodrigues, np.ndarray): rodrigues = np.array(rodrigues)
        halfangle = math.atan(np.linalg.norm(rodrigues))
        c = math.cos(halfangle)
        w = c
        x,y,z = c*rodrigues
        return cls([w,x,y,z])

    @staticmethod
    def fromAngleAxis(cls, angle, axis):
        if not isinstance(axis, np.ndarray): axis = np.array(axis)
        axis /= np.linalg.norm(axis)
        s = math.sin(angle / 2.0)
        w = math.cos(angle / 2.0)
        x = axis[0] * s
        y = axis[1] * s
        z = axis[2] * s
        return cls([w,x,y,z])

    @staticmethod
    def fromEulers(cls, eulers, type = 'Bunge'):
        cdef double c1,s1,c2,s2,c3,s3
        cdef double[4] q

        eulers *= 0.5                                                                  # reduce to half angles

        c1 = math.cos(eulers[0])
        s1 = math.sin(eulers[0])
        c2 = math.cos(eulers[1])
        s2 = math.sin(eulers[1])
        c3 = math.cos(eulers[2])
        s3 = math.sin(eulers[2])

        if type.lower() == 'bunge' or type.lower() == 'zxz':
            q[0] =   c1 * c2 * c3 - s1 * c2 * s3
            q[1] =   c1 * s2 * c3 + s1 * s2 * s3
            q[2] = - c1 * s2 * s3 + s1 * s2 * c3
            q[3] =   c1 * c2 * s3 + s1 * c2 * c3
        else:
            q[0] = c1 * c2 * c3 - s1 * s2 * s3
            q[1] = s1 * s2 * c3 + c1 * c2 * s3
            q[2] = s1 * c2 * c3 + c1 * s2 * s3
            q[3] = c1 * s2 * c3 - s1 * c2 * s3
        return Quaternion(q)

    ## Modified Method to calculate Quaternion from Orientation Matrix, Source: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

    @staticmethod
    def fromMatrix(cls, m):
        # This is a slow implementation
        if m.shape != (3,3) and np.prod(m.shape) == 9:
            m = m.reshape(3,3)

        tr=m[0,0]+m[1,1]+m[2,2]
        if tr > 0.00000001:
            s = math.sqrt(tr + 1.0)*2.0

            return Quaternion(
                [ s*0.25,
                  (m[2,1] - m[1,2])/s,
                  (m[0,2] - m[2,0])/s,
                  (m[1,0] - m[0,1])/s
                ])

        elif m[0,0] > m[1,1] and m[0,0] > m[2,2]:
            t = m[0,0] - m[1,1] - m[2,2] + 1.0
            s = 2.0*math.sqrt(t)

            return cls(
                [ (m[2,1] - m[1,2])/s,
                   s*0.25,
                   (m[0,1] + m[1,0])/s,
                   (m[2,0] + m[0,2])/s,
                ])

        elif m[1,1] > m[2,2]:
            t = -m[0,0] + m[1,1] - m[2,2] + 1.0
            s = 2.0*math.sqrt(t)

            return cls(
              [ (m[0,2] - m[2,0])/s,
                (m[0,1] + m[1,0])/s,
                s*0.25,
                (m[1,2] + m[2,1])/s,
              ])

        else:
            t = -m[0,0] - m[1,1] + m[2,2] + 1.0
            s = 2.0*math.sqrt(t)

            return cls(
              [ (m[1,0] - m[0,1])/s,
                (m[2,0] + m[0,2])/s,
                (m[1,2] + m[2,1])/s,
                s*0.25,
              ])

    @staticmethod
    def new_interpolate(cls, q1, q2, t):
        # see http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872_2007014421.pdf for (another?) way to interpolate quaternions

        assert isinstance(q1, Quaternion) and isinstance(q2, Quaternion)
        Q = Quaternion.fromIdentity()

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

##
# Define lattice_type to make it easier for future
# development
cdef enum lattice_type:
    NONE = 0
    ORTHORHOMBIC= 1
    TETRAGONAL = 2
    HEXAGONAL = 3
    CUBIC = 4
##
# Symmetry class
cdef class Symmetry:
    cdef public lattice_type lattice
    # cdef enum LATTICES:
    #     NONE = 0
    #     ORTHORHOMBIC= 1
    #     TETRAGONAL = 2
    #     HEXAGONAL = 3
    #     CUBIC = 4

    def __init__(self, symmetry):
        if symmetry == 0 or symmetry == None:
            self.lattice = NONE
        elif symmetry == 1 or symmetry == 'orthorhombic':
            self.lattice = ORTHORHOMBIC
        elif symmetry == 2 or symmetry == 'tetragonal':
            self.lattice = TETRAGONAL
        elif symmetry == 3 or symmetry == 'hexagonal':
            self.lattice = HEXAGONAL
        elif symmetry == 4 or symmetry == 'cubic':
            self.lattice = CUBIC
        else:
            self.lattice = NONE

    def __copy__(self):
        return self.__class__(self.lattice)

    copy = __copy__

    def __repr__(self):
        return '{}'.format(self.lattice)

    def __richcmp__(Symmetry self, Symmetry other, int op):
        cdef bint tmp

        tmp = self.lattice == other.lattice
        if op == 2:    #__eq__
            return tmp
        elif op ==3:   #__ne__
            return not tmp
        else:
            return NotImplemented

    def __cmp__(self,other):
        return cmp(self.lattice,other.lattice)

    def equivalentQuaternions(self,quaternion):
        '''
        List of symmetrically equivalent quaternions based on own symmetry.
        '''
        if self.lattice == CUBIC:
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
        elif self.lattice == HEXAGONAL:
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
        elif self.lattice == TETRAGONAL:
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
        elif self.lattice == ORTHORHOMBIC:
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
        # due to the use of list comprehension, the speed grain is quite
        # limited
        return [quaternion*Quaternion(q) for q in symQuats]

    def inFZ(self,R):
        '''
        Check whether given Rodrigues vector falls into fundamental zone of own symmetry.
        '''
        if isinstance(R, Quaternion): R = R.asRodrigues()       # translate accidentially passed quaternion
        R = abs(R)                                              # fundamental zone in Rodrigues space is point symmetric around origin
        if self.lattice == CUBIC:
            return     math.sqrt(2.0)-1.0 >= R[0] \
                   and math.sqrt(2.0)-1.0 >= R[1] \
                   and math.sqrt(2.0)-1.0 >= R[2] \
                   and 1.0 >= R[0] + R[1] + R[2]
        elif self.lattice == HEXAGONAL:
            return     1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2] \
                   and 2.0 >= math.sqrt(3)*R[0] + R[1] \
                   and 2.0 >= math.sqrt(3)*R[1] + R[0] \
                   and 2.0 >= math.sqrt(3) + R[2]
        elif self.lattice == TETRAGONAL:
            return     1.0 >= R[0] and 1.0 >= R[1] \
                   and math.sqrt(2.0) >= R[0] + R[1] \
                   and math.sqrt(2.0) >= R[2] + 1.0
        elif self.lattice == ORTHORHOMBIC:
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
        if isinstance(R, Quaternion): R = R.asRodrigues()       # translate accidentally passed quaternion

        cdef double epsilon = 0.0

        if self.lattice == CUBIC:
            return R[0] >= R[1]+epsilon                and R[1] >= R[2]+epsilon    and R[2] >= epsilon and self.inFZ(R)

        elif self.lattice == HEXAGONAL:
            return R[0] >= math.sqrt(3)*(R[1]+epsilon) and R[1] >= epsilon         and R[2] >= epsilon and self.inFZ(R)

        elif self.lattice == TETRAGONAL:
            return R[0] >= R[1]+epsilon                and R[1] >= epsilon         and R[2] >= epsilon and self.inFZ(R)

        elif self.lattice == ORTHORHOMBIC:
            return R[0] >= epsilon                     and R[1] >= epsilon         and R[2] >= epsilon and self.inFZ(R)

        else:
            return True

    def inSST(self,vector,color = False):
        '''
        Check whether given vector falls into standard stereographic triangle of own symmetry.
        Return inverse pole figure color if requested.
        '''
#     basis = {4 : np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                          [1.,0.,1.]/np.sqrt(2.),                        # direction of green
#                                          [1.,1.,1.]/np.sqrt(3.)]).transpose()),         # direction of blue
#              3 : np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                          [1.,0.,0.],                                    # direction of green
#                                          [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).transpose()),      # direction of blue
#              2 : np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                          [1.,0.,0.],                                    # direction of green
#                                          [1.,1.,0.]/np.sqrt(2.)]).transpose()),         # direction of blue
#              1 : np.linalg.inv(np.array([[0.,0.,1.],                                    # direction of red
#                                          [1.,0.,0.],                                    # direction of green
#                                          [0.,1.,0.]]).transpose()),                     # direction of blue
#             }
        if self.lattice == CUBIC:
            basis = np.array([ [-1.            ,  0.            ,  1. ],
                               [ np.sqrt(2.), -np.sqrt(2.),        0. ],
                               [ 0.            ,  np.sqrt(3.),     0. ] ])
        elif self.lattice == HEXAGONAL:
            basis = np.array([ [ 0.            ,  0.            ,  1. ],
                               [ 1.            , -np.sqrt(3.),     0. ],
                               [ 0.            ,  2.            ,  0. ] ])
        elif self.lattice == TETRAGONAL:
            basis = np.array([ [ 0.            ,  0.            ,  1. ],
                               [ 1.            , -1.            ,  0. ],
                               [ 0.            ,  np.sqrt(2.),     0. ] ])
        elif self.lattice == ORTHORHOMBIC:
            basis = np.array([ [ 0., 0., 1.],
                               [ 1., 0., 0.],
                               [ 0., 1., 0.] ])
        else:
            basis = np.zeros((3,3),dtype=float)

        if np.all(basis == 0.0):
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


##
# Orientation class is a composite class of Symmetry and Quaternion
cdef class Orientation:
    cdef public Quaternion quaternion
    cdef public Symmetry symmetry

    def __init__(self,
                 quaternion = Quaternion.fromIdentity(),
                 Rodrigues  = None,
                 angleAxis  = None,
                 matrix     = None,
                 Eulers     = None,
                 random     = False,                                                     # put any integer to have a fixed seed or True for real random
                 symmetry   = None,
                ):
        if random:                                                                         # produce random orientation
            if isinstance(random, bool ):
                self.quaternion = Quaternion.fromRandom()
            else:
                self.quaternion = Quaternion.fromRandom(randomSeed=random)
        elif isinstance(Eulers, np.ndarray) and Eulers.shape == (3,):                      # based on given Euler angles
            self.quaternion = Quaternion.fromEulers(Eulers,'bunge')
        elif isinstance(matrix, np.ndarray) :                                              # based on given rotation matrix
            self.quaternion = Quaternion.fromMatrix(matrix)
        elif isinstance(angleAxis, np.ndarray) and angleAxis.shape == (4,):                # based on given angle and rotation axis
            self.quaternion = Quaternion.fromAngleAxis(angleAxis[0],angleAxis[1:4])
        elif isinstance(Rodrigues, np.ndarray) and Rodrigues.shape == (3,):                # based on given Rodrigues vector
            self.quaternion = Quaternion.fromRodrigues(Rodrigues)
        elif isinstance(quaternion, Quaternion):                                           # based on given quaternion
            self.quaternion = quaternion.homomorph()
        elif isinstance(quaternion, np.ndarray) and quaternion.shape == (4,):              # based on given quaternion
            self.quaternion = Quaternion(quaternion).homomorph()

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

    def asAngleAxis(self):
        return self.quaternion.asAngleAxis()

    def asMatrix(self):
        return self.quaternion.asMatrix()

    def inFZ(self):
        return self.symmetry.inFZ(self.quaternion.asRodrigues())

    def equivalentQuaternions(self):
        return self.symmetry.equivalentQuaternions(self.quaternion)

    def equivalentOrientations(self):
        return map(lambda q: Orientation(quaternion=q,symmetry=self.symmetry.lattice),
                   self.equivalentQuaternions())


    def reduced(self):
        '''Transform orientation to fall into fundamental zone according to symmetry'''
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
                theQ = they * me
                breaker = lowerSymmetry.inDisorientationSST(theQ.asRodrigues()) #\
#                      or lowerSymmetry.inDisorientationSST(theQ.conjugated().asRodrigues())
                if breaker: break
            if breaker: break

        return Orientation(quaternion=theQ,symmetry=self.symmetry.lattice) #, me.conjugated(), they

    def inversePole(self,axis,SST = True):
        '''
        axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)
        '''

        if SST:                                                                                  # pole requested to be within SST
            for i,q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):            # test all symmetric equivalent quaternions
                pole = q.conjugated()*axis                                                           # align crystal direction to axis
                if self.symmetry.inSST(pole): break                                                  # found SST version
        else:
            pole = self.quaternion.conjugated()*axis                                               # align crystal direction to axis

        return pole

    def IPFcolor(self,axis):
        ''' TSL color of inverse pole figure for given axis '''
        color = np.zeros(3,'d')
        for q in self.symmetry.equivalentQuaternions(self.quaternion):
            pole = q.conjugated()*axis                                                             # align crystal direction to axis
            inSST,color = self.symmetry.inSST(pole,color=True)
            if inSST: break

        return color

    @staticmethod
    def getAverageOrientation(cls, orientationList):
        """
        RETURN THE AVERAGE ORIENTATION
        ref: F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman.
             Averaging Quaternions,
             Journal of Guidance, Control, and Dynamics, Vol. 30, No. 4 (2007), pp. 1193-1197.
             doi: 10.2514/1.28949
        sample usage:
             a = Orientation(Eulers=np.radians([10, 10, 0]), symmetry=3)
             b = Orientation(Eulers=np.radians([20, 0, 0]),  symmetry=3)
             avg = Orientation.getAverageOrientation([a,b])
        """

        if not all(isinstance(item, Orientation) for item in orientationList):
            raise TypeError("Only instances of Orientation can be averaged.")

        N = len(orientationList)
        M = orientationList.pop(0).quaternion.asM()
        for o in orientationList:
          M += o.quaternion.asM()
        eig, vec = np.linalg.eig(M/N)

        return Orientation(quaternion = Quaternion(quatArray = vec.T[eig.argmax()]))

    def related(self, relationModel, direction, targetSymmetry = None):

        if relationModel not in ['KS','GT','GTdash','NW','Pitsch','Bain']:  return None
        if int(direction) == 0:  return None

        # KS from S. Morito et al./Journal of Alloys and Compounds 5775 (2013) S587-S592
        # GT from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
        # GT' from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
        # NW from H. Kitahara et al./Materials Characterization 54 (2005) 378-386
        # Pitsch from Y. He et al./Acta Materialia 53 (2005) 1179-1190
        # Bain from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81

        variant  = int(abs(direction))-1
        (me,other)  = (0,1) if direction > 0 else (1,0)

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
                                  [[  1,  1, -1],[  0,  1,  1]],\
                                  [[  1,  1, -1],[  0,  1,  1]],\
                                  [[  1,  1, -1],[  0,  1,  1]],\
                                  [[  1,  1, -1],[  0,  1,  1]],\
                                  [[  1,  1, -1],[  0,  1,  1]],\
                                  [[  1,  1, -1],[  0,  1,  1]]]),
                  'GT': \
                        np.array([[[  1,  1,  1],[  1,  0,  1]],\
                                  [[  1,  1,  1],[  1,  1,  0]],\
                                  [[  1,  1,  1],[  0,  1,  1]],\
                                  [[ -1, -1,  1],[ -1,  0,  1]],\
                                  [[ -1, -1,  1],[ -1, -1,  0]],\
                                  [[ -1, -1,  1],[  0, -1,  1]],\
                                  [[ -1,  1,  1],[ -1,  0,  1]],\
                                  [[ -1,  1,  1],[ -1,  1,  0]],\
                                  [[ -1,  1,  1],[  0,  1,  1]],\
                                  [[  1, -1,  1],[  1,  0,  1]],\
                                  [[  1, -1,  1],[  1, -1,  0]],\
                                  [[  1, -1,  1],[  0, -1,  1]],\
                                  [[  1,  1,  1],[  1,  1,  0]],\
                                  [[  1,  1,  1],[  0,  1,  1]],\
                                  [[  1,  1,  1],[  1,  0,  1]],\
                                  [[ -1, -1,  1],[ -1, -1,  0]],\
                                  [[ -1, -1,  1],[  0, -1,  1]],\
                                  [[ -1, -1,  1],[ -1,  0,  1]],\
                                  [[ -1,  1,  1],[ -1,  1,  0]],\
                                  [[ -1,  1,  1],[  0,  1,  1]],\
                                  [[ -1,  1,  1],[ -1,  0,  1]],\
                                  [[  1, -1,  1],[  1, -1,  0]],\
                                  [[  1, -1,  1],[  0, -1,  1]],\
                                  [[  1, -1,  1],[  1,  0,  1]]]),
                  'GTdash': \
                        np.array([[[  7, 17, 17],[ 12,  5, 17]],\
                                  [[ 17,  7, 17],[ 17, 12,  5]],\
                                  [[ 17, 17,  7],[  5, 17, 12]],\
                                  [[ -7,-17, 17],[-12, -5, 17]],\
                                  [[-17, -7, 17],[-17,-12,  5]],\
                                  [[-17,-17,  7],[ -5,-17, 12]],\
                                  [[  7,-17,-17],[ 12, -5,-17]],\
                                  [[ 17, -7,-17],[ 17,-12, -5]],\
                                  [[ 17,-17, -7],[  5,-17,-12]],\
                                  [[ -7, 17,-17],[-12,  5,-17]],\
                                  [[-17,  7,-17],[-17, 12, -5]],\
                                  [[-17, 17, -7],[ -5, 17,-12]],\
                                  [[  7, 17, 17],[ 12, 17,  5]],\
                                  [[ 17,  7, 17],[  5, 12, 17]],\
                                  [[ 17, 17,  7],[ 17,  5, 12]],\
                                  [[ -7,-17, 17],[-12,-17,  5]],\
                                  [[-17, -7, 17],[ -5,-12, 17]],\
                                  [[-17,-17,  7],[-17, -5, 12]],\
                                  [[  7,-17,-17],[ 12,-17, -5]],\
                                  [[ 17, -7,-17],[ 5, -12,-17]],\
                                  [[ 17,-17,  7],[ 17, -5,-12]],\
                                  [[ -7, 17,-17],[-12, 17, -5]],\
                                  [[-17,  7,-17],[ -5, 12,-17]],\
                                  [[-17, 17, -7],[-17,  5,-12]]]),
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
                                  [[ -1, -1,  1],[  0,  1,  1]],\
                                  [[ -1, -1,  1],[  0,  1,  1]],\
                                  [[ -1, -1,  1],[  0,  1,  1]]]),
                  'Pitsch': \
                        np.array([[[  0,  1,  0],[ -1,  0,  1]],\
                                  [[  0,  0,  1],[  1, -1,  0]],\
                                  [[  1,  0,  0],[  0,  1, -1]],\
                                  [[  1,  0,  0],[  0, -1, -1]],\
                                  [[  0,  1,  0],[ -1,  0, -1]],\
                                  [[  0,  0,  1],[ -1, -1,  0]],\
                                  [[  0,  1,  0],[ -1,  0, -1]],\
                                  [[  0,  0,  1],[ -1, -1,  0]],\
                                  [[  1,  0,  0],[  0, -1, -1]],\
                                  [[  1,  0,  0],[  0, -1,  1]],\
                                  [[  0,  1,  0],[  1,  0, -1]],\
                                  [[  0,  0,  1],[ -1,  1,  0]]]),
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
                                  [[  1,  0, -1],[ -1, -1,  1]],\
                                  [[  1,  0, -1],[ -1,  1, -1]],\
                                  [[ -1, -1,  0],[ -1, -1,  1]],\
                                  [[ -1, -1,  0],[ -1,  1, -1]],\
                                  [[  0,  1,  1],[ -1, -1,  1]],\
                                  [[  0,  1,  1],[ -1,  1, -1]],\
                                  [[  0, -1,  1],[ -1, -1,  1]],\
                                  [[  0, -1,  1],[ -1,  1, -1]],\
                                  [[ -1,  0, -1],[ -1, -1,  1]],\
                                  [[ -1,  0, -1],[ -1,  1, -1]],\
                                  [[  1,  1,  0],[ -1, -1,  1]],\
                                  [[  1,  1,  0],[ -1,  1, -1]],\
                                  [[ -1,  1,  0],[ -1, -1,  1]],\
                                  [[ -1,  1,  0],[ -1,  1, -1]],\
                                  [[  0, -1, -1],[ -1, -1,  1]],\
                                  [[  0, -1, -1],[ -1,  1, -1]],\
                                  [[  1,  0,  1],[ -1, -1,  1]],\
                                  [[  1,  0,  1],[ -1,  1, -1]]]),
                  'GT': \
                        np.array([[[ -5,-12, 17],[-17, -7, 17]],\
                                  [[ 17, -5,-12],[ 17,-17, -7]],\
                                  [[-12, 17, -5],[ -7, 17,-17]],\
                                  [[  5, 12, 17],[ 17,  7, 17]],\
                                  [[-17,  5,-12],[-17, 17, -7]],\
                                  [[ 12,-17, -5],[  7,-17,-17]],\
                                  [[ -5, 12,-17],[-17,  7,-17]],\
                                  [[ 17,  5, 12],[ 17, 17,  7]],\
                                  [[-12,-17,  5],[ -7,-17, 17]],\
                                  [[  5,-12,-17],[ 17, -7,-17]],\
                                  [[-17, -5, 12],[-17,-17,  7]],\
                                  [[ 12, 17,  5],[  7, 17, 17]],\
                                  [[ -5, 17,-12],[-17, 17, -7]],\
                                  [[-12, -5, 17],[ -7,-17, 17]],\
                                  [[ 17,-12, -5],[ 17, -7,-17]],\
                                  [[  5,-17,-12],[ 17,-17, -7]],\
                                  [[ 12,  5, 17],[  7, 17, 17]],\
                                  [[-17, 12, -5],[-17,  7,-17]],\
                                  [[ -5,-17, 12],[-17,-17,  7]],\
                                  [[-12,  5,-17],[ -7, 17,-17]],\
                                  [[ 17, 12,  5],[ 17,  7, 17]],\
                                  [[  5, 17, 12],[ 17, 17,  7]],\
                                  [[ 12, -5,-17],[  7,-17,-17]],\
                                  [[-17,-12,  5],[-17,  7, 17]]]),
                  'GTdash': \
                        np.array([[[  0,  1, -1],[  1,  1, -1]],\
                                  [[ -1,  0,  1],[ -1,  1,  1]],\
                                  [[  1, -1,  0],[  1, -1,  1]],\
                                  [[  0, -1, -1],[ -1, -1, -1]],\
                                  [[  1,  0,  1],[  1, -1,  1]],\
                                  [[  1, -1,  0],[  1, -1, -1]],\
                                  [[  0,  1, -1],[ -1,  1, -1]],\
                                  [[  1,  0,  1],[  1,  1,  1]],\
                                  [[ -1, -1,  0],[ -1, -1,  1]],\
                                  [[  0, -1, -1],[  1, -1, -1]],\
                                  [[ -1,  0,  1],[ -1, -1,  1]],\
                                  [[ -1, -1,  0],[ -1, -1, -1]],\
                                  [[  0, -1,  1],[  1, -1,  1]],\
                                  [[  1,  0, -1],[  1,  1, -1]],\
                                  [[ -1,  1,  0],[ -1,  1,  1]],\
                                  [[  0,  1,  1],[ -1,  1,  1]],\
                                  [[ -1,  0, -1],[ -1, -1, -1]],\
                                  [[ -1,  1,  0],[ -1,  1, -1]],\
                                  [[  0, -1,  1],[ -1, -1,  1]],\
                                  [[ -1,  0, -1],[ -1,  1, -1]],\
                                  [[  1,  1,  0],[  1,  1,  1]],\
                                  [[  0,  1,  1],[  1,  1,  1]],\
                                  [[  1,  0, -1],[  1, -1, -1]],\
                                  [[  1,  1,  0],[  1,  1, -1]]]),
                  'NW': \
                        np.array([[[  2, -1, -1],[  0, -1,  1]],\
                                  [[ -1,  2, -1],[  0, -1,  1]],\
                                  [[ -1, -1,  2],[  0, -1,  1]],\
                                  [[ -2, -1, -1],[  0, -1,  1]],\
                                  [[  1,  2, -1],[  0, -1,  1]],\
                                  [[  1, -1,  2],[  0, -1,  1]],\
                                  [[  2,  1, -1],[  0, -1,  1]],\
                                  [[ -1, -2, -1],[  0, -1,  1]],\
                                  [[ -1,  1,  2],[  0, -1,  1]],\
                                  [[ -1,  2,  1],[  0, -1,  1]],\
                                  [[ -1,  2,  1],[  0, -1,  1]],\
                                  [[ -1, -1, -2],[  0, -1,  1]]]),
                  'Pitsch': \
                        np.array([[[  1,  0,  1],[  1, -1,  1]],\
                                  [[  1,  1,  0],[  1,  1, -1]],\
                                  [[  0,  1,  1],[ -1,  1,  1]],\
                                  [[  0,  1, -1],[ -1,  1, -1]],\
                                  [[ -1,  0,  1],[ -1, -1,  1]],\
                                  [[  1, -1,  0],[  1, -1, -1]],\
                                  [[  1,  0, -1],[  1, -1, -1]],\
                                  [[ -1,  1,  0],[ -1,  1, -1]],\
                                  [[  0, -1,  1],[ -1, -1,  1]],\
                                  [[  0,  1,  1],[ -1,  1,  1]],\
                                  [[  1,  0,  1],[  1, -1,  1]],\
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
        myMatrix  = np.array([myPlane,myNormal,np.cross(myPlane,myNormal)])

        otherPlane   = [float(i) for i in planes[relationModel][variant,other]]                         # map(float, planes[...]) does not work in python 3
        otherPlane  /= np.linalg.norm(otherPlane)
        otherNormal  = [float(i) for i in normals[relationModel][variant,other]]                        # map(float, planes[...]) does not work in python 3
        otherNormal /= np.linalg.norm(otherNormal)
        otherMatrix  = np.array([otherPlane,otherNormal,np.cross(otherPlane,otherNormal)])

        rot=np.dot(otherMatrix.T,myMatrix)

        return Orientation(matrix=np.dot(rot,self.asMatrix()))                                      # no symmetry information ??