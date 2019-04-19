# -*- coding: UTF-8 no BOM -*-

import numpy as np

P = -1                                                                                              # convention (see DOI:10.1088/0965-0393/23/8/083501)

####################################################################################################
class Quaternion:
    u"""
    Quaternion with basic operations
    
    q is the real part, p = (x, y, z) are the imaginary parts.
    Definition of multiplication depends on variable P, P ∈ {-1,1}.
    """                                                                              

    def __init__(self,
                 q    = 0.0,
                 p    = np.zeros(3,dtype=float)):
      """Initializes to identity unless specified"""
      self.q = float(q)
      self.p = np.array(p,dtype=float)


    def __copy__(self):
      """Copy"""
      return self.__class__(q=self.q,
                            p=self.p.copy())

    copy = __copy__


    def __iter__(self):
      """Components"""
      return iter(self.asList())

    def __repr__(self):
      """Readable string"""
      return 'Quaternion: (real={q:+.6f}, imag=<{p[0]:+.6f}, {p[1]:+.6f}, {p[2]:+.6f}>)'.format(q=self.q,p=self.p)
      

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
      else:
        return NotImplemented
      
    def __pos__(self):                                                                      
      """Unary positive operator"""                                                                       
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
      else:
        return NotImplemented
        
    def __neg__(self):                                                                      
      """Unary positive operator"""  
      self.q *= -1.0
      self.p *= -1.0                                                                     
      return self
      
      
    def __mul__(self, other):
      """Multiplication with quaternion or scalar"""
      if isinstance(other, Quaternion):
        return self.__class__(q=self.q*other.q - np.dot(self.p,other.p),
                              p=self.q*other.p + other.q*self.p + P * np.cross(self.p,other.p))
      elif isinstance(other, (int, float)):
        return self.__class__(q=self.q*other,
                              p=self.p*other)
      else:
        return NotImplemented

    def __imul__(self, other):
      """In-place multiplication with quaternion or scalar"""
      if isinstance(other, Quaternion):
        self.q = self.q*other.q - np.dot(self.p,other.p)
        self.p = self.q*other.p + other.q*self.p + P * np.cross(self.p,other.p)
        return self
      elif isinstance(other, (int, float)):
        self.q *= other
        self.p *= other
        return self
      else:
        return NotImplemented


    def __truediv__(self, other):
      """Divsion with quaternion or scalar"""
      if isinstance(other, Quaternion):
        s = other.conjugate()/abs(other)**2.
        return self.__class__(q=self.q * s,
                              p=self.p * s)
      elif isinstance(other, (int, float)):
        return self.__class__(q=self.q / other,
                              p=self.p / other)
      else:
        return NotImplemented

    def __itruediv__(self, other):
      """In-place divsion with quaternion or scalar"""
      if isinstance(other, Quaternion):
        s = other.conjugate()/abs(other)**2.
        self *= s
        return self
      elif isinstance(other, (int, float)):
        self.q /= other
        self.p /= other
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
        return self
      else:
        return NotImplemented
      

    def __abs__(self):
      """Norm"""
      return np.sqrt(self.q ** 2 + np.dot(self.p,self.p))

    magnitude = __abs__


    def __eq__(self,other):
      """Equal (sufficiently close) to each other"""
      return np.isclose(( self-other).magnitude(),0.0) \
          or np.isclose((-self-other).magnitude(),0.0)

    def __ne__(self,other):
      """Not equal (sufficiently close) to each other"""
      return not self.__eq__(other)


    def asM(self):
      """Intermediate representation useful for quaternion averaging (see F. Landis Markley et al.)"""
      return np.outer(self.asArray(),self.asArray())

    def asArray(self):
      """As numpy array"""
      return np.array((self.q,self.p[0],self.p[1],self.p[2]))
 
    def asList(self):
      """As list"""
      return [self.q]+list(self.p)


    def normalize(self):
      """Normalizes in-place"""
      d = self.magnitude()
      if d > 0.0: self /= d
      return self
      
    def normalized(self):
      """Returns normalized copy"""
      return self.copy().normalize()


    def conjugate(self):
      """Conjugates in-place"""
      self.p = -self.p
      return self
      
    def conjugated(self):
      """Returns conjugated copy"""
      return self.copy().conjugate()


    def homomorph(self):
      """Homomorphs in-place"""
      self.q = -self.q
      self.p = -self.p
      return self
      
    def homomorphed(self):
      """Returns homomorphed copy"""
      return self.copy().homomorph()
