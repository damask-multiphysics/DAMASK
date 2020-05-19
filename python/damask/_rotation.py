import numpy as np

from . import mechanics

_P = -1

# parameters for conversion from/to cubochoric
_sc   = np.pi**(1./6.)/6.**(1./6.)
_beta = np.pi**(5./6.)/6.**(1./6.)/2.
_R1   = (3.*np.pi/4.)**(1./3.)


class Rotation:
    u"""
    Orientation stored with functionality for conversion to different representations.

    References
    ----------
    D. Rowenhorst et al., Modelling and Simulation in Materials Science and Engineering 23:083501, 2015
    https://doi.org/10.1088/0965-0393/23/8/083501

    Conventions
    -----------
    Convention 1: Coordinate frames are right-handed.
    Convention 2: A rotation angle ω is taken to be positive for a counterclockwise rotation
                  when viewing from the end point of the rotation axis towards the origin.
    Convention 3: Rotations will be interpreted in the passive sense.
    Convention 4: Euler angle triplets are implemented using the Bunge convention,
                  with the angular ranges as [0, 2π],[0, π],[0, 2π].
    Convention 5: The rotation angle ω is limited to the interval [0, π].
    Convention 6: the real part of a quaternion is positive, Re(q) > 0
    Convention 7: P = -1 (as default).

    Usage
    -----
    Vector "a" (defined in coordinate system "A") is passively rotated
               resulting in new coordinates "b" when expressed in system "B".
    b = Q @ a
    b = np.dot(Q.as_matrix(),a)

    """

    __slots__ = ['quaternion']

    def __init__(self,quaternion = np.array([1.0,0.0,0.0,0.0])):
        """
        Initializes to identity unless specified.

        Parameters
        ----------
        quaternion : numpy.ndarray, optional
            Unit quaternion that follows the conventions. Use .from_quaternion to perform a sanity check.

        """
        self.quaternion = quaternion.copy()


    @property
    def shape(self):
        return self.quaternion.shape[:-1]


    def __copy__(self):
        """Copy."""
        return self.__class__(self.quaternion)

    copy = __copy__


    def __repr__(self):
        """Orientation displayed as unit quaternion, rotation matrix, and Bunge-Euler angles."""
        if self.quaternion.shape != (4,):
            raise NotImplementedError('Support for multiple rotations missing')
        return '\n'.join([
               'Quaternion: (real={:.3f}, imag=<{:+.3f}, {:+.3f}, {:+.3f}>)'.format(*(self.quaternion)),
               'Matrix:\n{}'.format(self.as_matrix()),
               'Bunge Eulers / deg: ({:3.2f}, {:3.2f}, {:3.2f})'.format(*self.as_Eulers(degrees=True)),
                ])


    def __mul__(self, other):
        """
        Multiplication.

        Parameters
        ----------
        other : numpy.ndarray or Rotation
            Vector, second or fourth order tensor, or rotation object that is rotated.

        Todo
        ----
        Document details active/passive)
        consider rotation of (3,3,3,3)-matrix

        """
        if self.quaternion.shape != (4,):
            raise NotImplementedError('Support for multiple rotations missing')
        if isinstance(other, Rotation):
            self_q  = self.quaternion[0]
            self_p  = self.quaternion[1:]
            other_q = other.quaternion[0]
            other_p = other.quaternion[1:]
            R = self.__class__(np.append(self_q*other_q - np.dot(self_p,other_p),
                                         self_q*other_p + other_q*self_p + _P * np.cross(self_p,other_p)))
            return R.standardize()
        elif isinstance(other, np.ndarray):
            if other.shape == (3,):
                A = self.quaternion[0]**2.0 - np.dot(self.quaternion[1:],self.quaternion[1:])
                B = 2.0 * np.dot(self.quaternion[1:],other)
                C = 2.0 * _P*self.quaternion[0]

                return A*other + B*self.quaternion[1:] + C * np.cross(self.quaternion[1:],other)

            elif other.shape == (3,3,):
                R = self.as_matrix()
                return np.dot(R,np.dot(other,R.T))
            elif other.shape == (3,3,3,3,):
                R = self.as_matrix()
                RR = np.outer(R, R)
                RRRR = np.outer(RR, RR).reshape(4 * (3,3))
                axes = ((0, 2, 4, 6), (0, 1, 2, 3))
                return np.tensordot(RRRR, other, axes)
            else:
                raise ValueError('Can only rotate vectors, 2nd order ternsors, and 4th order tensors')
        else:
            raise TypeError('Cannot rotate {}'.format(type(other)))


    def __matmul__(self, other):
        """
        Rotation.

        details to be discussed
        """
        if isinstance(other, Rotation):
            q_m = self.quaternion[...,0:1]
            p_m = self.quaternion[...,1:]
            q_o = other.quaternion[...,0:1]
            p_o = other.quaternion[...,1:]
            q = (q_m*q_o - np.einsum('...i,...i',p_m,p_o).reshape(self.shape+(1,)))
            p = q_m*p_o + q_o*p_m + _P * np.cross(p_m,p_o)
            return self.__class__(np.block([q,p])).standardize()

        elif isinstance(other,np.ndarray):
            if self.shape + (3,) == other.shape:
                q_m = self.quaternion[...,0]
                p_m = self.quaternion[...,1:]
                A = q_m**2.0 - np.einsum('...i,...i',p_m,p_m)
                B = 2.0 * np.einsum('...i,...i',p_m,other)
                C = 2.0 * _P * q_m
                return np.block([(A * other[...,i]).reshape(self.shape+(1,)) +
                                 (B * p_m[...,i]).reshape(self.shape+(1,)) +
                                 (C * (  p_m[...,(i+1)%3]*other[...,(i+2)%3]\
                                       - p_m[...,(i+2)%3]*other[...,(i+1)%3])).reshape(self.shape+(1,))
                                 for i in [0,1,2]])
            if self.shape + (3,3) == other.shape:
                R = self.as_matrix()
                return np.einsum('...im,...jn,...mn',R,R,other)
            if self.shape + (3,3,3,3) == other.shape:
                R = self.as_matrix()
                return np.einsum('...im,...jn,...ko,...lp,...mnop',R,R,R,R,other)
            else:
                raise ValueError('Can only rotate vectors, 2nd order ternsors, and 4th order tensors')
        else:
            raise TypeError('Cannot rotate {}'.format(type(other)))

    def inverse(self):
        """In-place inverse rotation/backward rotation."""
        self.quaternion[...,1:] *= -1
        return self

    def inversed(self):
        """Inverse rotation/backward rotation."""
        return self.copy().inverse()


    def standardize(self):
        """In-place quaternion representation with positive real part."""
        self.quaternion[self.quaternion[...,0] < 0.0] *= -1
        return self

    def standardized(self):
        """Quaternion representation with positive real part."""
        return self.copy().standardize()


    def misorientation(self,other):
        """
        Get Misorientation.

        Parameters
        ----------
        other : Rotation
            Rotation to which the misorientation is computed.

        """
        return other*self.inversed()


    def broadcast_to(self,shape):
        if self.shape == ():
            q = np.broadcast_to(self.quaternion,shape+(4,))
        else:
            q = np.block([np.broadcast_to(self.quaternion[...,0:1],shape+(1,)),
                          np.broadcast_to(self.quaternion[...,1:2],shape+(1,)),
                          np.broadcast_to(self.quaternion[...,2:3],shape+(1,)),
                          np.broadcast_to(self.quaternion[...,3:4],shape+(1,))])
        return self.__class__(q)


    def average(self,other):
        """
        Calculate the average rotation.

        Parameters
        ----------
        other : Rotation
            Rotation from which the average is rotated.

        """
        if self.quaternion.shape != (4,) or other.quaternion.shape != (4,):
            raise NotImplementedError('Support for multiple rotations missing')
        return Rotation.fromAverage([self,other])


    ################################################################################################
    # convert to different orientation representations (numpy arrays)

    def as_quaternion(self):
        """
        Unit quaternion [q, p_1, p_2, p_3].

        Parameters
        ----------
        quaternion : bool, optional
            return quaternion as DAMASK object.

        """
        return self.quaternion

    def as_Eulers(self,
                  degrees = False):
        """
        Bunge-Euler angles: (φ_1, ϕ, φ_2).

        Parameters
        ----------
        degrees : bool, optional
            return angles in degrees.

        """
        eu = Rotation.qu2eu(self.quaternion)
        if degrees: eu = np.degrees(eu)
        return eu

    def as_axis_angle(self,
                      degrees = False,
                      pair = False):
        """
        Axis angle representation [n_1, n_2, n_3, ω] unless pair == True: ([n_1, n_2, n_3], ω).

        Parameters
        ----------
        degrees : bool, optional
            return rotation angle in degrees.
        pair : bool, optional
            return tuple of axis and angle.

        """
        ax = Rotation.qu2ax(self.quaternion)
        if degrees: ax[...,3] = np.degrees(ax[...,3])
        return (ax[...,:3],ax[...,3]) if pair else ax

    def as_matrix(self):
        """Rotation matrix."""
        return Rotation.qu2om(self.quaternion)

    def as_Rodrigues(self,
                     vector = False):
        """
        Rodrigues-Frank vector representation [n_1, n_2, n_3, tan(ω/2)] unless vector == True: [n_1, n_2, n_3] * tan(ω/2).

        Parameters
        ----------
        vector : bool, optional
            return as actual Rodrigues--Frank vector, i.e. rotation axis scaled by tan(ω/2).

        """
        ro = Rotation.qu2ro(self.quaternion)
        return ro[...,:3]*ro[...,3] if vector else ro

    def as_homochoric(self):
        """Homochoric vector: (h_1, h_2, h_3)."""
        return Rotation.qu2ho(self.quaternion)

    def as_cubochoric(self):
        """Cubochoric vector: (c_1, c_2, c_3)."""
        return Rotation.qu2cu(self.quaternion)

    def M(self): # ToDo not sure about the name: as_M or M? we do not have a from_M
        """
        Intermediate representation supporting quaternion averaging.

        References
        ----------
        F. Landis Markley et al., Journal of Guidance, Control, and Dynamics 30(4):1193-1197, 2007
        https://doi.org/10.2514/1.28949

        """
        return np.einsum('...i,...j',self.quaternion,self.quaternion)


    ################################################################################################
    # Static constructors. The input data needs to follow the conventions, options allow to
    # relax the conventions.
    @staticmethod
    def from_quaternion(quaternion,
                        accept_homomorph = False,
                        P = -1):

        qu = np.array(quaternion,dtype=float)
        if qu.shape[:-2:-1] != (4,):
            raise ValueError('Invalid shape.')

        if P > 0: qu[...,1:4] *= -1                                                                 # convert from P=1 to P=-1
        if accept_homomorph:
            qu[qu[...,0] < 0.0] *= -1
        else:
            if np.any(qu[...,0] < 0.0):
                raise ValueError('Quaternion with negative first (real) component.')
        if not np.all(np.isclose(np.linalg.norm(qu,axis=-1), 1.0)):
            raise ValueError('Quaternion is not of unit length.')

        return Rotation(qu)

    @staticmethod
    def from_Eulers(eulers,
                    degrees = False):

        eu = np.array(eulers,dtype=float)
        if eu.shape[:-2:-1] != (3,):
            raise ValueError('Invalid shape.')

        eu = np.radians(eu) if degrees else eu
        if np.any(eu < 0.0) or np.any(eu > 2.0*np.pi) or np.any(eu[...,1] > np.pi):                 # ToDo: No separate check for PHI
            raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π].')

        return Rotation(Rotation.eu2qu(eu))

    @staticmethod
    def from_axis_angle(axis_angle,
                        degrees = False,
                        normalise = False,
                        P = -1):

        ax = np.array(axis_angle,dtype=float)
        if ax.shape[:-2:-1] != (4,):
            raise ValueError('Invalid shape.')

        if P > 0:     ax[...,0:3] *= -1                                                             # convert from P=1 to P=-1
        if degrees:   ax[...,  3]  = np.radians(ax[...,3])
        if normalise: ax[...,0:3] /= np.linalg.norm(ax[...,0:3],axis=-1)
        if np.any(ax[...,3] < 0.0) or np.any(ax[...,3] > np.pi):
            raise ValueError('Axis angle rotation angle outside of [0..π].')
        if not np.all(np.isclose(np.linalg.norm(ax[...,0:3],axis=-1), 1.0)):
            raise ValueError('Axis angle rotation axis is not of unit length.')

        return Rotation(Rotation.ax2qu(ax))

    @staticmethod
    def from_basis(basis,
                   orthonormal = True,
                   reciprocal = False):

        om = np.array(basis,dtype=float)
        if om.shape[:-3:-1] != (3,3):
            raise ValueError('Invalid shape.')

        if reciprocal:
            om = np.linalg.inv(mechanics.transpose(om)/np.pi)                                       # transform reciprocal basis set
            orthonormal = False                                                                     # contains stretch
        if not orthonormal:
            (U,S,Vh) = np.linalg.svd(om)                                                            # singular value decomposition
            om = np.einsum('...ij,...jl->...il',U,Vh)
        if not np.all(np.isclose(np.linalg.det(om),1.0)):
            raise ValueError('Orientation matrix has determinant ≠ 1.')
        if    not np.all(np.isclose(np.einsum('...i,...i',om[...,0],om[...,1]), 0.0)) \
           or not np.all(np.isclose(np.einsum('...i,...i',om[...,1],om[...,2]), 0.0)) \
           or not np.all(np.isclose(np.einsum('...i,...i',om[...,2],om[...,0]), 0.0)):
            raise ValueError('Orientation matrix is not orthogonal.')

        return Rotation(Rotation.om2qu(om))

    @staticmethod
    def from_matrix(om):

        return Rotation.from_basis(om)

    @staticmethod
    def from_Rodrigues(rodrigues,
                       normalise = False,
                       P = -1):

        ro = np.array(rodrigues,dtype=float)
        if ro.shape[:-2:-1] != (4,):
            raise ValueError('Invalid shape.')

        if P > 0:     ro[...,0:3] *= -1                                                             # convert from P=1 to P=-1
        if normalise: ro[...,0:3] /= np.linalg.norm(ro[...,0:3],axis=-1)
        if np.any(ro[...,3] < 0.0):
            raise ValueError('Rodrigues vector rotation angle not positive.')
        if not np.all(np.isclose(np.linalg.norm(ro[...,0:3],axis=-1), 1.0)):
            raise ValueError('Rodrigues vector rotation axis is not of unit length.')

        return Rotation(Rotation.ro2qu(ro))

    @staticmethod
    def from_homochoric(homochoric,
                        P = -1):

        ho = np.array(homochoric,dtype=float)
        if ho.shape[:-2:-1] != (3,):
            raise ValueError('Invalid shape.')

        if P > 0: ho *= -1                                                                          # convert from P=1 to P=-1

        if np.any(np.linalg.norm(ho,axis=-1) > (3.*np.pi/4.)**(1./3.)+1e-9):
            raise ValueError('Homochoric coordinate outside of the sphere.')

        return Rotation(Rotation.ho2qu(ho))

    @staticmethod
    def from_cubochoric(cubochoric,
                       P = -1):

        cu = np.array(cubochoric,dtype=float)
        if cu.shape[:-2:-1] != (3,):
            raise ValueError('Invalid shape.')

        if np.abs(np.max(cu))>np.pi**(2./3.) * 0.5+1e-9:
            raise ValueError('Cubochoric coordinate outside of the cube: {} {} {}.'.format(*cu))

        ho = Rotation.cu2ho(cu)
        if P > 0: ho *= -1                                                                          # convert from P=1 to P=-1

        return Rotation(Rotation.ho2qu(ho))


    @staticmethod
    def fromAverage(rotations,weights = None):
        """
        Average rotation.

        References
        ----------
        F. Landis Markley et al., Journal of Guidance, Control, and Dynamics 30(4):1193-1197, 2007
        https://doi.org/10.2514/1.28949

        Parameters
        ----------
        rotations : list of Rotations
            Rotations to average from
        weights : list of floats, optional
            Weights for each rotation used for averaging

        """
        if not all(isinstance(item, Rotation) for item in rotations):
            raise TypeError('Only instances of Rotation can be averaged.')

        N = len(rotations)
        if not weights:
            weights = np.ones(N,dtype='i')

        for i,(r,n) in enumerate(zip(rotations,weights)):
            M =          r.M() * n if i == 0 \
                else M + r.M() * n                                                                  # noqa add (multiples) of this rotation to average noqa
        eig, vec = np.linalg.eig(M/N)

        return Rotation.from_quaternion(np.real(vec.T[eig.argmax()]),accept_homomorph = True)

    @staticmethod
    def from_random(shape=None):
        if shape is None:
            r = np.random.random(3)
        elif hasattr(shape, '__iter__'):
            r = np.random.random(tuple(shape)+(3,))
        else:
            r = np.random.random((shape,3))

        A = np.sqrt(r[...,2])
        B = np.sqrt(1.0-r[...,2])
        q = np.stack([np.cos(2.0*np.pi*r[...,0])*A,
                      np.sin(2.0*np.pi*r[...,1])*B,
                      np.cos(2.0*np.pi*r[...,1])*B,
                      np.sin(2.0*np.pi*r[...,0])*A],axis=-1)

        return Rotation(q.reshape(r.shape[:-1]+(4,)) if shape is not None else q).standardize()


    # for compatibility (old names do not follow convention)
    asM            = M
    fromQuaternion = from_quaternion
    fromEulers     = from_Eulers
    asAxisAngle    = as_axis_angle
    asRodrigues    = as_Rodrigues


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
    #---------- Quaternion ----------
    @staticmethod
    def qu2om(qu):
        qq = qu[...,0:1]**2-(qu[...,1:2]**2 + qu[...,2:3]**2 + qu[...,3:4]**2)
        om = np.block([qq + 2.0*qu[...,1:2]**2,
                       2.0*(qu[...,2:3]*qu[...,1:2]-_P*qu[...,0:1]*qu[...,3:4]),
                       2.0*(qu[...,3:4]*qu[...,1:2]+_P*qu[...,0:1]*qu[...,2:3]),
                       2.0*(qu[...,1:2]*qu[...,2:3]+_P*qu[...,0:1]*qu[...,3:4]),
                       qq + 2.0*qu[...,2:3]**2,
                       2.0*(qu[...,3:4]*qu[...,2:3]-_P*qu[...,0:1]*qu[...,1:2]),
                       2.0*(qu[...,1:2]*qu[...,3:4]-_P*qu[...,0:1]*qu[...,2:3]),
                       2.0*(qu[...,2:3]*qu[...,3:4]+_P*qu[...,0:1]*qu[...,1:2]),
                       qq + 2.0*qu[...,3:4]**2,
                      ]).reshape(qu.shape[:-1]+(3,3))
        return om

    @staticmethod
    def qu2eu(qu):
        """Quaternion to Bunge-Euler angles."""
        q02   = qu[...,0:1]*qu[...,2:3]
        q13   = qu[...,1:2]*qu[...,3:4]
        q01   = qu[...,0:1]*qu[...,1:2]
        q23   = qu[...,2:3]*qu[...,3:4]
        q03_s = qu[...,0:1]**2+qu[...,3:4]**2
        q12_s = qu[...,1:2]**2+qu[...,2:3]**2
        chi = np.sqrt(q03_s*q12_s)

        eu = np.where(np.abs(q12_s) < 1.0e-8,
                np.block([np.arctan2(-_P*2.0*qu[...,0:1]*qu[...,3:4],qu[...,0:1]**2-qu[...,3:4]**2),
                          np.zeros(qu.shape[:-1]+(2,))]),
                      np.where(np.abs(q03_s) < 1.0e-8,
                          np.block([np.arctan2(   2.0*qu[...,1:2]*qu[...,2:3],qu[...,1:2]**2-qu[...,2:3]**2),
                                    np.broadcast_to(np.pi,qu.shape[:-1]+(1,)),
                                    np.zeros(qu.shape[:-1]+(1,))]),
                          np.block([np.arctan2((-_P*q02+q13)*chi, (-_P*q01-q23)*chi),
                                    np.arctan2( 2.0*chi,          q03_s-q12_s    ),
                                    np.arctan2(( _P*q02+q13)*chi, (-_P*q01+q23)*chi)])
                              )
                     )
        # reduce Euler angles to definition range
        eu[np.abs(eu)<1.e-6] = 0.0
        eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
        return eu

    @staticmethod
    def qu2ax(qu):
        """
        Quaternion to axis angle pair.

        Modified version of the original formulation, should be numerically more stable
        """
        with np.errstate(invalid='ignore',divide='ignore'):
            s = np.sign(qu[...,0:1])/np.sqrt(qu[...,1:2]**2+qu[...,2:3]**2+qu[...,3:4]**2)
            omega = 2.0 * np.arccos(np.clip(qu[...,0:1],-1.0,1.0))
            ax = np.where(np.broadcast_to(qu[...,0:1] < 1.0e-8,qu.shape),
                          np.block([qu[...,1:4],np.broadcast_to(np.pi,qu.shape[:-1]+(1,))]),
                          np.block([qu[...,1:4]*s,omega]))
        ax[np.isclose(qu[...,0],1.,rtol=0.0)] = [0.0, 0.0, 1.0, 0.0]
        return ax

    @staticmethod
    def qu2ro(qu):
        """Quaternion to Rodrigues-Frank vector."""
        with np.errstate(invalid='ignore',divide='ignore'):
            s  = np.linalg.norm(qu[...,1:4],axis=-1,keepdims=True)
            ro = np.where(np.broadcast_to(np.abs(qu[...,0:1]) < 1.0e-12,qu.shape),
                          np.block([qu[...,1:2], qu[...,2:3], qu[...,3:4], np.broadcast_to(np.inf,qu.shape[:-1]+(1,))]),
                          np.block([qu[...,1:2]/s,qu[...,2:3]/s,qu[...,3:4]/s,
                                    np.tan(np.arccos(np.clip(qu[...,0:1],-1.0,1.0)))
                                   ])
                       )
        ro[np.abs(s).squeeze(-1) < 1.0e-12] = [0.0,0.0,_P,0.0]
        return ro

    @staticmethod
    def qu2ho(qu):
        """Quaternion to homochoric vector."""
        with np.errstate(invalid='ignore'):
            omega = 2.0 * np.arccos(np.clip(qu[...,0:1],-1.0,1.0))
            ho = np.where(np.abs(omega) < 1.0e-12,
                          np.zeros(3),
                          qu[...,1:4]/np.linalg.norm(qu[...,1:4],axis=-1,keepdims=True) \
                          * (0.75*(omega - np.sin(omega)))**(1./3.))
        return ho

    @staticmethod
    def qu2cu(qu):
        """Quaternion to cubochoric vector."""
        return Rotation.ho2cu(Rotation.qu2ho(qu))


    #---------- Rotation matrix ----------
    @staticmethod
    def om2qu(om):
        """
        Rotation matrix to quaternion.

        This formulation is from  www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion.
        The original formulation had issues.
        """
        trace = om[...,0,0:1]+om[...,1,1:2]+om[...,2,2:3]

        with np.errstate(invalid='ignore',divide='ignore'):
            s = [
                 0.5 / np.sqrt( 1.0 + trace),
                 2.0 * np.sqrt( 1.0 + om[...,0,0:1] - om[...,1,1:2] - om[...,2,2:3]),
                 2.0 * np.sqrt( 1.0 + om[...,1,1:2] - om[...,2,2:3] - om[...,0,0:1]),
                 2.0 * np.sqrt( 1.0 + om[...,2,2:3] - om[...,0,0:1] - om[...,1,1:2] )
                 ]
            qu= np.where(trace>0,
                         np.block([0.25 / s[0],
                                  (om[...,2,1:2] - om[...,1,2:3] ) * s[0],
                                  (om[...,0,2:3] - om[...,2,0:1] ) * s[0],
                                  (om[...,1,0:1] - om[...,0,1:2] ) * s[0]]),
                         np.where(om[...,0,0:1] > np.maximum(om[...,1,1:2],om[...,2,2:3]),
                                  np.block([(om[...,2,1:2] - om[...,1,2:3]) / s[1],
                                            0.25 * s[1],
                                            (om[...,0,1:2] + om[...,1,0:1]) / s[1],
                                            (om[...,0,2:3] + om[...,2,0:1]) / s[1]]),
                                  np.where(om[...,1,1:2] > om[...,2,2:3],
                                           np.block([(om[...,0,2:3] - om[...,2,0:1]) / s[2],
                                                     (om[...,0,1:2] + om[...,1,0:1]) / s[2],
                                                     0.25 * s[2],
                                                     (om[...,1,2:3] + om[...,2,1:2]) / s[2]]),
                                           np.block([(om[...,1,0:1] - om[...,0,1:2]) / s[3],
                                                     (om[...,0,2:3] + om[...,2,0:1]) / s[3],
                                                     (om[...,1,2:3] + om[...,2,1:2]) / s[3],
                                                     0.25 * s[3]]),
                                          )
                                 )
                        )*np.array([1,_P,_P,_P])
            qu[qu[...,0]<0] *=-1
        return qu

    @staticmethod
    def om2eu(om):
        """Rotation matrix to Bunge-Euler angles."""
        with np.errstate(invalid='ignore',divide='ignore'):
            zeta = 1.0/np.sqrt(1.0-om[...,2,2:3]**2)
            eu = np.where(np.isclose(np.abs(om[...,2,2:3]),1.0,1e-9),
                          np.block([np.arctan2(om[...,0,1:2],om[...,0,0:1]),
                                    np.pi*0.5*(1-om[...,2,2:3]),
                                    np.zeros(om.shape[:-2]+(1,)),
                                   ]),
                          np.block([np.arctan2(om[...,2,0:1]*zeta,-om[...,2,1:2]*zeta),
                                    np.arccos( om[...,2,2:3]),
                                    np.arctan2(om[...,0,2:3]*zeta,+om[...,1,2:3]*zeta)
                                   ])
                          )
        eu[np.abs(eu)<1.e-8] = 0.0
        eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
        return eu

    @staticmethod
    def om2ax(om):
        """Rotation matrix to axis angle pair."""
        return Rotation.qu2ax(Rotation.om2qu(om)) # HOTFIX
        diag_delta = -_P*np.block([om[...,1,2:3]-om[...,2,1:2],
                                   om[...,2,0:1]-om[...,0,2:3],
                                   om[...,0,1:2]-om[...,1,0:1]
                                 ])
        t = 0.5*(om.trace(axis2=-2,axis1=-1) -1.0).reshape(om.shape[:-2]+(1,))
        w,vr = np.linalg.eig(om)
        # mask duplicated real eigenvalues
        w[np.isclose(w[...,0],1.0+0.0j),1:] = 0.
        w[np.isclose(w[...,1],1.0+0.0j),2:] = 0.
        vr = np.swapaxes(vr,-1,-2)
        ax = np.where(np.abs(diag_delta)<1e-12,
                             np.real(vr[np.isclose(w,1.0+0.0j)]).reshape(om.shape[:-2]+(3,)),
                      np.abs(np.real(vr[np.isclose(w,1.0+0.0j)]).reshape(om.shape[:-2]+(3,))) \
                      *np.sign(diag_delta))
        ax = np.block([ax,np.arccos(np.clip(t,-1.0,1.0))])
        ax[np.abs(ax[...,3])<1.e-8] = [ 0.0, 0.0, 1.0, 0.0]
        return ax

    @staticmethod
    def om2ro(om):
        """Rotation matrix to Rodrigues-Frank vector."""
        return Rotation.eu2ro(Rotation.om2eu(om))

    @staticmethod
    def om2ho(om):
        """Rotation matrix to homochoric vector."""
        return Rotation.ax2ho(Rotation.om2ax(om))

    @staticmethod
    def om2cu(om):
        """Rotation matrix to cubochoric vector."""
        return Rotation.ho2cu(Rotation.om2ho(om))


    #---------- Bunge-Euler angles ----------
    @staticmethod
    def eu2qu(eu):
        """Bunge-Euler angles to quaternion."""
        ee = 0.5*eu
        cPhi = np.cos(ee[...,1:2])
        sPhi = np.sin(ee[...,1:2])
        qu = np.block([    cPhi*np.cos(ee[...,0:1]+ee[...,2:3]),
                       -_P*sPhi*np.cos(ee[...,0:1]-ee[...,2:3]),
                       -_P*sPhi*np.sin(ee[...,0:1]-ee[...,2:3]),
                       -_P*cPhi*np.sin(ee[...,0:1]+ee[...,2:3])])
        qu[qu[...,0]<0.0]*=-1
        return qu

    @staticmethod
    def eu2om(eu):
        """Bunge-Euler angles to rotation matrix."""
        c = np.cos(eu)
        s = np.sin(eu)
        om = np.block([+c[...,0:1]*c[...,2:3]-s[...,0:1]*s[...,2:3]*c[...,1:2],
                       +s[...,0:1]*c[...,2:3]+c[...,0:1]*s[...,2:3]*c[...,1:2],
                       +s[...,2:3]*s[...,1:2],
                       -c[...,0:1]*s[...,2:3]-s[...,0:1]*c[...,2:3]*c[...,1:2],
                       -s[...,0:1]*s[...,2:3]+c[...,0:1]*c[...,2:3]*c[...,1:2],
                       +c[...,2:3]*s[...,1:2],
                       +s[...,0:1]*s[...,1:2],
                       -c[...,0:1]*s[...,1:2],
                       +c[...,1:2]
                       ]).reshape(eu.shape[:-1]+(3,3))
        om[np.abs(om)<1.e-12] = 0.0
        return om

    @staticmethod
    def eu2ax(eu):
        """Bunge-Euler angles to axis angle pair."""
        t = np.tan(eu[...,1:2]*0.5)
        sigma = 0.5*(eu[...,0:1]+eu[...,2:3])
        delta = 0.5*(eu[...,0:1]-eu[...,2:3])
        tau   = np.linalg.norm(np.block([t,np.sin(sigma)]),axis=-1,keepdims=True)
        alpha = np.where(np.abs(np.cos(sigma))<1.e-12,np.pi,2.0*np.arctan(tau/np.cos(sigma)))
        with np.errstate(invalid='ignore',divide='ignore'):
            ax = np.where(np.broadcast_to(np.abs(alpha)<1.0e-12,eu.shape[:-1]+(4,)),
                          [0.0,0.0,1.0,0.0],
                          np.block([-_P/tau*t*np.cos(delta),
                                    -_P/tau*t*np.sin(delta),
                                    -_P/tau*  np.sin(sigma),
                                     alpha
                                    ]))
        ax[(alpha<0.0).squeeze()] *=-1
        return ax

    @staticmethod
    def eu2ro(eu):
        """Bunge-Euler angles to Rodrigues-Frank vector."""
        ax = Rotation.eu2ax(eu)
        ro = np.block([ax[...,:3],np.tan(ax[...,3:4]*.5)])
        ro[ax[...,3]>=np.pi,3] = np.inf
        ro[np.abs(ax[...,3])<1.e-16] = [ 0.0, 0.0, _P, 0.0 ]
        return ro

    @staticmethod
    def eu2ho(eu):
        """Bunge-Euler angles to homochoric vector."""
        return Rotation.ax2ho(Rotation.eu2ax(eu))

    @staticmethod
    def eu2cu(eu):
        """Bunge-Euler angles to cubochoric vector."""
        return Rotation.ho2cu(Rotation.eu2ho(eu))


    #---------- Axis angle pair ----------
    @staticmethod
    def ax2qu(ax):
        """Axis angle pair to quaternion."""
        c = np.cos(ax[...,3:4]*.5)
        s = np.sin(ax[...,3:4]*.5)
        qu = np.where(np.abs(ax[...,3:4])<1.e-6,[1.0, 0.0, 0.0, 0.0],np.block([c, ax[...,:3]*s]))
        return qu

    @staticmethod
    def ax2om(ax):
        """Axis angle pair to rotation matrix."""
        c = np.cos(ax[...,3:4])
        s = np.sin(ax[...,3:4])
        omc = 1. -c
        om = np.block([c+omc*ax[...,0:1]**2,
                       omc*ax[...,0:1]*ax[...,1:2] + s*ax[...,2:3],
                       omc*ax[...,0:1]*ax[...,2:3] - s*ax[...,1:2],
                       omc*ax[...,0:1]*ax[...,1:2] - s*ax[...,2:3],
                       c+omc*ax[...,1:2]**2,
                       omc*ax[...,1:2]*ax[...,2:3] + s*ax[...,0:1],
                       omc*ax[...,0:1]*ax[...,2:3] + s*ax[...,1:2],
                       omc*ax[...,1:2]*ax[...,2:3] - s*ax[...,0:1],
                       c+omc*ax[...,2:3]**2]).reshape(ax.shape[:-1]+(3,3))
        return om if _P < 0.0 else np.swapaxes(om,-1,-2)

    @staticmethod
    def ax2eu(ax):
        """Rotation matrix to Bunge Euler angles."""
        return Rotation.om2eu(Rotation.ax2om(ax))

    @staticmethod
    def ax2ro(ax):
        """Axis angle pair to Rodrigues-Frank vector."""
        ro = np.block([ax[...,:3],
                       np.where(np.isclose(ax[...,3:4],np.pi,atol=1.e-15,rtol=.0),
                                np.inf,
                                np.tan(ax[...,3:4]*0.5))
                      ])
        ro[np.abs(ax[...,3])<1.e-6] = [.0,.0,_P,.0]
        return ro

    @staticmethod
    def ax2ho(ax):
        """Axis angle pair to homochoric vector."""
        f = (0.75 * ( ax[...,3:4] - np.sin(ax[...,3:4]) ))**(1.0/3.0)
        ho = ax[...,:3] * f
        return ho

    @staticmethod
    def ax2cu(ax):
        """Axis angle pair to cubochoric vector."""
        return Rotation.ho2cu(Rotation.ax2ho(ax))


    #---------- Rodrigues-Frank vector ----------
    @staticmethod
    def ro2qu(ro):
        """Rodrigues-Frank vector to quaternion."""
        return Rotation.ax2qu(Rotation.ro2ax(ro))

    @staticmethod
    def ro2om(ro):
        """Rodgrigues-Frank vector to rotation matrix."""
        return Rotation.ax2om(Rotation.ro2ax(ro))

    @staticmethod
    def ro2eu(ro):
        """Rodrigues-Frank vector to Bunge-Euler angles."""
        return Rotation.om2eu(Rotation.ro2om(ro))

    @staticmethod
    def ro2ax(ro):
        """Rodrigues-Frank vector to axis angle pair."""
        with np.errstate(invalid='ignore',divide='ignore'):
            ax = np.where(np.isfinite(ro[...,3:4]),
                 np.block([ro[...,0:3]*np.linalg.norm(ro[...,0:3],axis=-1,keepdims=True),2.*np.arctan(ro[...,3:4])]),
                 np.block([ro[...,0:3],np.broadcast_to(np.pi,ro[...,3:4].shape)]))
        ax[np.abs(ro[...,3]) < 1.e-6]  = np.array([ 0.0, 0.0, 1.0, 0.0 ])
        return ax

    @staticmethod
    def ro2ho(ro):
        """Rodrigues-Frank vector to homochoric vector."""
        f = np.where(np.isfinite(ro[...,3:4]),2.0*np.arctan(ro[...,3:4]) -np.sin(2.0*np.arctan(ro[...,3:4])),np.pi)
        ho = np.where(np.broadcast_to(np.sum(ro[...,0:3]**2.0,axis=-1,keepdims=True) < 1.e-6,ro[...,0:3].shape),
                      np.zeros(3), ro[...,0:3]* (0.75*f)**(1.0/3.0))
        return ho

    @staticmethod
    def ro2cu(ro):
        """Rodrigues-Frank vector to cubochoric vector."""
        return Rotation.ho2cu(Rotation.ro2ho(ro))


    #---------- Homochoric vector----------
    @staticmethod
    def ho2qu(ho):
        """Homochoric vector to quaternion."""
        return Rotation.ax2qu(Rotation.ho2ax(ho))

    @staticmethod
    def ho2om(ho):
        """Homochoric vector to rotation matrix."""
        return Rotation.ax2om(Rotation.ho2ax(ho))

    @staticmethod
    def ho2eu(ho):
        """Homochoric vector to Bunge-Euler angles."""
        return Rotation.ax2eu(Rotation.ho2ax(ho))

    @staticmethod
    def ho2ax(ho):
        """Homochoric vector to axis angle pair."""
        tfit = np.array([+1.0000000000018852,      -0.5000000002194847,
                         -0.024999992127593126,    -0.003928701544781374,
                         -0.0008152701535450438,   -0.0002009500426119712,
                         -0.00002397986776071756,  -0.00008202868926605841,
                         +0.00012448715042090092,  -0.0001749114214822577,
                         +0.0001703481934140054,   -0.00012062065004116828,
                         +0.000059719705868660826, -0.00001980756723965647,
                         +0.000003953714684212874, -0.00000036555001439719544])
        hmag_squared = np.sum(ho**2.,axis=-1,keepdims=True)
        hm = hmag_squared.copy()
        s = tfit[0] + tfit[1] * hmag_squared
        for i in range(2,16):
            hm *= hmag_squared
            s  += tfit[i] * hm
        with np.errstate(invalid='ignore'):
            ax = np.where(np.broadcast_to(np.abs(hmag_squared)<1.e-6,ho.shape[:-1]+(4,)),
                          [ 0.0, 0.0, 1.0, 0.0 ],
                          np.block([ho/np.sqrt(hmag_squared),2.0*np.arccos(np.clip(s,-1.0,1.0))]))
        return ax

    @staticmethod
    def ho2ro(ho):
        """Axis angle pair to Rodrigues-Frank vector."""
        return Rotation.ax2ro(Rotation.ho2ax(ho))

    @staticmethod
    def ho2cu(ho):
        """
        Homochoric vector to cubochoric vector.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        rs = np.linalg.norm(ho,axis=-1,keepdims=True)

        xyz3 = np.take_along_axis(ho,Rotation._get_pyramid_order(ho,'forward'),-1)

        with np.errstate(invalid='ignore',divide='ignore'):
            # inverse M_3
            xyz2 = xyz3[...,0:2] * np.sqrt( 2.0*rs/(rs+np.abs(xyz3[...,2:3])) )
            qxy = np.sum(xyz2**2,axis=-1,keepdims=True)

            q2 = qxy + np.max(np.abs(xyz2),axis=-1,keepdims=True)**2
            sq2 = np.sqrt(q2)
            q = (_beta/np.sqrt(2.0)/_R1) * np.sqrt(q2*qxy/(q2-np.max(np.abs(xyz2),axis=-1,keepdims=True)*sq2))
            tt = np.clip((np.min(np.abs(xyz2),axis=-1,keepdims=True)**2\
                +np.max(np.abs(xyz2),axis=-1,keepdims=True)*sq2)/np.sqrt(2.0)/qxy,-1.0,1.0)
            T_inv = np.where(np.abs(xyz2[...,1:2]) <= np.abs(xyz2[...,0:1]),
                                np.block([np.ones_like(tt),np.arccos(tt)/np.pi*12.0]),
                                np.block([np.arccos(tt)/np.pi*12.0,np.ones_like(tt)]))*q
            T_inv[xyz2<0.0] *= -1.0
            T_inv[np.broadcast_to(np.isclose(qxy,0.0,rtol=0.0,atol=1.0e-12),T_inv.shape)] = 0.0
            cu = np.block([T_inv, np.where(xyz3[...,2:3]<0.0,-np.ones_like(xyz3[...,2:3]),np.ones_like(xyz3[...,2:3])) \
                                  * rs/np.sqrt(6.0/np.pi),
                          ])/ _sc

        cu[np.isclose(np.sum(np.abs(ho),axis=-1),0.0,rtol=0.0,atol=1.0e-16)] = 0.0
        cu = np.take_along_axis(cu,Rotation._get_pyramid_order(ho,'backward'),-1)

        return cu

    #---------- Cubochoric ----------
    @staticmethod
    def cu2qu(cu):
        """Cubochoric vector to quaternion."""
        return Rotation.ho2qu(Rotation.cu2ho(cu))

    @staticmethod
    def cu2om(cu):
        """Cubochoric vector to rotation matrix."""
        return Rotation.ho2om(Rotation.cu2ho(cu))

    @staticmethod
    def cu2eu(cu):
        """Cubochoric vector to Bunge-Euler angles."""
        return Rotation.ho2eu(Rotation.cu2ho(cu))

    @staticmethod
    def cu2ax(cu):
        """Cubochoric vector to axis angle pair."""
        return Rotation.ho2ax(Rotation.cu2ho(cu))

    @staticmethod
    def cu2ro(cu):
        """Cubochoric vector to Rodrigues-Frank vector."""
        return Rotation.ho2ro(Rotation.cu2ho(cu))

    @staticmethod
    def cu2ho(cu):
        """
        Cubochoric vector to homochoric vector.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        with np.errstate(invalid='ignore',divide='ignore'):
            # get pyramide and scale by grid parameter ratio
            XYZ = np.take_along_axis(cu,Rotation._get_pyramid_order(cu,'forward'),-1) * _sc
            order = np.abs(XYZ[...,1:2]) <= np.abs(XYZ[...,0:1])
            q = np.pi/12.0 * np.where(order,XYZ[...,1:2],XYZ[...,0:1]) \
                           / np.where(order,XYZ[...,0:1],XYZ[...,1:2])
            c = np.cos(q)
            s = np.sin(q)
            q = _R1*2.0**0.25/_beta/ np.sqrt(np.sqrt(2.0)-c) \
              * np.where(order,XYZ[...,0:1],XYZ[...,1:2])

            T = np.block([ (np.sqrt(2.0)*c - 1.0), np.sqrt(2.0) * s]) * q

            # transform to sphere grid (inverse Lambert)
            c = np.sum(T**2,axis=-1,keepdims=True)
            s = c *         np.pi/24.0 /XYZ[...,2:3]**2
            c = c * np.sqrt(np.pi/24.0)/XYZ[...,2:3]
            q = np.sqrt( 1.0 - s)

            ho = np.where(np.isclose(np.sum(np.abs(XYZ[...,0:2]),axis=-1,keepdims=True),0.0,rtol=0.0,atol=1.0e-16),
                          np.block([np.zeros_like(XYZ[...,0:2]),np.sqrt(6.0/np.pi) *XYZ[...,2:3]]),
                          np.block([np.where(order,T[...,0:1],T[...,1:2])*q,
                                    np.where(order,T[...,1:2],T[...,0:1])*q,
                                    np.sqrt(6.0/np.pi) * XYZ[...,2:3] - c])
                          )

        ho[np.isclose(np.sum(np.abs(cu),axis=-1),0.0,rtol=0.0,atol=1.0e-16)] = 0.0
        ho = np.take_along_axis(ho,Rotation._get_pyramid_order(cu,'backward'),-1)

        return ho


    @staticmethod
    def _get_pyramid_order(xyz,direction=None):
        """
        Get order of the coordinates.

        Depending on the pyramid in which the point is located, the order need to be adjusted.

        Parameters
        ----------
        xyz : numpy.ndarray
           coordinates of a point on a uniform refinable grid on a ball or
           in a uniform refinable cubical grid.

        References
        ----------
        D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
        https://doi.org/10.1088/0965-0393/22/7/075013

        """
        order = {'forward': np.array([[0,1,2],[1,2,0],[2,0,1]]),
                 'backward':np.array([[0,1,2],[2,0,1],[1,2,0]])}

        p = np.where(np.maximum(np.abs(xyz[...,0]),np.abs(xyz[...,1])) <= np.abs(xyz[...,2]),0,
                     np.where(np.maximum(np.abs(xyz[...,1]),np.abs(xyz[...,2])) <= np.abs(xyz[...,0]),1,2))

        return order[direction][p]
