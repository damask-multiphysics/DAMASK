import numpy as np

from ._Lambert import ball_to_cube, cube_to_ball

P = -1

def iszero(a):
    return np.isclose(a,0.0,atol=1.0e-12,rtol=0.0)


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
    b = Q * a
    b = np.dot(Q.asMatrix(),a)

    """

    __slots__ = ['quaternion']

    def __init__(self,quaternion = np.array([1.0,0.0,0.0,0.0])):
        """
        Initializes to identity unless specified.

        Parameters
        ----------
        quaternion : numpy.ndarray, optional
            Unit quaternion that follows the conventions. Use .fromQuaternion to perform a sanity check.

        """
        self.quaternion = quaternion.copy()

    def __copy__(self):
        """Copy."""
        return self.__class__(self.quaternion)

    copy = __copy__


    def __repr__(self):
        """Orientation displayed as unit quaternion, rotation matrix, and Bunge-Euler angles."""
        return '\n'.join([
               'Quaternion: (real={:.3f}, imag=<{:+.3f}, {:+.3f}, {:+.3f}>)'.format(*(self.quaternion)),
               'Matrix:\n{}'.format(self.asMatrix()),
               'Bunge Eulers / deg: ({:3.2f}, {:3.2f}, {:3.2f})'.format(*self.asEulers(degrees=True)),
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
        considere rotation of (3,3,3,3)-matrix

        """
        if isinstance(other, Rotation):                                                             # rotate a rotation
            self_q  = self.quaternion[0]
            self_p  = self.quaternion[1:]
            other_q = other.quaternion[0]
            other_p = other.quaternion[1:]
            R = self.__class__(np.append(self_q*other_q - np.dot(self_p,other_p),
                                         self_q*other_p + other_q*self_p + P * np.cross(self_p,other_p)))
            return R.standardize()
        elif isinstance(other, (tuple,np.ndarray)):
            if isinstance(other,tuple) or other.shape == (3,):                                      # rotate a single (3)-vector or meshgrid
                A = self.quaternion[0]**2.0 - np.dot(self.quaternion[1:],self.quaternion[1:])
                B = 2.0 * (  self.quaternion[1]*other[0]
                           + self.quaternion[2]*other[1]
                           + self.quaternion[3]*other[2])
                C = 2.0 * P*self.quaternion[0]

                return np.array([
                  A*other[0] + B*self.quaternion[1] + C*(self.quaternion[2]*other[2] - self.quaternion[3]*other[1]),
                  A*other[1] + B*self.quaternion[2] + C*(self.quaternion[3]*other[0] - self.quaternion[1]*other[2]),
                  A*other[2] + B*self.quaternion[3] + C*(self.quaternion[1]*other[1] - self.quaternion[2]*other[0]),
                  ])
            elif other.shape == (3,3,):                                                             # rotate a single (3x3)-matrix
                return np.dot(self.asMatrix(),np.dot(other,self.asMatrix().T))
            elif other.shape == (3,3,3,3,):
                raise NotImplementedError
            else:
                return NotImplemented
        else:
            return NotImplemented


    def inverse(self):
        """In-place inverse rotation/backward rotation."""
        self.quaternion[1:] *= -1
        return self

    def inversed(self):
        """Inverse rotation/backward rotation."""
        return self.copy().inverse()


    def standardize(self):
        """In-place quaternion representation with positive q."""
        if self.quaternion[0] < 0.0: self.quaternion*=-1
        return self

    def standardized(self):
        """Quaternion representation with positive q."""
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


    def average(self,other):
        """
        Calculate the average rotation.

        Parameters
        ----------
        other : Rotation
            Rotation from which the average is rotated.

        """
        return Rotation.fromAverage([self,other])


    ################################################################################################
    # convert to different orientation representations (numpy arrays)

    def asQuaternion(self):
        """
        Unit quaternion [q, p_1, p_2, p_3] unless quaternion == True: damask.quaternion object.

        Parameters
        ----------
        quaternion : bool, optional
            return quaternion as DAMASK object.

        """
        return self.quaternion

    def asEulers(self,
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

    def asAxisAngle(self,
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
        if degrees: ax[3] = np.degrees(ax[3])
        return (ax[:3],np.degrees(ax[3])) if pair else ax

    def asMatrix(self):
        """Rotation matrix."""
        return Rotation.qu2om(self.quaternion)

    def asRodrigues(self,
                    vector = False):
        """
        Rodrigues-Frank vector representation [n_1, n_2, n_3, tan(ω/2)] unless vector == True: [n_1, n_2, n_3] * tan(ω/2).

        Parameters
        ----------
        vector : bool, optional
            return as actual Rodrigues--Frank vector, i.e. rotation axis scaled by tan(ω/2).

        """
        ro = Rotation.qu2ro(self.quaternion)
        return ro[:3]*ro[3] if vector else ro

    def asHomochoric(self):
        """Homochoric vector: (h_1, h_2, h_3)."""
        return Rotation.qu2ho(self.quaternion)

    def asCubochoric(self):
        """Cubochoric vector: (c_1, c_2, c_3)."""
        return Rotation.qu2cu(self.quaternion)

    def asM(self):
        """
        Intermediate representation supporting quaternion averaging.

        References
        ----------
        F. Landis Markley et al., Journal of Guidance, Control, and Dynamics 30(4):1193-1197, 2007
        https://doi.org/10.2514/1.28949

        """
        return np.outer(self.quaternion,self.quaternion)


    ################################################################################################
    # static constructors. The input data needs to follow the convention, options allow to
    # relax these convections
    @staticmethod
    def fromQuaternion(quaternion,
                       acceptHomomorph = False,
                       P = -1):

        qu =   quaternion if isinstance(quaternion,np.ndarray) and quaternion.dtype == np.dtype(float) \
                          else np.array(quaternion,dtype=float)
        if P > 0: qu[1:4] *= -1                                                                     # convert from P=1 to P=-1
        if qu[0] < 0.0:
            if acceptHomomorph:
                qu *= -1.
            else:
                raise ValueError('Quaternion has negative first component: {}.'.format(qu[0]))
        if not np.isclose(np.linalg.norm(qu), 1.0):
            raise ValueError('Quaternion is not of unit length: {} {} {} {}.'.format(*qu))

        return Rotation(qu)

    @staticmethod
    def fromEulers(eulers,
                   degrees = False):

        eu = eulers if isinstance(eulers, np.ndarray) and eulers.dtype == np.dtype(float) \
                    else np.array(eulers,dtype=float)
        eu = np.radians(eu) if degrees else eu
        if np.any(eu < 0.0) or np.any(eu > 2.0*np.pi) or eu[1] > np.pi:
            raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π]: {} {} {}.'.format(*eu))

        return Rotation(Rotation.eu2qu(eu))

    @staticmethod
    def fromAxisAngle(angleAxis,
                      degrees = False,
                      normalise = False,
                      P = -1):

        ax = angleAxis if isinstance(angleAxis, np.ndarray) and angleAxis.dtype == np.dtype(float) \
                       else np.array(angleAxis,dtype=float)
        if P > 0:     ax[0:3] *= -1                                                                 # convert from P=1 to P=-1
        if degrees:   ax[  3]  = np.radians(ax[3])
        if normalise: ax[0:3] /= np.linalg.norm(ax[0:3])
        if ax[3] < 0.0 or ax[3] > np.pi:
            raise ValueError('Axis angle rotation angle outside of [0..π]: {}.'.format(ax[3]))
        if not np.isclose(np.linalg.norm(ax[0:3]), 1.0):
            raise ValueError('Axis angle rotation axis is not of unit length: {} {} {}.'.format(*ax[0:3]))

        return Rotation(Rotation.ax2qu(ax))

    @staticmethod
    def fromBasis(basis,
                  orthonormal = True,
                  reciprocal = False,
                 ):

        om = basis if isinstance(basis, np.ndarray) else np.array(basis).reshape(3,3)
        if reciprocal:
            om = np.linalg.inv(om.T/np.pi)                                                          # transform reciprocal basis set
            orthonormal = False                                                                     # contains stretch
        if not orthonormal:
            (U,S,Vh) = np.linalg.svd(om)                                                            # singular value decomposition
            om = np.dot(U,Vh)
        if not np.isclose(np.linalg.det(om),1.0):
            raise ValueError('matrix is not a proper rotation: {}.'.format(om))
        if    not np.isclose(np.dot(om[0],om[1]), 0.0) \
           or not np.isclose(np.dot(om[1],om[2]), 0.0) \
           or not np.isclose(np.dot(om[2],om[0]), 0.0):
            raise ValueError('matrix is not orthogonal: {}.'.format(om))

        return Rotation(Rotation.om2qu(om))

    @staticmethod
    def fromMatrix(om,
                  ):

        return Rotation.fromBasis(om)

    @staticmethod
    def fromRodrigues(rodrigues,
                      normalise = False,
                      P = -1):

        ro = rodrigues if isinstance(rodrigues, np.ndarray) and rodrigues.dtype == np.dtype(float) \
                       else np.array(rodrigues,dtype=float)
        if P > 0:     ro[0:3] *= -1                                                                 # convert from P=1 to P=-1
        if normalise: ro[0:3] /= np.linalg.norm(ro[0:3])
        if not np.isclose(np.linalg.norm(ro[0:3]), 1.0):
            raise ValueError('Rodrigues rotation axis is not of unit length: {} {} {}.'.format(*ro[0:3]))
        if ro[3] < 0.0:
            raise ValueError('Rodrigues rotation angle not positive: {}.'.format(ro[3]))

        return Rotation(Rotation.ro2qu(ro))

    @staticmethod
    def fromHomochoric(homochoric,
                       P = -1):

        ho = homochoric if isinstance(homochoric, np.ndarray) and homochoric.dtype == np.dtype(float) \
                        else np.array(homochoric,dtype=float)
        if P > 0: ho *= -1                                                                          # convert from P=1 to P=-1

        return Rotation(Rotation.ho2qu(ho))

    @staticmethod
    def fromCubochoric(cubochoric,
                       P = -1):

        cu = cubochoric if isinstance(cubochoric, np.ndarray) and cubochoric.dtype == np.dtype(float) \
                        else np.array(cubochoric,dtype=float)
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
            M =          r.asM() * n if i == 0 \
                else M + r.asM() * n                                                                # noqa add (multiples) of this rotation to average noqa
        eig, vec = np.linalg.eig(M/N)

        return Rotation.fromQuaternion(np.real(vec.T[eig.argmax()]),acceptHomomorph = True)


    @staticmethod
    def fromRandom():
        r = np.random.random(3)
        A = np.sqrt(r[2])
        B = np.sqrt(1.0-r[2])
        return Rotation(np.array([np.cos(2.0*np.pi*r[0])*A,
                                  np.sin(2.0*np.pi*r[1])*B,
                                  np.cos(2.0*np.pi*r[1])*B,
                                  np.sin(2.0*np.pi*r[0])*A])).standardize()


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
        if len(qu.shape) == 1:
            """Quaternion to rotation matrix."""
            qq = qu[0]**2-(qu[1]**2 + qu[2]**2 + qu[3]**2)
            om = np.diag(qq + 2.0*np.array([qu[1],qu[2],qu[3]])**2)

            om[1,0] = 2.0*(qu[2]*qu[1]+qu[0]*qu[3])
            om[0,1] = 2.0*(qu[1]*qu[2]-qu[0]*qu[3])
            om[2,1] = 2.0*(qu[3]*qu[2]+qu[0]*qu[1])
            om[1,2] = 2.0*(qu[2]*qu[3]-qu[0]*qu[1])
            om[0,2] = 2.0*(qu[1]*qu[3]+qu[0]*qu[2])
            om[2,0] = 2.0*(qu[3]*qu[1]-qu[0]*qu[2])
            return om if P > 0.0 else om.T
        else:
            qq = qu[...,0:1]**2-(qu[...,1:2]**2 + qu[...,2:3]**2 + qu[...,3:4]**2)
            om = np.block([qq + 2.0*qu[...,1:2]**2,
                           2.0*(qu[...,2:3]*qu[...,1:2]+qu[...,0:1]*qu[...,3:4]),
                           2.0*(qu[...,3:4]*qu[...,1:2]-qu[...,0:1]*qu[...,2:3]),
                           2.0*(qu[...,1:2]*qu[...,2:3]-qu[...,0:1]*qu[...,3:4]),
                           qq + 2.0*qu[...,2:3]**2,
                           2.0*(qu[...,3:4]*qu[...,2:3]+qu[...,0:1]*qu[...,1:2]),
                           2.0*(qu[...,1:2]*qu[...,3:4]+qu[...,0:1]*qu[...,2:3]),
                           2.0*(qu[...,2:3]*qu[...,3:4]-qu[...,0:1]*qu[...,1:2]),
                           qq + 2.0*qu[...,3:4]**2,
                          ]).reshape(qu.shape[:-1]+(3,3))
            return om # TODO: TRANSPOSE FOR P = 1

    @staticmethod
    def qu2eu(qu):
        """Quaternion to Bunge-Euler angles."""
        if len(qu.shape) == 1:
            q03 = qu[0]**2+qu[3]**2
            q12 = qu[1]**2+qu[2]**2
            chi = np.sqrt(q03*q12)
            if   np.abs(q03)< 1.e-6:
                eu = np.array([np.arctan2(-P*2.0*qu[0]*qu[3],qu[0]**2-qu[3]**2), 0.0,   0.0])
            elif np.abs(q12)< 1.e-6:
                eu = np.array([np.arctan2(   2.0*qu[1]*qu[2],qu[1]**2-qu[2]**2), np.pi, 0.0])
            else:
                eu = np.array([np.arctan2((-P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-P*qu[0]*qu[1]-qu[2]*qu[3])*chi ),
                               np.arctan2( 2.0*chi, q03-q12 ),
                               np.arctan2(( P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-P*qu[0]*qu[1]+qu[2]*qu[3])*chi )])
        else:
            q02   = qu[...,0:1]*qu[...,2:3]
            q13   = qu[...,1:2]*qu[...,3:4]
            q01   = qu[...,0:1]*qu[...,1:2]
            q23   = qu[...,2:3]*qu[...,3:4]
            q03_s = qu[...,0:1]**2+qu[...,3:4]**2
            q12_s = qu[...,1:2]**2+qu[...,2:3]**2
            chi = np.sqrt(q03_s*q12_s)

            eu = np.where(np.abs(q12_s) < 1.0e-6,
                          np.block([np.arctan2(-P*2.0*qu[...,0:1]*qu[...,3:4],qu[...,0:1]**2-qu[...,3:4]**2),
                                    np.zeros(qu.shape[:-1]+(2,))]),
                          np.block([np.arctan2((-P*q02+q13)*chi, (-P*q01-q23)*chi),
                                    np.arctan2( 2.0*chi,          q03_s-q12_s    ),
                                    np.arctan2(( P*q02+q13)*chi, (-P*q01+q23)*chi)])
                         )
            eu = np.where(np.abs(q03_s) < 1.0e-6,
                          np.block([np.arctan2(   2.0*qu[...,1:2]*qu[...,2:3],qu[...,1:2]**2-qu[...,2:3]**2),
                                    np.ones( qu.shape[:-1]+(1,))*np.pi,
                                    np.zeros(qu.shape[:-1]+(1,))]),
                          eu)                                                                       # TODO: Where not needed
        # reduce Euler angles to definition range, i.e a lower limit of 0.0
        eu[np.abs(eu)<1.e-6] = 0.0
        eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
        return eu

    @staticmethod
    def qu2ax(qu):
        """
        Quaternion to axis angle pair.

        Modified version of the original formulation, should be numerically more stable
        """
        if len(qu.shape) == 1:
            if iszero(np.sum(qu[1:4]**2)):                                                          # set axis to [001] if the angle is 0/360
                ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
            elif np.abs(qu[0]) > 1.e-6:
                s = np.sign(qu[0])/np.sqrt(qu[1]**2+qu[2]**2+qu[3]**2)
                omega = 2.0 * np.arccos(np.clip(qu[0],-1.0,1.0))
                ax = ax = np.array([ qu[1]*s, qu[2]*s, qu[3]*s, omega ])
            else:
                ax = ax = np.array([ qu[1], qu[2], qu[3], np.pi])
        else:
            with np.errstate(divide='ignore'):
                s = np.sign(qu[...,0:1])/np.sqrt(qu[...,1:2]**2+qu[...,2:3]**2+qu[...,3:4]**2)
            omega = 2.0 * np.arccos(np.clip(qu[...,0:1],-1.0,1.0))

            ax = np.where(qu[...,0:1] < 1.0e-6,
                          np.block([qu[...,1:4],np.ones(qu.shape[:-1]+(1,))*np.pi]),
                          np.block([qu[...,1:4]*s,omega]))
            ax = np.where(np.expand_dims(np.sum(np.abs(qu[:,1:4])**2,axis=-1) < 1.0e-6,-1),
                          [0.0, 0.0, 1.0, 0.0], ax)                                                 # TODO: Where not needed
        return ax


    @staticmethod
    def qu2ro(qu):
        """Quaternion to Rodrigues-Frank vector."""
        if len(qu.shape) == 1:
            if iszero(qu[0]):
                ro = np.array([qu[1], qu[2], qu[3], np.inf])
            else:
                s = np.linalg.norm([qu[1],qu[2],qu[3]])
                ro = np.array([0.0,0.0,P,0.0] if iszero(s) else \
                              [ qu[1]/s,  qu[2]/s,  qu[3]/s, np.tan(np.arccos(np.clip(qu[0],-1.0,1.0)))])
        else:
            s  = np.expand_dims(np.linalg.norm(qu[...,1:4],axis=1),-1)
            ro = np.where(np.abs(s) < 1.0e-12,
                          [0.0,0.0,P,0.0],
                          np.block([qu[...,1:2]/s,qu[...,2:3]/s,qu[...,3:4]/s,
                                    np.tan(np.arccos(np.clip(qu[:,0:1],-1.0,1.0)))
                                   ])
                           )
            ro = np.where(np.abs(qu[...,0:1]) < 1.0e-12,
                          np.block([qu[...,1:2], qu[...,2:3], qu[...,3:4], np.ones(qu.shape[:-1]+(1,))*np.inf]),ro) # TODO: Where not needed
        return ro

    @staticmethod
    def qu2ho(qu):
        """Quaternion to homochoric vector."""
        if len(qu.shape) == 1:
            if np.isclose(qu[0],1.0):
                ho = np.zeros(3)
            else:
                omega = 2.0 * np.arccos(np.clip(qu[0],-1.0,1.0))
                ho = np.array([qu[1], qu[2], qu[3]])
                f  = 0.75 * ( omega - np.sin(omega) )
                ho = ho/np.linalg.norm(ho) * f**(1./3.)
        else:
            omega = 2.0 * np.arccos(np.clip(qu[...,0:1],-1.0,1.0))
            ho = np.where(np.abs(omega) < 1.0e-12,
                          np.zeros(3),
                          qu[...,1:4]/np.linalg.norm(qu[...,1:4],axis=1).reshape(qu.shape[:-1]+(1,)) * (0.75*(omega - np.sin(omega)))**(1./3.))
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

        The original formulation (direct conversion) had (numerical?) issues
        """
        return Rotation.eu2qu(Rotation.om2eu(om))

    @staticmethod
    def om2eu(om):
        """Rotation matrix to Bunge-Euler angles."""
        if not np.isclose(np.abs(om[2,2]),1.0,1.e-4):
            zeta = 1.0/np.sqrt(1.0-om[2,2]**2)
            eu = np.array([np.arctan2(om[2,0]*zeta,-om[2,1]*zeta),
                           np.arccos(om[2,2]),
                           np.arctan2(om[0,2]*zeta, om[1,2]*zeta)])
        else:
            eu = np.array([np.arctan2( om[0,1],om[0,0]), np.pi*0.5*(1-om[2,2]),0.0])                # following the paper, not the reference implementation

        # reduce Euler angles to definition range, i.e a lower limit of 0.0
        eu[np.abs(eu)<1.e-6] = 0.0
        eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
        return eu

    @staticmethod
    def om2ax(om):
        """Rotation matrix to axis angle pair."""
        ax=np.empty(4)

        # first get the rotation angle
        t = 0.5*(om.trace() -1.0)
        ax[3] = np.arccos(np.clip(t,-1.0,1.0))

        if np.abs(ax[3])<1.e-6:
            ax = np.array([ 0.0, 0.0, 1.0, 0.0])
        else:
            w,vr = np.linalg.eig(om)
            # next, find the eigenvalue (1,0j)
            i = np.where(np.isclose(w,1.0+0.0j))[0][0]
            ax[0:3] = np.real(vr[0:3,i])
            diagDelta = np.array([om[1,2]-om[2,1],om[2,0]-om[0,2],om[0,1]-om[1,0]])
            ax[0:3] = np.where(np.abs(diagDelta)<1.e-6, ax[0:3],np.abs(ax[0:3])*np.sign(-P*diagDelta))
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
        if len(eu.shape) == 1:
            ee = 0.5*eu
            cPhi = np.cos(ee[1])
            sPhi = np.sin(ee[1])
            qu = np.array([     cPhi*np.cos(ee[0]+ee[2]),
                             -P*sPhi*np.cos(ee[0]-ee[2]),
                             -P*sPhi*np.sin(ee[0]-ee[2]),
                             -P*cPhi*np.sin(ee[0]+ee[2]) ])
            if qu[0] < 0.0: qu*=-1
        else:
            ee = 0.5*eu
            cPhi = np.cos(ee[...,1:2])
            sPhi = np.sin(ee[...,1:2])
            qu = np.block([  cPhi*np.cos(ee[...,0:1]+ee[...,2:3]),
                          -P*sPhi*np.cos(ee[...,0:1]-ee[...,2:3]),
                          -P*sPhi*np.sin(ee[...,0:1]-ee[...,2:3]),
                          -P*cPhi*np.sin(ee[...,0:1]+ee[...,2:3])])
            qu[qu[...,0]<0.0]*=-1
        return qu


    @staticmethod
    def eu2om(eu):
        """Bunge-Euler angles to rotation matrix."""
        if len(eu.shape) == 1:
            c = np.cos(eu)
            s = np.sin(eu)

            om = np.array([[+c[0]*c[2]-s[0]*s[2]*c[1], +s[0]*c[2]+c[0]*s[2]*c[1], +s[2]*s[1]],
                           [-c[0]*s[2]-s[0]*c[2]*c[1], -s[0]*s[2]+c[0]*c[2]*c[1], +c[2]*s[1]],
                           [+s[0]*s[1],                -c[0]*s[1],                +c[1]     ]])
        else:
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
        if len(eu.shape) == 1:
            t = np.tan(eu[1]*0.5)
            sigma = 0.5*(eu[0]+eu[2])
            delta = 0.5*(eu[0]-eu[2])
            tau   = np.linalg.norm([t,np.sin(sigma)])
            alpha = np.pi if iszero(np.cos(sigma)) else \
                    2.0*np.arctan(tau/np.cos(sigma))

            if np.abs(alpha)<1.e-6:
                ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
            else:
                ax = -P/tau * np.array([ t*np.cos(delta), t*np.sin(delta), np.sin(sigma) ])         # passive axis angle pair so a minus sign in front
                ax = np.append(ax,alpha)
                if alpha < 0.0: ax *= -1.0                                                          # ensure alpha is positive
        else:
            t = np.tan(eu[...,1:2]*0.5)
            sigma = 0.5*(eu[...,0:1]+eu[...,2:3])
            delta = 0.5*(eu[...,0:1]-eu[...,2:3])
            tau   = np.linalg.norm(np.block([t,np.sin(sigma)]),axis=-1).reshape(-1,1)
            alpha = np.where(np.abs(np.cos(sigma))<1.e-12,np.pi,2.0*np.arctan(tau/np.cos(sigma)))
            ax = np.where(np.broadcast_to(np.abs(alpha)<1.0e-12,eu.shape[:-1]+(4,)),
                          [0.0,0.0,1.0,0.0],
                          np.block([-P/tau*t*np.cos(delta),
                                    -P/tau*t*np.sin(delta),
                                    -P/tau*  np.sin(sigma),
                                     alpha
                                    ]))
            ax[(alpha<0.0).squeeze()] *=-1
        return ax

    @staticmethod
    def eu2ro(eu):
        """Bunge-Euler angles to Rodrigues-Frank vector."""
        if len(eu.shape) == 1:
            ro = Rotation.eu2ax(eu)                                                                     # convert to axis angle pair representation
            if ro[3] >= np.pi:                                                                          # Differs from original implementation. check convention 5
                ro[3] = np.inf
            elif iszero(ro[3]):
                ro = np.array([ 0.0, 0.0, P, 0.0 ])
            else:
                ro[3] = np.tan(ro[3]*0.5)
        else:
            ax = Rotation.eu2ax(eu)
            ro = np.block([ax[:,:3],np.tan(ax[:,3:4]*.5)])
            ro[ax[:,3]>=np.pi,3] = np.inf
            ro[np.abs(ax[:,3])<1.e-16] = [ 0.0, 0.0, P, 0.0 ]
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
        if len(ax.shape) == 1:
            if np.abs(ax[3])<1.e-6:
                qu = np.array([ 1.0, 0.0, 0.0, 0.0 ])
            else:
                c = np.cos(ax[3]*0.5)
                s = np.sin(ax[3]*0.5)
                qu = np.array([ c, ax[0]*s, ax[1]*s, ax[2]*s ])
            return qu
        else:
            c = np.cos(ax[...,3:4]*.5)
            s = np.sin(ax[...,3:4]*.5)
            qu = np.where(np.abs(ax[...,3:4])<1.e-12,[1.0, 0.0, 0.0, 0.0],np.block([c, ax[...,:3]*s]))
            return qu

    @staticmethod
    def ax2om(ax):
        """Axis angle pair to rotation matrix."""
        if len(ax.shape) == 1:
            c = np.cos(ax[3])
            s = np.sin(ax[3])
            omc = 1.0-c
            om=np.diag(ax[0:3]**2*omc + c)

            for idx in [[0,1,2],[1,2,0],[2,0,1]]:
                q = omc*ax[idx[0]] * ax[idx[1]]
                om[idx[0],idx[1]] = q + s*ax[idx[2]]
                om[idx[1],idx[0]] = q - s*ax[idx[2]]
            return om if P < 0.0 else om.T
        else:
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
            return om # TODO: TRANSPOSE FOR P = 1

    @staticmethod
    def ax2eu(ax):
        """Rotation matrix to Bunge Euler angles."""
        return Rotation.om2eu(Rotation.ax2om(ax))

    @staticmethod
    def ax2ro(ax):
        """Axis angle pair to Rodrigues-Frank vector."""
        if len(ax.shape) == 1:
            if np.abs(ax[3])<1.e-6:
                ro = [ 0.0, 0.0, P, 0.0 ]
            else:
                ro = [ax[0], ax[1], ax[2]]
                # 180 degree case
                ro += [np.inf] if np.isclose(ax[3],np.pi,atol=1.0e-15,rtol=0.0) else \
                      [np.tan(ax[3]*0.5)]
            return np.array(ro)
        else:
            ro = np.block([ax[...,:3],
                           np.where(np.isclose(ax[...,3:4],np.pi,atol=1.e-15,rtol=.0),
                                    np.inf,
                                    np.tan(ax[...,3:4]*0.5))
                          ])
            ro[np.abs(ax[...,3])<1.e-6] = [.0,.0,P,.0]
            return ro

    @staticmethod
    def ax2ho(ax):
        """Axis angle pair to homochoric vector."""
        if len(ax.shape) == 1:
            f = (0.75 * ( ax[3] - np.sin(ax[3]) ))**(1.0/3.0)
            ho = ax[0:3] * f
            return ho
        else:
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

    @staticmethod
    def ro2ho(ro):
        """Rodrigues-Frank vector to homochoric vector."""
        if iszero(np.sum(ro[0:3]**2.0)):
            ho = [ 0.0, 0.0, 0.0 ]
        else:
            f = 2.0*np.arctan(ro[3]) -np.sin(2.0*np.arctan(ro[3])) if np.isfinite(ro[3]) else np.pi
            ho = ro[0:3] * (0.75*f)**(1.0/3.0)
        return np.array(ho)

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

    @staticmethod
    def ho2ro(ho):
        """Axis angle pair to Rodrigues-Frank vector."""
        return Rotation.ax2ro(Rotation.ho2ax(ho))

    @staticmethod
    def ho2cu(ho):
        """Homochoric vector to cubochoric vector."""
        return ball_to_cube(ho)


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
        """Cubochoric vector to homochoric vector."""
        return cube_to_ball(cu)
