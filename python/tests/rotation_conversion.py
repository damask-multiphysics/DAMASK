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

import numpy as np

from damask import _rotation

_P = _rotation._P

# parameters for conversion from/to cubochoric
_sc   = _rotation._sc
_beta = _rotation._beta
_R1   = _rotation._R1

def iszero(a):
    return np.isclose(a,0.0,atol=1.0e-12,rtol=0.0)

#---------- Quaternion ----------
def qu2om(qu):
    """Quaternion to rotation matrix."""
    qq = qu[0]**2-(qu[1]**2 + qu[2]**2 + qu[3]**2)
    om = np.diag(qq + 2.0*np.array([qu[1],qu[2],qu[3]])**2)

    om[0,1] = 2.0*(qu[2]*qu[1]+qu[0]*qu[3])
    om[1,0] = 2.0*(qu[1]*qu[2]-qu[0]*qu[3])
    om[1,2] = 2.0*(qu[3]*qu[2]+qu[0]*qu[1])
    om[2,1] = 2.0*(qu[2]*qu[3]-qu[0]*qu[1])
    om[2,0] = 2.0*(qu[1]*qu[3]+qu[0]*qu[2])
    om[0,2] = 2.0*(qu[3]*qu[1]-qu[0]*qu[2])
    return om if _P < 0.0 else np.swapaxes(om,-1,-2)

def qu2eu(qu):
    """Quaternion to Bunge-Euler angles."""
    q03 = qu[0]**2+qu[3]**2
    q12 = qu[1]**2+qu[2]**2
    chi = np.sqrt(q03*q12)
    if   np.abs(q12) < 1.e-8:
        eu = np.array([np.arctan2(-_P*2.0*qu[0]*qu[3],qu[0]**2-qu[3]**2), 0.0,   0.0])
    elif np.abs(q03) < 1.e-8:
        eu = np.array([np.arctan2(   2.0*qu[1]*qu[2],qu[1]**2-qu[2]**2), np.pi, 0.0])
    else:
        eu = np.array([np.arctan2((-_P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-_P*qu[0]*qu[1]-qu[2]*qu[3])*chi ),
                       np.arctan2( 2.0*chi, q03-q12 ),
                       np.arctan2(( _P*qu[0]*qu[2]+qu[1]*qu[3])*chi, (-_P*qu[0]*qu[1]+qu[2]*qu[3])*chi )])
    # reduce Euler angles to definition range
    eu[np.abs(eu)<1.e-6] = 0.0
    eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
    return eu

def qu2ax(qu):
    """
    Quaternion to axis angle pair.

    Modified version of the original formulation, should be numerically more stable
    """
    if np.abs(np.sum(qu[1:4]**2)) < 1.e-6:                                                          # set axis to [001] if the angle is 0/360
        ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
    elif qu[0] > 1.e-6:
        s = np.sign(qu[0])/np.sqrt(qu[1]**2+qu[2]**2+qu[3]**2)
        omega = 2.0 * np.arccos(np.clip(qu[0],-1.0,1.0))
        ax = ax = np.array([ qu[1]*s, qu[2]*s, qu[3]*s, omega ])
    else:
        ax = ax = np.array([ qu[1], qu[2], qu[3], np.pi])
    return ax

def qu2ro(qu):
    """Quaternion to Rodrigues-Frank vector."""
    if iszero(qu[0]):
        ro = np.array([qu[1], qu[2], qu[3], np.inf])
    else:
        s = np.linalg.norm(qu[1:4])
        ro = np.array([0.0,0.0,_P,0.0] if iszero(s) else \
                      [ qu[1]/s,  qu[2]/s,  qu[3]/s, np.tan(np.arccos(np.clip(qu[0],-1.0,1.0)))])
    return ro

def qu2ho(qu):
    """Quaternion to homochoric vector."""
    omega = 2.0 * np.arccos(np.clip(qu[0],-1.0,1.0))
    if np.abs(omega) < 1.0e-12:
        ho = np.zeros(3)
    else:
        ho = np.array([qu[1], qu[2], qu[3]])
        f  = 0.75 * ( omega - np.sin(omega) )
        ho = ho/np.linalg.norm(ho) * f**(1./3.)
    return ho


#---------- Rotation matrix ----------
def om2qu(a):
    trace = a[0,0] + a[1,1] + a[2,2]
    if trace > 0:
        s = 0.5 / np.sqrt(trace+ 1.0)
        qu = np.array([0.25 / s,( a[2,1] - a[1,2] ) * s,( a[0,2] - a[2,0] ) * s,( a[1,0] - a[0,1] ) * s])
    else:
        if ( a[0,0] > a[1,1] and a[0,0] > a[2,2] ):
            s = 2.0 * np.sqrt( 1.0 + a[0,0] - a[1,1] - a[2,2])
            qu = np.array([ (a[2,1] - a[1,2]) / s,0.25 * s,(a[0,1] + a[1,0]) / s,(a[0,2] + a[2,0]) / s])
        elif (a[1,1] > a[2,2]):
            s = 2.0 * np.sqrt( 1.0 + a[1,1] - a[0,0] - a[2,2])
            qu = np.array([ (a[0,2] - a[2,0]) / s,(a[0,1] + a[1,0]) / s,0.25 * s,(a[1,2] + a[2,1]) / s])
        else:
            s = 2.0 * np.sqrt( 1.0 + a[2,2] - a[0,0] - a[1,1] )
            qu = np.array([ (a[1,0] - a[0,1]) / s,(a[0,2] + a[2,0]) / s,(a[1,2] + a[2,1]) / s,0.25 * s])
    return qu*np.array([1.,_P,_P,_P])

def om2eu(om):
    """Rotation matrix to Bunge-Euler angles."""
    if not np.isclose(np.abs(om[2,2]),1.0,1.e-4):
        zeta = 1.0/np.sqrt(1.0-om[2,2]**2)
        eu = np.array([np.arctan2(om[2,0]*zeta,-om[2,1]*zeta),
                       np.arccos(om[2,2]),
                       np.arctan2(om[0,2]*zeta, om[1,2]*zeta)])
    else:
        eu = np.array([np.arctan2( om[0,1],om[0,0]), np.pi*0.5*(1-om[2,2]),0.0])                    # following the paper, not the reference implementation
    eu[np.abs(eu)<1.e-6] = 0.0
    eu = np.where(eu<0, (eu+2.0*np.pi)%np.array([2.0*np.pi,np.pi,2.0*np.pi]),eu)
    return eu

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
        diagDelta = -_P*np.array([om[1,2]-om[2,1],om[2,0]-om[0,2],om[0,1]-om[1,0]])
        diagDelta[np.abs(diagDelta)<1.e-6] = 1.0
        ax[0:3] = np.where(np.abs(diagDelta)<0, ax[0:3],np.abs(ax[0:3])*np.sign(diagDelta))
    return ax

#---------- Bunge-Euler angles ----------
def eu2qu(eu):
    """Bunge-Euler angles to quaternion."""
    ee = 0.5*eu
    cPhi = np.cos(ee[1])
    sPhi = np.sin(ee[1])
    qu = np.array([    cPhi*np.cos(ee[0]+ee[2]),
                   -_P*sPhi*np.cos(ee[0]-ee[2]),
                   -_P*sPhi*np.sin(ee[0]-ee[2]),
                   -_P*cPhi*np.sin(ee[0]+ee[2]) ])
    if qu[0] < 0.0: qu*=-1
    return qu

def eu2om(eu):
    """Bunge-Euler angles to rotation matrix."""
    c = np.cos(eu)
    s = np.sin(eu)

    om = np.array([[+c[0]*c[2]-s[0]*s[2]*c[1], +s[0]*c[2]+c[0]*s[2]*c[1], +s[2]*s[1]],
                   [-c[0]*s[2]-s[0]*c[2]*c[1], -s[0]*s[2]+c[0]*c[2]*c[1], +c[2]*s[1]],
                   [+s[0]*s[1],                -c[0]*s[1],                +c[1]     ]])
    om[np.abs(om)<1.e-12] = 0.0
    return om

def eu2ax(eu):
    """Bunge-Euler angles to axis angle pair."""
    t = np.tan(eu[1]*0.5)
    sigma = 0.5*(eu[0]+eu[2])
    delta = 0.5*(eu[0]-eu[2])
    tau   = np.linalg.norm([t,np.sin(sigma)])
    alpha = np.pi if iszero(np.cos(sigma)) else \
            2.0*np.arctan(tau/np.cos(sigma))

    if np.abs(alpha)<1.e-6:
        ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
    else:
        ax = -_P/tau * np.array([ t*np.cos(delta), t*np.sin(delta), np.sin(sigma) ])                # passive axis angle pair so a minus sign in front
        ax = np.append(ax,alpha)
        if alpha < 0.0: ax *= -1.0                                                                  # ensure alpha is positive
    return ax

def eu2ro(eu):
    """Bunge-Euler angles to Rodrigues-Frank vector."""
    ro = eu2ax(eu)                                                                                  # convert to axis angle pair representation
    if ro[3] >= np.pi:                                                                              # Differs from original implementation. check convention 5
        ro[3] = np.inf
    elif iszero(ro[3]):
        ro = np.array([ 0.0, 0.0, _P, 0.0 ])
    else:
        ro[3] = np.tan(ro[3]*0.5)
    return ro

#---------- Axis angle pair ----------
def ax2qu(ax):
    """Axis angle pair to quaternion."""
    if np.abs(ax[3])<1.e-6:
        qu = np.array([ 1.0, 0.0, 0.0, 0.0 ])
    else:
        c = np.cos(ax[3]*0.5)
        s = np.sin(ax[3]*0.5)
        qu = np.array([ c, ax[0]*s, ax[1]*s, ax[2]*s ])
    return qu

def ax2om(ax):
    """Axis angle pair to rotation matrix."""
    c = np.cos(ax[3])
    s = np.sin(ax[3])
    omc = 1.0-c
    om=np.diag(ax[0:3]**2*omc + c)

    for idx in [[0,1,2],[1,2,0],[2,0,1]]:
        q = omc*ax[idx[0]] * ax[idx[1]]
        om[idx[0],idx[1]] = q + s*ax[idx[2]]
        om[idx[1],idx[0]] = q - s*ax[idx[2]]
    return om if _P < 0.0 else np.swapaxes(om,(-1,-2))

def ax2ro(ax):
    """Axis angle pair to Rodrigues-Frank vector."""
    if np.abs(ax[3])<1.e-6:
        ro = [ 0.0, 0.0, _P, 0.0 ]
    else:
        ro = [ax[0], ax[1], ax[2]]
        # 180 degree case
        ro += [np.inf] if np.isclose(ax[3],np.pi,atol=1.0e-15,rtol=0.0) else \
              [np.tan(ax[3]*0.5)]
    ro = np.array(ro)
    return ro

def ax2ho(ax):
    """Axis angle pair to homochoric vector."""
    f = (0.75 * ( ax[3] - np.sin(ax[3]) ))**(1.0/3.0)
    ho = ax[0:3] * f
    return ho


#---------- Rodrigues-Frank vector ----------
def ro2ax(ro):
    """Rodrigues-Frank vector to axis angle pair."""
    if np.abs(ro[3]) < 1.e-6:
        ax = np.array([ 0.0, 0.0, 1.0, 0.0 ])
    elif not np.isfinite(ro[3]):
        ax = np.array([ ro[0], ro[1], ro[2], np.pi ])
    else:
        angle = 2.0*np.arctan(ro[3])
        ta = np.linalg.norm(ro[0:3])
        ax = np.array([ ro[0]*ta, ro[1]*ta, ro[2]*ta, angle ])
    return ax

def ro2ho(ro):
    """Rodrigues-Frank vector to homochoric vector."""
    if np.sum(ro[0:3]**2.0) < 1.e-6:
        ho = np.zeros(3)
    else:
        f = 2.0*np.arctan(ro[3]) -np.sin(2.0*np.arctan(ro[3])) if np.isfinite(ro[3]) else np.pi
        ho = ro[0:3] * (0.75*f)**(1.0/3.0)
    return ho

#---------- Homochoric vector----------
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

def ho2cu(ho):
    """
    Homochoric vector to cubochoric vector.

    References
    ----------
    D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
    https://doi.org/10.1088/0965-0393/22/7/075013

    """
    rs = np.linalg.norm(ho)

    if np.allclose(ho,0.0,rtol=0.0,atol=1.0e-16):
        cu = np.zeros(3)
    else:
        xyz3 = ho[_get_pyramid_order(ho,'forward')]

        # inverse M_3
        xyz2 = xyz3[0:2] * np.sqrt( 2.0*rs/(rs+np.abs(xyz3[2])) )

        # inverse M_2
        qxy = np.sum(xyz2**2)

        if np.isclose(qxy,0.0,rtol=0.0,atol=1.0e-16):
            Tinv = np.zeros(2)
        else:
            q2 = qxy + np.max(np.abs(xyz2))**2
            sq2 = np.sqrt(q2)
            q = (_beta/np.sqrt(2.0)/_R1) * np.sqrt(q2*qxy/(q2-np.max(np.abs(xyz2))*sq2))
            tt = np.clip((np.min(np.abs(xyz2))**2+np.max(np.abs(xyz2))*sq2)/np.sqrt(2.0)/qxy,-1.0,1.0)
            Tinv = np.array([1.0,np.arccos(tt)/np.pi*12.0]) if np.abs(xyz2[1]) <= np.abs(xyz2[0]) else \
                   np.array([np.arccos(tt)/np.pi*12.0,1.0])
            Tinv = q * np.where(xyz2<0.0,-Tinv,Tinv)

        # inverse M_1
        cu = np.array([ Tinv[0], Tinv[1],  (-1.0 if xyz3[2] < 0.0 else 1.0) * rs / np.sqrt(6.0/np.pi) ]) /_sc
        cu = cu[_get_pyramid_order(ho,'backward')]
    return cu

#---------- Cubochoric ----------
def cu2ho(cu):
    """
    Cubochoric vector to homochoric vector.

    References
    ----------
    D. Roşca et al., Modelling and Simulation in Materials Science and Engineering 22:075013, 2014
    https://doi.org/10.1088/0965-0393/22/7/075013

    """
    # transform to the sphere grid via the curved square, and intercept the zero point
    if np.allclose(cu,0.0,rtol=0.0,atol=1.0e-16):
        ho = np.zeros(3)
    else:
        # get pyramide and scale by grid parameter ratio
        XYZ = cu[_get_pyramid_order(cu,'forward')] * _sc

        # intercept all the points along the z-axis
        if np.allclose(XYZ[0:2],0.0,rtol=0.0,atol=1.0e-16):
            ho = np.array([0.0, 0.0, np.sqrt(6.0/np.pi) * XYZ[2]])
        else:
            order = [1,0] if np.abs(XYZ[1]) <= np.abs(XYZ[0]) else [0,1]
            q = np.pi/12.0 * XYZ[order[0]]/XYZ[order[1]]
            c = np.cos(q)
            s = np.sin(q)
            q = _R1*2.0**0.25/_beta * XYZ[order[1]] / np.sqrt(np.sqrt(2.0)-c)
            T = np.array([ (np.sqrt(2.0)*c - 1.0), np.sqrt(2.0) * s]) * q

            # transform to sphere grid (inverse Lambert)
            # note that there is no need to worry about dividing by zero, since XYZ[2] can not become zero
            c = np.sum(T**2)
            s = c *         np.pi/24.0 /XYZ[2]**2
            c = c * np.sqrt(np.pi/24.0)/XYZ[2]

            q = np.sqrt( 1.0 - s )
            ho = np.array([ T[order[1]] * q, T[order[0]] * q, np.sqrt(6.0/np.pi) * XYZ[2] - c ])

        ho = ho[_get_pyramid_order(cu,'backward')]
    return ho

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
    order = {'forward':np.array([[0,1,2],[1,2,0],[2,0,1]]),
             'backward':np.array([[0,1,2],[2,0,1],[1,2,0]])}
    if   np.maximum(abs(xyz[0]),abs(xyz[1])) <= xyz[2] or \
         np.maximum(abs(xyz[0]),abs(xyz[1])) <=-xyz[2]:
        p = 0
    elif np.maximum(abs(xyz[1]),abs(xyz[2])) <= xyz[0] or \
         np.maximum(abs(xyz[1]),abs(xyz[2])) <=-xyz[0]:
        p = 1
    elif np.maximum(abs(xyz[2]),abs(xyz[0])) <= xyz[1] or \
         np.maximum(abs(xyz[2]),abs(xyz[0])) <=-xyz[1]:
        p = 2
    return order[direction][p]
