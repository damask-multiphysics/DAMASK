####################################################################################################
# Code below available according to below conditions on https://github.com/MarDiehl/3Drotations
####################################################################################################
# Copyright (c) 2017-2019, Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
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

sc   = np.pi**(1./6.)/6.**(1./6.)
beta = np.pi**(5./6.)/6.**(1./6.)/2.
R1   = (3.*np.pi/4.)**(1./3.)

def CubeToBall(cube):
  
    if np.abs(np.max(cube))>np.pi**(2./3.) * 0.5:
      raise ValueError
  
    # transform to the sphere grid via the curved square, and intercept the zero point
    if np.allclose(cube,0.0,rtol=0.0,atol=1.0e-300):
      ball = np.zeros(3)
    else:
      # get pyramide and scale by grid parameter ratio
      p = GetPyramidOrder(cube)
      XYZ = cube[p] * sc

      # intercept all the points along the z-axis
      if np.allclose(XYZ[0:2],0.0,rtol=0.0,atol=1.0e-300):
        ball = np.array([0.0, 0.0, np.sqrt(6.0/np.pi) * XYZ[2]])
      else:
        order = [1,0] if np.abs(XYZ[1]) <= np.abs(XYZ[0]) else [0,1]
        q = np.pi/12.0 * XYZ[order[0]]/XYZ[order[1]]
        c = np.cos(q)
        s = np.sin(q)
        q = R1*2.0**0.25/beta * XYZ[order[1]] / np.sqrt(np.sqrt(2.0)-c)
        T = np.array([ (np.sqrt(2.0)*c - 1.0), np.sqrt(2.0) * s]) * q

        # transform to sphere grid (inverse Lambert)
        # note that there is no need to worry about dividing by zero, since XYZ[2] can not become zero
        c = np.sum(T**2)
        s = c *         np.pi/24.0 /XYZ[2]**2
        c = c * np.sqrt(np.pi/24.0)/XYZ[2]
        q = np.sqrt( 1.0 - s )
        ball = np.array([ T[order[1]] * q, T[order[0]] * q, np.sqrt(6.0/np.pi) * XYZ[2] - c ])
    
      # reverse the coordinates back to the regular order according to the original pyramid number
      ball = ball[p]

    return ball


def BallToCube(ball):
  
    rs = np.linalg.norm(ball)
    if rs > R1: 
      raise ValueError
  
    if np.allclose(ball,0.0,rtol=0.0,atol=1.0e-300):
      cube = np.zeros(3)
    else:
      p = GetPyramidOrder(ball)
      xyz3 = ball[p] 

      # inverse M_3
      xyz2 = xyz3[0:2] * np.sqrt( 2.0*rs/(rs+np.abs(xyz3[2])) )
      
      # inverse M_2
      qxy = np.sum(xyz2**2)
      
      if np.isclose(qxy,0.0,rtol=0.0,atol=1.0e-300):
        Tinv = np.zeros(2)
      else:
        q2 = qxy + np.max(np.abs(xyz2))**2
        sq2 = np.sqrt(q2)
        q = (beta/np.sqrt(2.0)/R1) * np.sqrt(q2*qxy/(q2-np.max(np.abs(xyz2))*sq2))
        tt = np.clip((np.min(np.abs(xyz2))**2+np.max(np.abs(xyz2))*sq2)/np.sqrt(2.0)/qxy,-1.0,1.0)
        Tinv = np.array([1.0,np.arccos(tt)/np.pi*12.0]) if np.abs(xyz2[1]) <= np.abs(xyz2[0]) else \
               np.array([np.arccos(tt)/np.pi*12.0,1.0])
        Tinv = q * np.where(xyz2<0.0,-Tinv,Tinv)
      
      # inverse M_1
      cube = np.array([ Tinv[0], Tinv[1],  (-1.0 if xyz3[2] < 0.0 else 1.0) * rs / np.sqrt(6.0/np.pi) ]) /sc

      # reverst the coordinates back to the regular order according to the original pyramid number
      cube = cube[p]
      
    return cube

def GetPyramidOrder(xyz):
 
    if   (abs(xyz[0])<= xyz[2]) and (abs(xyz[1])<= xyz[2]) or \
         (abs(xyz[0])<=-xyz[2]) and (abs(xyz[1])<=-xyz[2]):
      return [0,1,2]
    elif (abs(xyz[2])<= xyz[0]) and (abs(xyz[1])<= xyz[0]) or \
         (abs(xyz[2])<=-xyz[0]) and (abs(xyz[1])<=-xyz[0]):
      return [1,2,0]
    elif (abs(xyz[0])<= xyz[1]) and (abs(xyz[2])<= xyz[1]) or \
         (abs(xyz[0])<=-xyz[1]) and (abs(xyz[2])<=-xyz[1]):
      return [2,0,1]
