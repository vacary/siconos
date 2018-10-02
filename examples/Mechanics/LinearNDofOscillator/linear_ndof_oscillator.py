#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2018 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

from numpy.linalg import norm
from siconos.kernel import LagrangianLinearTIDS, LagrangianDS,  NewtonImpactNSL,\
    LagrangianLinearTIR, Interaction, NonSmoothDynamicalSystem, MoreauJeanOSI,\
    TimeDiscretisation, LCP, TimeStepping, SiconosVector

from siconos.kernel import SimpleMatrix, getMatrix, SPARSE, DENSE
import numpy as np
import math

from numpy import eye, empty, float64, zeros

# User-defined main parameters
nDof = 2  #degrees of freedom for the beam
t0 = 1e-8    #initial computation time
T = 50.0                  # final computation time
h = 1e-2                # time step

position_init = 0.00000      # initial position
velocity_init =  0.0      # initial velocity

epsilon = 0.5#1e-1
theta = 1/2.0 + epsilon              # theta for MoreauJeanOSI integrator
#theta = 1.0

omega=10
e=0.0

class ForcedNdofOscillator(LagrangianDS):
    def __init__(self,x, v, nDof, omega):
        self._nDof = nDof
        self._omega= omega
        from scipy.sparse import csr_matrix
        K = csr_matrix((nDof, nDof))
        M = csr_matrix((nDof, nDof))
        
        K[0,0] = 1.
        K[0,1] = -1.
        M[0,0] = 1.
        
        for i in range(1,nDof-1):
            K[i,i] = 2.
            K[i,i-1] = -1.
            K[i,i+1] = -1.
            M[i,i] = 1.0
            
        K[nDof-1,nDof-2] = -1.
        K[nDof-1,nDof-1] = 1.
        M[nDof-1,nDof-1] = 1.0

        print(K.toarray())
        print(M.toarray())

        self._K_numpy = K.toarray()
        LagrangianDS.__init__(self,x,v,M.toarray())
        self.setKPtr(K.toarray())
        self.setFExtPtr(SiconosVector(nDof))


    # def jacobianqForces(self):
    #     return self._K_numpy
        
        

    def computeFExt(self, time, fExt=None):
        print('######### call computeFExt ########## with time =', time, "and fext = ", fExt  )
        #self.display(False);
        A=1.0
        if fExt is None :
            fExt = self._fExt
        if isinstance(fExt,SiconosVector):
            fExt.zero()
            fExt.setValue( self._nDof -1,  A*math.cos(self._omega*time))
        else:
            fExt[:]=0
            fExt[self._nDof -1]=  A*math.cos(self._omega*time)


# class ContactLastMass(Interaction):
#      def __init__(self, e):
#          self._e =e
#          H = np.zeros((1,nDof))
#          H[0,0]=1.
         
#          nslaw = NewtonImpactNSL(self._e)
#          relation = LagrangianLinearTIR(H)
#          Interaction.__init__(self,nslaw, relation)
         
# from scipy.sparse import csr_matrix

q0 = np.full((nDof), position_init)
v0 = np.full((nDof), velocity_init)

ndofOscillator = ForcedNdofOscillator(q0,v0,nDof, omega)
# inter = ContactLastMass(e)


# -------------
# --- Model ---
# -------------
impactingNdofOscillator = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system in the non smooth dynamical system
impactingNdofOscillator.insertDynamicalSystem(ndofOscillator);

# # link the interaction and the dynamical system
# impactingNdofOscillator.link(inter,ndofOscillator);

# ------------------
# --- Simulation ---
# ------------------

# -- (1) OneStepIntegrators --
OSI = MoreauJeanOSI(theta,0.5)

# -- (2) Time discretisation --
t = TimeDiscretisation(t0,h)

# -- (3) one step non smooth problem
osnspb = LCP()

s = TimeStepping(impactingNdofOscillator, t,OSI,osnspb)

k =0

N = int((T-t0)/h)
dataPlot = np.zeros((N+1, 8))

q = ndofOscillator.q()
v = ndofOscillator.velocity()
p = ndofOscillator.p(1)
#lambda_ = inter.lambda_(1)

# time loop
while s.hasNextEvent():
    s.computeOneStep()

    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0]/h
    #dataPlot[k, 4] = lambda_[0]
    dataPlot[k, 5] = q[nDof-1]
    dataPlot[k, 6] = v[nDof-1]
    dataPlot[k, 7] = p[nDof-1]/h

    k += 1
    s.nextStep()

dataPlot.resize(k,8)

import matplotlib.pyplot as plt

fig_size = [14, 14]
plt.rcParams["figure.figsize"] = fig_size

plt.subplot(411)
plt.title('displacement')
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title('velocity')
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title('reaction')
plt.grid()
plt.subplot(414)
plt.plot(dataPlot[:, 0], dataPlot[:, 4])
plt.title('lambda')
plt.grid()

plt.figure()
plt.subplot(411)
plt.title('displacement')
plt.plot(dataPlot[:, 0], dataPlot[:, 5])
plt.grid()
plt.subplot(412)
plt.title('velocity')
plt.plot(dataPlot[:, 0], dataPlot[:, 6])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 7])
plt.title('reaction')
plt.grid()

plt.show()
