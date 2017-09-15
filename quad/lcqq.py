#!/usr/bin/python3
"""
Legendre-Chebyshev quadrangular quadrature (P_N T_N) quadrature set on the unit sphere

This Python3 code is a inherited class of quadrature sets on the unit
sphere, i.e.

..math: 
   \\int_{S^2} f(\\Omega) d\\Omega \\approx \\sum_{i=1}^M w_i f(\\Omega_i),

where 

..math: 
   \\Omega = (\\mu, \\eta, \\xi)\in S^2 := \\{\\mu^2+\\eta^2+\\xi^2=1\\}.

Author
______

Pedro H A Konzen - UFRGS - Aug/2017
Ana Carolina Bof - UFRGS - Aug/2017

License
_______

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

References
__________

[1] Longoni, G., Haghighat, A. (2001). Development of new quadrature sets with the ordinate splitting technique. 
Proceedings the 2001 American Nuclear Society International Meeting on Mathematical Methods for Nuclear 
Applications (M&C 2001), Salt Lake City, UT.

"""

import quadrature

import numpy as np
from numpy import polynomial as poly
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class Lcqq(quadrature.Quadrature):

    def __init__(self,N=0):
        quadrature.Quadrature.__init__(self,N)

        #verify if N is even
        assert (N%2 == 0), "N must be even."

        #num. nodes of first octant
        self.nnodes = int(N**2/4)

        #nodes and weights on first octant
        self.mu = np.resize(self.mu, self.nnodes)
        self.eta = np.resize(self.eta, self.nnodes)
        self.xi = np.resize(self.xi, self.nnodes)
        self.w = np.resize(self.w, self.nnodes)

        #build quadrature set on first octant
        self.buildFirstOctantSet()
 
        #build quadrature set on first octant
        self.buildFirstQuadrantSet()

    def __repr__(self):
        note =  "\nLegendre-Chebyshev quadrangular quadrature set\n"
        note += "Author:\n" +\
                "\tPedro H A Konzen - UFRGS - 2017\n" +\
                "\tAna Carolina Bof\n" +\
                "This program comes with ABSOLUTELY NO WARRANTY.\n\n" +\
                "This program is free software: you can redistribute" +\
                "it and/or modify\n" +\
                "it under the terms of the GNU General Public" +\
                "License as published by\n" +\
                "the Free Software Foundation, either version 3" +\
                "of the License, or\n" +\
                "(at your option) any later version.\n\n"
        return note

    def buildFirstOctantSet(self):
        '''
        Build the quadrature set on the first octant.
        '''
        #Gauss-Legendre quadrature (x,w) - (points, weights)
        x,wi = poly.legendre.leggauss(self.N);
        #just the positive values
        x0 = x[int(self.N/2):] 
        w0 = wi[int(self.N/2):]

        #wbar
        wbar = np.zeros(self.N, dtype='double')
        for j in np.arange(1,self.N+1):
            wbar[j-1] = np.pi/2 * (1 - (self.N-2*j+1)/self.N)
            
        #Quadratute points on the first octant (3D computations)
        c = 0
        for i in np.arange(int(self.N/2)):
            for j in np.arange(int(self.N/2)):
                aux = np.sqrt(1-x0[i]**2)*np.cos(wbar[j])
                if (aux > 0.0):
                    self.mu[c] = aux
                    self.eta[c] = np.sqrt(1 - self.mu[c]**2-x0[i]**2)
                    self.xi[c] = np.sqrt(1 - self.mu[c]**2 - self.eta[c]**2)
                    self.w[c] = np.pi * w0[i]/self.N
                    c += 1
        
    def buildFirstQuadrantSet(self):
        '''
        Build the quadrature set on the first quadrant.
        '''
        #Gauss-Legendre quadrature (x,w) - (points, weights)
        x,wi = poly.legendre.leggauss(self.N);
        #just the positive values
        x0 = x[int(self.N/2):] 
        w0 = wi[int(self.N/2):]

        #wbar
        wbar = np.zeros(self.N, dtype='double')
        for j in np.arange(1,self.N+1):
            wbar[j-1] = np.pi/2 * (1 - (self.N-2*j+1)/self.N)
        
        #Quadratute points on the first quadrant (2D computations)
        c = 0
        for i in np.arange(int(self.N/2)):
            for j in np.arange(int(self.N/2)):
                aux = np.sqrt(1-x0[i]**2)*np.cos(wbar[j])
                if (aux > 0.0):
                    self.mu[c] = aux
                    self.eta[c] = np.sqrt(1 - self.mu[c]**2-x0[i]**2)
                    self.w[c] = 2*np.pi * w0[i]/self.N
                    c += 1
                            
