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

[1] colocar a ref.

"""

import quadrature

import numpy as np
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

        print("to be continued ...")

    def build3d(self):
        '''
        Build the quadrature set on the unit sphere.
        '''
        self.nnodes3 = 8*self.nnodes
        self.mu3 = np.resize(self.mu3, self.nnodes3)
        self.eta3 = np.resize(self.eta3, self.nnodes3)
        self.xi3 = np.resize(self.xi3, self.nnodes3)
        self.w3 = np.resize(self.w3, self.nnodes3)
        node3 = 0
        for i in [1,-1]:
            for j in [1,-1]:
                for k in [1,-1]:
                    for node in range(self.nnodes):
                        self.mu3[node3] = i*self.mu[node]
                        self.eta3[node3] = j*self.eta[node]
                        self.xi3[node3] = k*self.xi[node]
                        self.w3[node3] = self.w[node]
                        node3 += 1
                        
