#!/usr/bin/python3
"""
Quadrature sets on the unit sphere

This Python3 code is a base class of quadrature sets on the unit
sphere, i.e.

..math: 
   \\int_{S^2} f(\\Omega) d\\Omega \\approx \\sum_{i=1}^M w_i f(\\Omega_i),

where 

..math: 
   \\Omega = (\\mu, \\eta, \\xi)\in S^2 := \\{\\mu^2+\\eta^2+\\xi^2=1\\}.

Author
______

Pedro H A Konzen - UFRGS - Mar/2017

Licence
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

"""

import numpy as np

class Quadrature:

    def __init__(self,N=0):
        assert N!=0, "You must provide a non zero quadrature order."

        #quadrature order
        self.N = N
        #quadrature set on first octant
        self.nnodes = 0
        self.mu = np.empty(self.nnodes, dtype='double')
        self.eta = np.empty(self.nnodes, dtype='double')
        self.xi = np.empty(self.nnodes, dtype='double')
        self.w = np.empty(self.nnodes, dtype='double')
        #complete quadrature set
        self.nnodes3 = 0
        self.mu3 = np.empty(self.nnodes3, dtype='double')
        self.eta3 = np.empty(self.nnodes3, dtype='double')
        self.xi3 = np.empty(self.nnodes3, dtype='double')
        self.w3 = np.empty(self.nnodes3, dtype='double')

    def __name__(self):
        return "\nBase class for quadrature sets on the unit sphere\n\
        Author: Pedro H A Konzen - UFRGS - 2017\n\n\
        This program comes with ABSOLUTELY NO WARRANTY.\n\n\
        This program is free software: you can redistribute\
        it and/or modify\n\
        it under the terms of the GNU General Public License\
        as published by\n\
        the Free Software Foundation, either version 3 of the\
        License, or\n\
        (at your option) any later version.\n\n"

    def printFirstOctantSet(self):
        print("%s" % ((4*15+3)*"*"))
        print("Quadrature set on first octant")
        print("mu              eta             xi              w")
        for i in range(self.nnodes):
            print("%1.9e %1.9e %1.9e %1.9e" % (self.mu[i],\
                                               self.eta[i],\
                                               self.xi[i],\
                                               self.w[i]))
        print("\nNum. nodes on first octant: %d\n" % self.nnodes)
        print("%s\n" % ((4*15+3)*"*"))

        
    def diagnostics(self):
        self.zerothMomentError()

    def zerothMomentError(self):
        s = 0.0
        for i in range(self.nnodes3):
            s += self.w3[i]
        print("Zeroth moment error: %1.2e" % np.fabs(s-4*np.pi))
        
