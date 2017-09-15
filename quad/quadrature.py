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
Ana Carolina Bof - UFRGS - Jun/2017

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

"""

import numpy as np
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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
        
        #Quadratute points on the first quadrant (2D computations)
        self.nnodes2 = 0
        self.mu2 = np.empty(self.nnodes2, dtype='double')
        self.eta2 = np.empty(self.nnodes2, dtype='double')
        self.w2 = np.empty(self.nnodes2, dtype='double')
        
        #complete quadrature set
        self.nnodes3 = 0
        self.mu3 = np.empty(self.nnodes3, dtype='double')
        self.eta3 = np.empty(self.nnodes3, dtype='double')
        self.xi3 = np.empty(self.nnodes3, dtype='double')
        self.w3 = np.empty(self.nnodes3, dtype='double')

    def __name__(self):
        return "\nBase class for quadrature sets on the unit sphere\n\
        Authors: \n\
        \tPedro H A Konzen - UFRGS - 2017\n\
        \tAna Carolina Bof- UFRGS - 2017\n\n\
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

    def plotFirstOctant(self,
                        fname="plot",
                        extension="png",
                        showAxis=True,
                        show=False):
        #font letter
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif', size=12)
        

        fig = plt.figure(figsize=(4,3), dpi=300, 
                         linewidth=0.0, facecolor="white")
        ax = plt.subplot(1,1,1, projection='3d')

        ax.set_xlabel("x")
        ax.set_xticks([0,1])
        ax.set_ylabel("y")
        ax.set_yticks([0])
        ax.set_zlabel("z")
        ax.set_zticks([1])

        ax.view_init(15,45)
        ax.view_init(15,45)

        #plot's list
        p = []

        #external spherical triangle
        Npts = 30
        xx = np.linspace(1,0,Npts)
        xx = np.append(xx,np.zeros(Npts))
        xx = np.append(xx,np.linspace(0,1,Npts))

        yy = np.linspace(0,1,Npts)
        yy = np.append(yy,np.linspace(1,0,Npts))
        yy = np.append(yy,np.zeros(Npts))

        zz = np.zeros(Npts)
        zz = np.append(zz,np.linspace(0,1,Npts))
        zz = np.append(zz,np.linspace(1,0,Npts))

        for i in range(xx.size):
            factor = np.sqrt(xx[i]**2 + yy[i]**2 + zz[i]**2)
            xx[i] = xx[i]/factor
            yy[i] = yy[i]/factor
            zz[i] = zz[i]/factor
        p.append(ax.plot(xx,yy,zz,color="black",
                          linestyle="-",linewidth=0.75))

        #plot quadrature nodes
        for i in range(self.nnodes):
            p.append(ax.plot([self.mu[i]],[self.eta[i]],[self.xi[i]],
                             marker=".", markersize=3,
                             color="red"))

        plt.savefig(fname+"."+extension,
                    bbox_inches="tight")
        if (show):
            plt.show()

            

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

                        
    def build2d(self):
        '''
        Build the quadrature set on the unit circle.
        '''
        self.nnodes2 = 4*self.nnodes
        self.mu2 = np.resize(self.mu2, self.nnodes2)
        self.eta2 = np.resize(self.eta2, self.nnodes2)
        self.w2 = np.resize(self.w2, self.nnodes2)
        
        node2 = 0

        for i in [1,-1]:
            for j in [1,-1]:
                for node in range(self.nnodes):
                    self.mu2[node2] = i*self.mu[node]
                    self.eta2[node2] = j*self.eta[node]
                    self.w2[node2] = self.w[node]
                    node2 += 1

    def printFirstQuadrantSet(self):
        print("%s" % ((4*15+3)*"*"))
        print("Quadrature set on first quadrant")
        print("mu              eta             w")
        for i in range(self.nnodes):
            print("%1.9e %1.9e %1.9e" % (self.mu[i],\
                                         self.eta[i],\
                                         self.w[i]))
        print("\nNum. nodes on first quadrant: %d\n" % self.nnodes)
        print("%s\n" % ((4*15+3)*"*"))
        
    def diagnostics(self):
        self.zerothMomentError()

    def zerothMomentError(self):
        s = 0.0
        for i in range(self.nnodes3):
            s += self.w3[i]
        print("Zeroth moment error: %1.2e" % np.fabs(s-4*np.pi))
