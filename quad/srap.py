#!/usr/bin/python3
"""
Sherical Rings Arithmetic Progression (SRAP_N) quadrature set 
on the unit sphere.

This Python3 code is a inherited class of quadrature sets on the unit
sphere, i.e.

..math: 
   \\int_{S^2} f(\\Omega) d\\Omega \\approx \\sum_{i=1}^M w_i f(\\Omega_i),

where 

..math: 
   \\Omega = (\\mu, \\eta, \\xi)\in S^2 := \\{\\mu^2+\\eta^2+\\xi^2=1\\}.

Author
______

Pedro H A Konzen - UFRGS - Sep/2017

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

[1] Thurgood, C. P., A. Pollard, and H. A. Becker. "The TN quadrature set for the discrete ordinates method." TRANSACTIONS-AMERICAN SOCIETY OF MECHANICAL ENGINEERS JOURNAL OF HEAT TRANSFER 117 (1995): 1068-1069.

"""

import quadrature

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class SRAP(quadrature.Quadrature):

    def __init__(self,N=0):
        quadrature.Quadrature.__init__(self,N)
        assert (self.N>=2), "Quadrature order N must be >= 2."
        #num. nodes of first octant
        self.nnodes = int(self.N*(self.N+3)/2)
        #nodes and weights on first octant
        self.mu = np.resize(self.mu, self.nnodes)
        self.eta = np.resize(self.eta, self.nnodes)
        self.xi = np.resize(self.xi, self.nnodes)
        self.w = np.resize(self.w, self.nnodes)
        #centroids
        self.centroid = []
        self.level = []
        #build quadrature set on first octant
        self.buildFirstOctantSet()
        

    def __repr__(self):
        note =  "\nSRAP_N quadrature set\n"
        note += "Author:\n" +\
                "\tPedro H A Konzen - UFRGS - 2017\n" +\
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
        M = int(self.N*(self.N+3)/2)
    
        self.level = []
        for l in np.arange(self.N):
            self.level.append({})
            m = l+2
            dphi = np.pi/2/m
            self.level[l]['ne'] = m
            self.level[l]['phi'] = []
            
            for i in np.arange(m+1):
                self.level[l]['phi'].append(i*np.pi/(2*m))

            if (l == 0):
                t1 = 0
                t2 = np.arccos(1-2/M)
            else:
                t1 = self.level[l-1]['theta'][1]
                t2 = np.arccos(np.cos(t1)-(l+2)/M)

            self.level[l]['theta'] = [t1,t2]
 
            #centroid of the solid angle
            self.level[l]['p'] = []
            for e in np.arange(m):
                t1 = self.level[l]['theta'][0]
                t2 = self.level[l]['theta'][1]
                p1 = self.level[l]['phi'][e]
                p2 = self.level[l]['phi'][e+1]
                mu = (np.sin(p2)-np.sin(p1)) * \
                  (t2/2-np.sin(2*t2)/4 -t1/2 + np.sin(2*t1)/4)
                eta = (np.cos(p1)-np.cos(p2)) * \
                  (t2/2-np.sin(2*t2)/4 -t1/2 + np.sin(2*t1)/4)
                xi = 1.0/2 * (p2-p1) * (np.sin(t2)**2-np.sin(t1)**2)

                f = np.sqrt(mu**2+eta**2+xi**2)
                self.level[l]['p'].append([mu/f,eta/f,xi/f])

        assert (abs(self.level[self.N-1]['theta'][1] - np.pi/2) < 1e-14)

        surf = np.pi/2/self.nnodes

        c=0
        for l in np.arange(self.N):
            for e in np.arange(self.level[l]['ne']):
                self.mu[c] = self.level[l]['p'][e][0]
                self.eta[c] = self.level[l]['p'][e][1]
                self.xi[c] = self.level[l]['p'][e][2]
                self.w[c] = surf
                c += 1
                
            
        
                      
                        
    def plotFirstOctant(self,
                        fname="plot",
                        extension="png",
                        showAxis=True,
                        show=False):
        #font letter
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif', size=12)
        

        fig = plt.figure(figsize=(8,4), dpi=300, \
                             linewidth=0.0, facecolor="white")
        ax1 = plt.subplot(1,1,1, projection='3d')
        #ax.axis('off')

        ax1.set_xlabel("x")
        ax1.set_xticks([0,1])
        ax1.set_ylabel("y")
        ax1.set_yticks([0])
        ax1.set_zlabel("z")
        ax1.set_zticks([1])
                
        ax1.view_init(15,45)

        #plot's list
        nps = 30
        for l in np.arange(self.N):
            for e in np.arange(self.level[l]['ne']):
                #point
                ax1.plot([self.level[l]['p'][e][0]],
                        [self.level[l]['p'][e][1]],
                        [self.level[l]['p'][e][2]],
                        color="blue", linestyle="-", linewidth=0.75)
        
                ti = self.level[l]['theta'][0]
                tf = self.level[l]['theta'][1]
                dt = (tf-ti)/nps

                pi = self.level[l]['phi'][e]
                pf = self.level[l]['phi'][e+1]
                dp = (pf-pi)/nps

                #line
                xx = []
                yy = []
                zz = []
                for i in np.arange(nps):
                    theta = ti+i*dt
                    phi = pi
                    mu = np.sin(theta)*np.cos(phi)
                    eta = np.sin(theta)*np.sin(phi)
                    xi = np.cos(theta)
                    xx.append(mu)
                    yy.append(eta)
                    zz.append(xi)
                ax1.plot(xx,yy,zz,color='black',linestyle='-',linewidth=0.75)

                #line
                xx = []
                yy = []
                zz = []
                for i in np.arange(nps):
                    theta = ti+i*dt
                    phi = pf
                    mu = np.sin(theta)*np.cos(phi)
                    eta = np.sin(theta)*np.sin(phi)
                    xi = np.cos(theta)
                    xx.append(mu)
                    yy.append(eta)
                    zz.append(xi)
                ax1.plot(xx,yy,zz,color='black',linestyle='-',linewidth=0.75)

                #line
                xx = []
                yy = []
                zz = []
                for i in np.arange(nps):
                    theta = ti
                    phi = pi+i*dp
                    mu = np.sin(theta)*np.cos(phi)
                    eta = np.sin(theta)*np.sin(phi)
                    xi = np.cos(theta)
                    xx.append(mu)
                    yy.append(eta)
                    zz.append(xi)
                ax1.plot(xx,yy,zz,color='black',linestyle='-',linewidth=0.75)

                #line
                xx = []
                yy = []
                zz = []
                for i in np.arange(nps):
                    theta = tf
                    phi = pi+i*dp
                    mu = np.sin(theta)*np.cos(phi)
                    eta = np.sin(theta)*np.sin(phi)
                    xi = np.cos(theta)
                    xx.append(mu)
                    yy.append(eta)
                    zz.append(xi)
                ax1.plot(xx,yy,zz,color='black',linestyle='-',linewidth=0.75)

        for node in np.arange(self.nnodes):
            ax1.plot([self.mu[node]],\
                         self.eta[node],\
                         self.xi[node],\
                         color='blue', marker='.')
            #numbering
            ax1.text(self.mu[node],self.eta[node]+0.05,\
                         self.xi[node],str(node+1),fontsize=8)


                        
        plt.savefig(fname+"."+extension,
                    bbox_inches="tight")
