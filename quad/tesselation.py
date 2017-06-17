#!/usr/bin/python3
"""
Tesselation (T_N) quadrature set on the unit sphere

This Python3 code is a inherited class of quadrature sets on the unit
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

[1] Thurgood, C. P., A. Pollard, and H. A. Becker. "The TN quadrature set for the discrete ordinates method." TRANSACTIONS-AMERICAN SOCIETY OF MECHANICAL ENGINEERS JOURNAL OF HEAT TRANSFER 117 (1995): 1068-1069.

"""

import quadrature

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class Tesselation(quadrature.Quadrature):

    def __init__(self,N=0):
        quadrature.Quadrature.__init__(self,N)
        #num. nodes of first octant
        self.nnodes = N**2
        #notes and weights on first octant
        self.mu = np.resize(self.mu, self.nnodes)
        self.eta = np.resize(self.eta, self.nnodes)
        self.xi = np.resize(self.xi, self.nnodes)
        self.w = np.resize(self.w, self.nnodes)
        #centroids
        self.centroid = []
        self.strT = []
        self.sphT = []
        #build quadrature set on first octant
        self.buildFirstOctantSet()
        

    def __repr__(self):
        note =  "\nTesselation quadrature set\n"
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
        #vertices on (projected) straight triangle
        nv = 3*self.N
        #vertical step
        h = 1.0/self.N
        #level list of dictionaries
        level = []
        for l in range(self.N+1):
            #add new level
            level.append({})
            #num verteces
            nv = self.N+1-l
            level[l]['nv'] = nv
            #inicial vertice
            iv = [1-l*h,0,l*h]
            level[l]['iv'] = iv
            #final vertice
            fv = [0,1-l*h,l*h]
            level[l]['fv'] = fv
            #vertices
            level[l]['v'] = []
            for i in range(nv):
                level[l]['v'].append([iv[0]-i*h,iv[1]+i*h,iv[2]])

        #tesselate the straight triangle
        t = 0
        self.strT = np.zeros((self.nnodes,3,3),dtype='double')
        #loop over levels
        for l in range(self.N):
            for v in range(level[l]['nv']-1):
                self.strT[t,0,:] = level[l]['v'][v]
                self.strT[t,1,:] = [level[l]['v'][v][0]-h,
                               level[l]['v'][v][1]+h,
                               level[l]['v'][v][2]]
                self.strT[t,2,:] = [level[l]['v'][v][0]-h,
                               level[l]['v'][v][1],
                               level[l]['v'][v][2]+h]
                t += 1
            for v in range(1,level[l]['nv']-1):
                self.strT[t,0,:] = level[l]['v'][v]
                self.strT[t,1,:] = [level[l]['v'][v][0]-h,
                               level[l]['v'][v][1],
                               level[l]['v'][v][2]+h]
                self.strT[t,2,:] = [level[l]['v'][v][0],
                               level[l]['v'][v][1]-h,
                               level[l]['v'][v][2]+h]
                t += 1
        #compute straight triangle centroids
        self.centroid = np.empty((self.nnodes,3), dtype='double')
        for t in range(self.nnodes):
            for c in range(3):
                self.centroid[t,c] = (self.strT[t,0,c]+self.strT[t,1,c]+self.strT[t,2,c])/3

        #project straight triangles onto spherical triangles
        #also computes the quadrature nodes (first octant)
        self.sphT = np.empty((self.nnodes,3,3), dtype='double')
        for t in range(self.nnodes):
            for p in range(3):
                factor = np.sqrt(self.strT[t,p,0]**2 +
                                 self.strT[t,p,1]**2 +
                                 self.strT[t,p,2]**2)
                for c in range(3):
                    self.sphT[t,p,c] = self.strT[t,p,c]/factor
        
            factor = np.sqrt(self.centroid[t,0]**2 +
                             self.centroid[t,1]**2 +
                             self.centroid[t,2]**2)
            self.mu[t] = self.centroid[t,0]/factor
            self.eta[t] = self.centroid[t,1]/factor
            self.xi[t] = self.centroid[t,2]/factor

        #compute the quadrature weights (first octant), i.e.
        #the surface of each spherical triangle
        for t in range(self.nnodes):
            #compute arc lengths
            a = 0
            b = 0
            c = 0
            for cc in range(3):
                a += self.sphT[t,1,cc]*self.sphT[t,2,cc]
                b += self.sphT[t,0,cc]*self.sphT[t,2,cc]
                c += self.sphT[t,0,cc]*self.sphT[t,1,cc]
            a = np.arccos(a)
            b = np.arccos(b)
            c = np.arccos(c)
            #compute vertices angles
            A = np.arccos((np.cos(a) - np.cos(b)*np.cos(c))/\
                        (np.sin(b)*np.sin(c)))
            B = np.arccos((np.cos(b) - np.cos(c)*np.cos(a))/\
                        (np.sin(c)*np.sin(a)))
            C = np.arccos((np.cos(c) - np.cos(a)*np.cos(b))/\
                        (np.sin(a)*np.sin(b)))
            #compute the area as the spherical excess
            self.w[t] = A+B+C - np.pi

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
                        
                        
    def plotFirstOctant(self,
                        fname="plot",
                        extension="png",
                        showAxis=True,
                        show=False):
        #font letter
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif', size=12)
        

        fig = plt.figure(figsize=(8,4), dpi=300, 
                         linewidth=0.0, facecolor="white")
        ax1 = plt.subplot(1,2,1, projection='3d')
        ax2 = plt.subplot(1,2,2, projection='3d')
        #ax.axis('off')

        ax1.set_xlabel("x")
        ax1.set_xticks([0,1])
        ax1.set_ylabel("y")
        ax1.set_yticks([0])
        ax1.set_zlabel("z")
        ax1.set_zticks([1])
        
        ax2.set_xlabel("x")
        ax2.set_xticks([0,1])
        ax2.set_ylabel("y")
        ax2.set_yticks([0])
        ax2.set_zlabel("z")
        ax2.set_zticks([1])
        
        ax1.view_init(15,45)
        ax2.view_init(15,45)

        #plot's list
        p = []
        for t in range(self.nnodes):
            #straigth triangles
            p.append(ax1.plot([self.strT[t,0,0],self.strT[t,1,0],self.strT[t,2,0],self.strT[t,0,0]],
                              [self.strT[t,0,1],self.strT[t,1,1],self.strT[t,2,1],self.strT[t,0,1]],
                              [self.strT[t,0,2],self.strT[t,1,2],self.strT[t,2,2],self.strT[t,0,2]],
                              color="black",linestyle="--",linewidth=0.75))        

            #straigth triangles centroids
            p.append(ax1.plot([self.centroid[t,0]],
                              [self.centroid[t,1]],
                              [self.centroid[t,2]], 'k.', color="blue"))

            #numbering
            p.append(ax1.text(self.centroid[t,0],self.centroid[t,1]+0.05,
                              self.centroid[t,2],str(t),fontsize=8))

            #spherical triangles
            Npts = 30
            xx = np.linspace(self.strT[t,0,0],self.strT[t,1,0],Npts)
            xx = np.append(xx,np.linspace(self.strT[t,1,0],self.strT[t,2,0],Npts))
            xx = np.append(xx,np.linspace(self.strT[t,2,0],self.strT[t,0,0],Npts))

            yy = np.linspace(self.strT[t,0,1],self.strT[t,1,1],Npts)
            yy = np.append(yy,np.linspace(self.strT[t,1,1],self.strT[t,2,1],Npts))
            yy = np.append(yy,np.linspace(self.strT[t,2,1],self.strT[t,0,1],Npts))

            zz = np.linspace(self.strT[t,0,2],self.strT[t,1,2],Npts)
            zz = np.append(zz,np.linspace(self.strT[t,1,2],self.strT[t,2,2],Npts))
            zz = np.append(zz,np.linspace(self.strT[t,2,2],self.strT[t,0,2],Npts))

            for i in range(xx.size):
                factor = np.sqrt(xx[i]**2 + yy[i]**2 + zz[i]**2)
                xx[i] = xx[i]/factor
                yy[i] = yy[i]/factor
                zz[i] = zz[i]/factor

            p.append(ax2.plot(xx,yy,zz,color="black",
                          linestyle="--",linewidth=0.75))

            #spherical triangles centroids
            factor = np.sqrt(self.centroid[t,0]**2+self.centroid[t,1]**2+self.centroid[t,2]**2)
            p.append(ax2.plot([self.centroid[t,0]/factor+0.025],
                          [self.centroid[t,1]/factor],
                          [self.centroid[t,2]/factor], 'k.', color="blue"))

            #numbering
            p.append(ax2.text(self.centroid[t,0]/factor,self.centroid[t,1]/factor,
                      self.centroid[t,2]/factor,str(t),fontsize=8))

        #external straight triangle
        p.append(ax1.plot([1,0,0,1],[0,1,0,0],[0,0,1,0],"k-"))

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
        p.append(ax2.plot(xx,yy,zz,color="black",
                          linestyle="-",linewidth=0.75))

                        
        plt.savefig(fname+"."+extension,
                    bbox_inches="tight")
