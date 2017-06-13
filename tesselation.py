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

import quadrature

import numpy as np

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
        #build quadrature set on first octant
        self.buildFirstOctantSet()
        

    def __repr__(self):
        note =  "\nTesselation quadrature set\n"
        note += "Author: Pedro H A Konzen - UFRGS - 2017\n\n" +\
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
        strT = np.zeros((self.nnodes,3,3),dtype='double')
        #loop over levels
        for l in range(self.N):
            for v in range(level[l]['nv']-1):
                strT[t,0,:] = level[l]['v'][v]
                strT[t,1,:] = [level[l]['v'][v][0]-h,
                               level[l]['v'][v][1]+h,
                               level[l]['v'][v][2]]
                strT[t,2,:] = [level[l]['v'][v][0]-h,
                               level[l]['v'][v][1],
                               level[l]['v'][v][2]+h]
                t += 1
            for v in range(1,level[l]['nv']-1):
                strT[t,0,:] = level[l]['v'][v]
                strT[t,1,:] = [level[l]['v'][v][0]-h,
                               level[l]['v'][v][1],
                               level[l]['v'][v][2]+h]
                strT[t,2,:] = [level[l]['v'][v][0],
                               level[l]['v'][v][1]-h,
                               level[l]['v'][v][2]+h]
                t += 1
        #compute straight triangle centroids
        centroid = np.empty((self.nnodes,3), dtype='double')
        for t in range(self.nnodes):
            for c in range(3):
                centroid[t,c] = (strT[t,0,c]+strT[t,1,c]+strT[t,2,c])/3

        #project straight triangles onto spherical triangles
        #also computes the quadrature nodes (first octant)
        sphT = np.empty((self.nnodes,3,3), dtype='double')
        for t in range(self.nnodes):
            for p in range(3):
                factor = np.sqrt(strT[t,p,0]**2 +
                                 strT[t,p,1]**2 +
                                 strT[t,p,2]**2)
                for c in range(3):
                    sphT[t,p,c] = strT[t,p,c]/factor
        
            factor = np.sqrt(centroid[t,0]**2 +
                             centroid[t,1]**2 +
                             centroid[t,2]**2)
            self.mu[t] = centroid[t,0]/factor
            self.eta[t] = centroid[t,1]/factor
            self.xi[t] = centroid[t,2]/factor

        #compute the quadrature weights (first octant), i.e.
        #the surface of each spherical triangle
        for t in range(self.nnodes):
            #compute arc lengths
            a = 0
            b = 0
            c = 0
            for cc in range(3):
                a += sphT[t,1,cc]*sphT[t,2,cc]
                b += sphT[t,0,cc]*sphT[t,2,cc]
                c += sphT[t,0,cc]*sphT[t,1,cc]
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
                        
