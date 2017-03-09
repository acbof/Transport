"""
Legendre-Chebyshev triangular quadrature
________________________________________

This Python code provides points and weights of the Legendre-Chebyshev 
triangular quadrature for integrals on the unit sphere, i.e.:

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
from numpy import polynomial as poly

print("Legendre-Chebyshev triangular quadrature\n\
Author: Pedro H A Konzen - UFRGS - 2017\n\n\
This program comes with ABSOLUTELY NO WARRANTY.\n\n\
This program is free software: you can redistribute it and/or modify\n\
it under the terms of the GNU General Public License as published by\n\
the Free Software Foundation, either version 3 of the License, or\n\
(at your option) any later version.\n\n\
")

#order (must be even)
N = 4

#check parity
if (N % 2):
    raise Exception('N must be even!')

#Gauss-Legendre quadrature (x,w) - (points, weights)
x,w = poly.legendre.leggauss(N);
#just the positive values
xi = x[N/2:] 
w0 = w[N/2:]

#wbar
wbar = np.zeros((N,N),dtype='double')
for i in np.arange(1,N/2+1):
    for j in np.arange(1,N-2*i+2+1):
        wbar[i-1,j-1] = np.pi/2 * (1 - (N-2*j-2*i+3)/(N-2*i+2))

#Quadratute points on the first octant (3D computations)
M = N*(N+2)/8
mu3= np.zeros(M,dtype='double')
eta3= np.zeros(M,dtype='double')
xi3= np.zeros(M,dtype='double')
ww3= np.zeros(M,dtype='double')

c = 0
for i in np.arange(1,N/2+1):
    for j in np.arange(1,N-2*i+2+1):
        aux = np.sqrt(1-xi[i-1]**2)*np.cos(wbar[i-1,j-1])
        if (aux > 0.0):
            mu3[c] = aux
            eta3[c] = np.sqrt(1-mu3[c]**2-xi[i-1]**2)
            xi3[c] = np.sqrt(1 - mu3[c]**2 - eta3[c]**2)
            ww3[c] = np.pi * w0[i-1]/(N-2*i+2)
            c += 1

#print quadrature table
print("Quadrature table for 3D computations")
print("[mu, eta, xi, w]")
for i in np.arange(M):
    print(mu3[i],eta3[i],xi3[i],ww3[i])

#Testing
#build the complete table
n = 8*M
mu = np.zeros(n, dtype="double")
eta = np.zeros(n, dtype="double")
xi = np.zeros(n, dtype="double")
ww = np.zeros(n, dtype="double")

c=0
for m in np.arange(M):
    for i in [1,-1]:
        for j in [1,-1]:
            for k in [1,-1]:
                mu[c] = i*mu3[m]
                eta[c] = j*eta3[m]
                xi[c] = k*xi3[m]
                ww[c] = ww3[m]
                c += 1

#Zeroth moment
print("\n\n Zeroth moment:")
e = 4*np.pi
print("\tExpected: %1.5e" % e)
s = np.sum(ww)
print("\tComputed: %1.5e" % s)
print("\tError: %1.5e" % np.fabs(e-s))

#First moment
print("\n\n First moment:")
e = np.zeros(3)
print("\tExpected: ", e)
s = np.zeros(3)
for i in np.arange(n):
    s += ww[i]*np.array([mu[i],eta[i],xi[i]])
print("\tComputed: ", s)
print("\tError: " , np.fabs(e-s))

#Second moment
print("\n\n Second moment:")
e = 4*np.pi/3*np.eye(3)
print("\tExpected: \n", e)
s = 0
for i in np.arange(n):
    s += ww[i] * np.multiply(
        np.array([[mu[i]],[eta[i]],[xi[i]]]),
        np.array([mu[i],eta[i],xi[i]]))
print("\tComputed: \n", s)
print("\tError: \n" , np.fabs(s - e))


"""
*********************************************
"""

#print quadrature table
print("\n\nQuadrature table for 2D computations")
mu2 = np.empty(1)
eta2 = np.empty(1)
ww2 = np.empty(1)

c=0
for m in np.arange(M):
    not_there = True
    for i in np.arange(0,c):
        if ((np.fabs(mu3[m] - mu2[i]) < 1e-15) and
            (np.fabs(eta3[m] - eta2[i]) < 1e-15)):
            not_there = False
            ww2[i] += ww3[m]
            break
    if ((c==0) or (not_there)):
        if (c==0):
            mu2[c] = mu3[m]
            eta2[c] = eta3[m]
            ww2[c] = 2*ww3[m]
        else:
            mu2 = np.append(mu2,mu3[m])
            eta2 = np.append(eta2,eta3[m])
            ww2 = np.append(ww2,2*ww3[m])
        c += 1

M2 = c

print("[mu, eta, w]")
for i in np.arange(M2):
    print(mu2[i],eta2[i],ww2[i])

#Testing
#build the complete table
n = 4*M2
mu = np.zeros(n, dtype="double")
eta = np.zeros(n, dtype="double")
ww = np.zeros(n, dtype="double")

c=0
for m in np.arange(M2):
    for i in [1,-1]:
        for j in [1,-1]:
            mu[c] = i*mu2[m]
            eta[c] = j*eta2[m]
            ww[c] = ww2[m]
            c += 1

#Zeroth moment
print("\n\n Zeroth moment:")
e = 4*np.pi
print("\tExpected: %1.5e" % e)
s = np.sum(ww)
print("\tComputed: %1.5e" % s)
print("\tError: %1.5e" % np.fabs(e-s))

#First moment
print("\n\n First moment:")
e = np.zeros(2)
print("\tExpected: ", e)
s = np.zeros(2)
for i in np.arange(n):
    s += ww[i]*np.array([mu[i],eta[i]])
print("\tComputed: ", s)
print("\tError: " , np.fabs(e-s))

#Second moment
print("\n\n Second moment:")
e = 4*np.pi/3*np.eye(2)
print("\tExpected: \n", e)
s = 0
for i in np.arange(n):
    s += ww[i] * np.multiply(
        np.array([[mu[i]],[eta[i]]]),
        np.array([mu[i],eta[i]]))
print("\tComputed: \n", s)
print("\tError: \n" , np.fabs(s - e))


"""
*********************************************
"""

#print quadrature table
print("\n\nQuadrature table for 1D computations")
mu1 = np.empty(1)
ww1 = np.empty(1)

c=0
for m in np.arange(M):
    not_there = True
    for i in np.arange(0,c):
        if (np.fabs(mu3[m] - mu1[i]) < 1e-15):
            not_there = False
            ww1[i] += ww3[m]
            break
    if ((c==0) or (not_there)):
        if (c==0):
            mu1[c] = mu3[m]
            ww1[c] = 4*ww3[m]
        else:
            mu1 = np.append(mu1,mu3[m])
            ww1 = np.append(ww1,4*ww3[m])
        c += 1

M1 = c

print("[mu, w]")
for i in np.arange(M1):
    print(mu1[i],ww1[i])

#Testing
#build the complete table
n = 2*M1
mu = np.zeros(n, dtype="double")
ww = np.zeros(n, dtype="double")

c=0
for m in np.arange(M2):
    for i in [1,-1]:
        mu[c] = i*mu1[m]
        ww[c] = ww1[m]
        c += 1

#Zeroth moment
print("\n\n Zeroth moment:")
e = 4*np.pi
print("\tExpected: %1.5e" % e)
s = np.sum(ww)
print("\tComputed: %1.5e" % s)
print("\tError: %1.5e" % np.fabs(e-s))

#First moment
print("\n\n First moment:")
e = 0.0
print("\tExpected: ", e)
s = 0.0
for i in np.arange(n):
    s += ww[i]*mu[i]
print("\tComputed: ", s)
print("\tError: " , np.fabs(e-s))

#Second moment
print("\n\n Second moment:")
e = 4*np.pi/3
print("\tExpected: \n", e)
s = 0
for i in np.arange(n):
    s += ww[i] * mu[i]**2
print("\tComputed: \n", s)
print("\tError: \n" , np.fabs(s - e))

        
        





        











