'''
Gauss-Chebyshev quadrangular quadrature

Pedro H A Konzen - UFRGS - Mar/2017
'''

import numpy as np
from numpy import polynomial as poly

#order (must be even)
N = 4

#check eveness
if (N % 2):
    raise Exception('N must be even!')

#Gauss-Legendre quadrature (x,w) - (points, weights)
x,w = poly.legendre.leggauss(N);
xi = x[N/2:] #just the positive values
w0 = w[N/2:]
#w = w/np.sum(w)

#wbar?
wbar = np.zeros(N,dtype='double')
for j in np.arange(1,N+1):
    wbar[j-1] = np.pi/2 * (1 - (N-2*j+1)/N)

#mu
mu0= np.zeros(N**2/4,dtype='double')
eta0= np.zeros(N**2/4,dtype='double')
xi0= np.zeros(N**2/4,dtype='double')
ww0= np.zeros(N**2/4,dtype='double')

for i in np.arange(N/2):
    for j in np.arange(N/2):
        mu0[i+j*N/2] = np.sqrt(1-xi[i]**2)*np.cos(wbar[j])
        eta0[i+j*N/2] = np.sqrt(1-mu0[i+j*N/2]**2-xi[i]**2)
        xi0[i+j*N/2] = np.sqrt(1 - mu0[i+j*N/2]**2 - eta0[i+j*N/2]**2)
        ww0[i+j*N/2] = np.pi * w[i]/N

#print table
print("[mu, eta, xi, w]")
for i in np.arange(N**2/4):
    print(mu0[i],eta0[i],xi0[i],ww0[0])

#build complete 2d table
n = 4*np.size(mu0)
mu = np.zeros(n)
eta = np.zeros(n)
ww = np.zeros(n)

c = 0
for i in [1,-1]:
    for j in [1,-1]:
        for k in np.arange(n/4):
            mu[c] = i*mu0[k]
            eta[c] = j*eta0[k]
            ww[c] = 2*ww0[k]
            c += 1

#Zeroth moment
print("\n\n Zeroth moment:")
print("\tExpected: %1.5e" % (4*np.pi))
print("\tComputed: %1.5e" % (np.sum(ww)))
print("\tError: %1.5e" % (np.fabs(4*np.pi - np.sum(ww))))

#First moment
s = np.zeros(2)
for i in np.arange(n):
    s += ww[i]*np.array([mu[i],eta[i]])
print("\n\n First moment:")
print("\tExpected: ", np.zeros(2))
print("\tComputed: ", s)
print("\tError: " , np.fabs(s - np.zeros(2)))

#Second moment
s = np.zeros((2,2))
for i in np.arange(n):
    s += ww[i]*np.multiply(np.array([[mu[i]],[eta[i]]]),np.array([[mu[i],eta[i]]]))
print("\n\n Second moment:")
e = 4*np.pi/3 * np.eye(2)
print("\tExpected: ", e)
print("\tComputed: ", s)
print("\tError: " , np.fabs(s - e))

            
        
        





        











