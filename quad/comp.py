# -*- coding: utf-8 -*-
"""
Comparative analisys for precision of quadrature sets on the unit sphere.

This Python3 code is a inherited class of quadrature sets on the unit
sphere, i.e.

..math: 
   \\int_{S^2} f(\\Omega) d\\Omega \\approx \\sum_{i=1}^M w_i f(\\Omega_i),

where 

..math: 
   \\Omega = (\\mu, \\eta, \\xi)\in S^2 := \\{\\mu^2+\\eta^2+\\xi^2=1\\}.

Author
______

Ana Carolina Bof - UFRGS - Jul/2018

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

"""


import numpy as np
from scipy import integrate
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lcqq
import lctq
import srap
import tesselation

note = "Author:\n" +\
       "\tAna Carolina Bof - UFRGS - 2018\n\n" +\
       "This program comes with ABSOLUTELY NO WARRANTY.\n\n" +\
       "This program is free software: you can redistribute" +\
       "it and/or modify\n" +\
       "it under the terms of the GNU General Public" +\
       "License as published by\n" +\
       "the Free Software Foundation, either version 3" +\
       "of the License, or (at\n" +\
       "your option) any later version.\n\n"

print(note)

N = 4
N = int(N)

TOL = 1e-7

assert (N%2 == 0), "N must be even."

lcq = lcqq.Lcqq(N)
lcq.build1d()
lct = lctq.Lctq(N)
lct.build1d()
sra = srap.SRAP(N)
sra.build1d()
tess = tesselation.Tesselation(N)
tess.build1d()


print('\nN = %i\n' % (N))

vgauss = []
verr1 = []
verr2 = []
verr3 = []
verr4 = []

nt=50
print("s              LCQQ               LCTQ               SRAP               Tesselation        Exact value")
for n in np.arange(nt):            
    vesp=2*np.pi/(n+1)
    verr1.append(lcq.nthMomentError(1,n))
        
    verr2.append(lct.nthMomentError(1,n))

    verr3.append(sra.nthMomentError(1,n))

    verr4.append(tess.nthMomentError(1,n))

    vgauss.append(2*np.pi*integrate.fixed_quad(lambda x: x**n,0,1,n=N)[0])
    vgauss[n] = np.fabs(vgauss[n]-vesp)/np.fabs(vesp)

    print("%d            %1.2e         %1.2e         %1.2e         %1.2e         %1.2e         %1.2e" % \
         (n,\
          verr1[n],\
          verr2[n],\
          verr3[n],\
          verr4[n],\
          vgauss[n],\
          vesp))

plt.plot(np.arange(nt),verr1,color="red",linestyle="-",label="LCQQ")
plt.plot(np.arange(nt),verr2,color="red",linestyle="--",label="LCTQ")
plt.plot(np.arange(nt),verr3,color="blue",linestyle="-",label="SRAP")
plt.plot(np.arange(nt),verr4,color="blue",linestyle="--",label="TESS")
plt.plot(np.arange(nt),vgauss,color="green",linestyle="--",label="TESS")
plt.legend()
plt.yscale('log')
plt.show()
