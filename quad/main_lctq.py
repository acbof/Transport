#!/usr/bin/python3
"""
Example use of the Legendre-Chebyshev quadrangular quadrature set.

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

"""

import lctq as quad

#quadrature order
N = 6


#initiate the quadrature
q = quad.Lctq(N)
print(q)

#build quadrature set on first octant
q.buildFirstOctantSet()

#print quadrature set on first octant
q.printFirstOctantSet()

#build the complete quadrature set
q.build3d()

#diagnotics
q.diagnostics()

#plot
q.plotFirstOctant(show=True)

#build quadrature set on first quadrant
q.buildFirstQuadrantSet()

#print quadrature set on first quadrant
q.printFirstQuadrantSet()

#build the complete quadrature set (2D)
q.build2d()

#diagnotics
q.diagnostics()

#print quadrature set on radius
q.buildRadiusSet()

#print quadrature set on radius
q.printRadiusSet()

#build the complete quadrature set (1D)
q.build1d()

#diagnotics
q.diagnostics()