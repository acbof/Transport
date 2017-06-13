"""
Example use of the tesselation quadrature set.

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

import tesselation as quad

#quadrature order
N = 20

#initiate the quadrature
q = quad.Tesselation(N)
print(q)

#print quadrature set on first octant
q.printFirstOctantSet()

#build the complete quadrature set
q.build3d()

#diagnotics
q.diagnostics()
