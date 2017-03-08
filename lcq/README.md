#lcq - Legendre-Chebyshev quadratures

This folder contains some useful codes for dealing with Legendre-Chebyshev 
quadratures for integrals over the unit sphere, i.e.:

..math  ::
  \\int_{S^2} f(\\Omega) d\\Omega \\approx \\sum_{i=1}^M w_i f(\\Omega_i)$,

where 

..math :: \\Omega = (\\mu, \\eta, \\xi)\in S^2 := \\{\\mu^2+\\eta^2+\\xi^2=1\\}$.

Author
______

Pedro H A Konzen - UFRGS - Mar/2017

Contents
________

* lcqq - Legendre-Chebyshev quadrangular quadrature
* lctq - Legendre-Chebyshev triangular quadrature

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
