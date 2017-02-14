## Copyright (C) 2016 Kalyan
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} compute_f_matrix (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Kalyan <kalyan@kalyan-XPS-L501X>
## Created: 2016-05-21

function F = compute_f_matrix (x_prime, y_prime, x, y)

  [coords_prime coords T_prime T] = compute_dlt(x_prime, y_prime, x, y);
  
  x_prime_n = coords_prime(1, :)';
  y_prime_n = coords_prime(2, :)';
  x_n = coords(1, :)';
  y_n = coords(2, :)';
  
  A = [(x_prime_n.*x_n) (x_prime_n.*y_n) x_prime_n (y_prime_n.*x_n) (y_prime_n.*y_n) y_prime_n x_n y_n ones(length(y_prime_n), 1)];
  
  [U, S, V] = svd(A);
  f = V(:, max(size(V)));
  Fp = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
  [Uf Sf Vf] = svd(Fp);
  Sf(3, 3) = 0;
  Fprime = Uf * Sf * Vf';
  
  F = T_prime'*Fprime*T;

endfunction
