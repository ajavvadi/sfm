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
## @deftypefn {Function File} {@var{retval} =} triangulate_test (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Kalyan <kalyan@kalyan-XPS-L501X>
## Created: 2016-05-22

function coords = triangulate_test (P1, P2, m1, m2)

  A = [m1(1)*P1(3, :)-P1(1, :); m1(2)*P1(3, :)-P1(2, :); m2(1)*P2(3, :) - P2(1, :); m2(2)*P2(3, :)-P2(2, :)];
  
  [u_a s_a v_a] = svd(A);
  
  coords = v_a(:, end);
  coords = coords(1:3)/coords(4);

endfunction
