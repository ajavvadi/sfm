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
## @deftypefn {Function File} {@var{retval} =} SfM2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Kalyan <kalyan@kalyan-XPS-L501X>
## Created: 2016-05-22

function SfM2 (path_im1, path_im2, dir_path)

%clear all;
%clc;

vl_setup;

im1 = imread(path_im1);
im2 = imread(path_im2);

tot_mat = dlmread(dir_path);
i_calib = tot_mat(1:3, 1:3)

%im1 = imread('/home/kalyan/Downloads/rectify-tutorial/desk1.png');
%im2 = imread('/home/kalyan/Downloads/rectify-tutorial/desk2.png');
%im1 = imread('/home/kalyan/Desktop/assignment4/images_folder/fountain/1.jpg');
%im2 = imread('/home/kalyan/Desktop/assignment4/images_folder/fountain/2.jpg');

%i_calib = [720 0 360; 0 720 288; 0 0 1];

[r1 c1 t1] = size(im1);
[r2 c2 t2] = size(im2);

i1_or = im1;
i2_or = im2;

if(t1 > 1)
  im1 = single(rgb2gray(im1));
endif

if(t2 > 1)
  im2 = single(rgb2gray(im2));
endif

[f_i1 d_i1] = vl_sift(im1);
[f_i2 d_i2] = vl_sift(im2);

d_i1 = double(d_i1);
d_i2 = double(d_i2);

n_d_i1 = diag(sqrt((d_i1')*d_i1));
n_d_i2 = diag(sqrt((d_i2')*d_i2));

for n_iterator = 1:size(d_i1, 2)
  norm_d_i1(:, n_iterator) = d_i1(:, n_iterator)./n_d_i1(n_iterator);
endfor

for n_iterator = 1:size(d_i2, 2)
  norm_d_i2(:, n_iterator) = d_i2(:, n_iterator)./n_d_i2(n_iterator);
endfor

cp_mat = zeros(size(max(size(d_i1, 2), size(d_i2, 2))), 2);

%[matches, scores] = vl_ubcmatch(d_i1, d_i2);

num_match_p = 1;
for d_iterator = 1:size(d_i1, 2)

  sym_test_mat = (norm_d_i1(:, d_iterator)' * norm_d_i2);
  [sorted_mat, indices_mat] = sort(acos(sym_test_mat), 'ascend');
  ratio_cp = sorted_mat(1)/sorted_mat(2);
  if (ratio_cp < 0.999)
    cp_mat(num_match_p, 1) = d_iterator;
    cp_mat(num_match_p, 2) = indices_mat(1);
    num_match_p = num_match_p + 1;
  else
  
  end

endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_prime - first image xcord
%x - second image xcord
%y_prime - first image ycord
%y - secnd image ycord
%Left Image - I1
%Right Image - I2
%%%RANSAC%%%%%%%%%%%%%%%%%%%%%

iter = 200;
p_to_c = min(size(cp_mat, 1), 10);
past_error = 100000;
past_n_inliers = 0;
f_ransac = zeros(3, 3);
%inliers_ransac = zeros(p_to_c, );
k = 20000000;

for i = 1:iter

  a = 1;
  b = size(cp_mat, 1);
  rand_ind = (b - a).*rand(p_to_c, 1) + a;
  rand_ind = round(rand_ind);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %Computing fundamental matrix
  %Destination: given by Right
  %Source: Left
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  x = f_i1(1, cp_mat(rand_ind, 1))';
  y = f_i1(2, cp_mat(rand_ind, 1))';
  
  x_prime = f_i2(1, cp_mat(rand_ind, 2))';
  y_prime = f_i2(2, cp_mat(rand_ind, 2))';

  f_mat = compute_f_matrix(x_prime, y_prime, x, y);
  
  inlier_ind = 1;
  error_term = 0;
  
  inlier_mat = zeros(size(cp_mat, 1), 1);
  counter = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %finding inliers
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  for j = 1:b
  
  x_cord_r = f_i2(1, cp_mat(j, 2));
  y_cord_r = f_i2(2, cp_mat(j, 2));
  
  x_cord_l = f_i1(1, cp_mat(j, 1));
  y_cord_l = f_i1(2, cp_mat(j, 1));
  
  excess_term = [x_cord_r y_cord_r 1] * f_mat * [x_cord_l y_cord_l 1]';
  
  %F_x1 = f_mat*
  
  if(excess_term < 0.006)
  %if(sqrt((floor(x_cord_d) - floor(x_cord_dash))^2 + (floor(y_cord_d) - floor(y_cord_dash))^2) < 36)
    inlier_mat(j) = 1;
    inlier_ind = inlier_ind + 1;
  endif
  
  endfor
  
  error_term;
  
  if(past_n_inliers < inlier_ind)
    past_n_inliers = inlier_ind;
    f_ransac = f_mat;
    if(inlier_ind > 1)
      inliers_ransac = inlier_mat;
      
    endif
  endif

endfor

f_ransac = f_ransac/f_ransac(3, 3)

non_zero_inliers = find(inliers_ransac);
x = f_i1(1, cp_mat(non_zero_inliers, 1))';
y = f_i1(2, cp_mat(non_zero_inliers, 1))';
x_prime = f_i2(1, cp_mat(non_zero_inliers, 2))';
y_prime = f_i2(2, cp_mat(non_zero_inliers, 2))';

f_mat_final = compute_f_matrix(x_prime, y_prime, x, y);

f_mat_final = f_mat_final./f_mat_final(3, 3);

e_mat = i_calib' * f_mat_final * i_calib;

[u_e d_e v_e] = svd(e_mat);
d_e(2, 2) = d_e(1, 1);
d_e(3, 3) = 0;
%e_mat_final = u_e * d_e * v_e';
e_mat_final = e_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute T and R from Essesntial matrix
%% This requires triagulation also!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimating the T and R
W_e = [0 -1 0; 1 0 0; 0 0 1];
Z_e = [0 1 0; -1 0 0; 0 0 0];
[U_e D_e V_e] = svd(e_mat_final);

t_e = U_e(:, 3)./U_e(3, 3);

%U_e = U_e * division_fact;
%V_e = V_e * division_fact;
%C(:, :, 1) = -U_e*Z*U_e';
%C(:, :, 2) = U_e*Z*U_e';
%C(:, :, 3) = U_e*W'*V_e';
C(:, :, 4) = U_e*W_e*V_e';
U_e(:, 3)

P1 = i_calib * [eye(3) zeros(3, 1)];
P2_big_mat = zeros(3, 4, 4);
%P2_big_mat(:, :, 1) = i_calib * [C(:, :, 4) U_e(:, 3)];
%P2_big_mat(:, :, 2) = i_calib * [C(:, :, 4) -U_e(:, 3)];
%P2_big_mat(:, :, 3) = i_calib * [U_e*W_e'*V_e' U_e(:, 3)];
%P2_big_mat(:, :, 4) = i_calib * [U_e*W_e'*V_e' -U_e(:, 3)];

P2_big_mat(:, :, 1) = i_calib * [C(:, :, 4) t_e];
P2_big_mat(:, :, 2) = i_calib * [C(:, :, 4) -t_e];
P2_big_mat(:, :, 3) = i_calib * [U_e*W_e'*V_e' t_e];
P2_big_mat(:, :, 4) = i_calib * [U_e*W_e'*V_e' -t_e];

pinv(P1)
for r_t_iterator = 1:4
  P2 = P2_big_mat(:, :, r_t_iterator);
  points_1 = [x';y'];
  points_2 = [x_prime'; y_prime'];
  
  id = fopen(strcat("data", num2str(r_t_iterator), ".ply"), "w");
  fprintf(id, 'ply\n');
  fprintf(id, 'format ascii 1.0\n');
  fprintf(id, cstrcat("element vertex", " ", num2str(size(points_1, 2)), "\n"));
  fprintf(id, 'property float x\n');
  fprintf(id, 'property float y\n');
  fprintf(id, 'property float z\n');
  fprintf(id, 'property uchar red\n');
  fprintf(id, 'property uchar green\n');
  fprintf(id, 'property uchar blue\n');
  fprintf(id, 'end_header\n');
  
  
  for p_iterator = 1:size(points_1, 2)
    yosh(p_iterator, :, r_t_iterator) = triangulate_test(P1, P2, points_1(:, p_iterator), points_2(:, p_iterator));
    r = i1_or(round(y(p_iterator)), round(x(p_iterator)), 1);
    g = i1_or(round(y(p_iterator)), round(x(p_iterator)), 2);
    b = i1_or(round(y(p_iterator)), round(x(p_iterator)), 3);
    fprintf(id, '%0.2f %0.2f %0.2f %d %d %d\n', yosh(p_iterator, 1, r_t_iterator), yosh(p_iterator, 2, r_t_iterator), yosh(p_iterator, 3, r_t_iterator), r, g, b);
  endfor
  
  fclose(id);
endfor

endfunction
