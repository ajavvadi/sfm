%clear all;
%clc;

im1_path = '/home/kalyan/Desktop/assignment4/images_folder/medusa/1.jpg';
im2_path = '/home/kalyan/Desktop/assignment4/images_folder/medusa/2.jpg';
dir_path = '/home/kalyan/Desktop/assignment4/images_folder/medusa/intrinsic.new';

SfM2(im1_path, im2_path, dir_path);