clear;
close all;
clc;

[X, Y] = meshgrid(-1: 1, -1: 1);
Y = -Y;

[R, Theta] = cart2pol(X, Y);

[Z, ~, ~, z3, ~, ~] = zernike_xy_jc(X, Y, 1:4, [1 1 1 1]);