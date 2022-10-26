clear;
close all;
clc;

addpath('functions/');

[X, Y] = meshgrid([-1 0 1], [-1 0 1]);
Y = -Y;
X
Y

[Z, Zx, Zy, z3, zx3, zy3, m, n]  = legendre_xy_jc(X, Y, 1: 3, [1,1]);
z3