clear;
close all;
clc;

addpath('functions/');

% %%
% [X, Y] = meshgrid([-1 0 1], [-1 0 1]);
% Y = -Y;
% 
% 
% [Z, Zx, Zy, z3, zx3, zy3, m, n]  = legendre_xy_jc(X, Y, 1: 3, [1,1]);
% Z
% z3

%% test fit
load('../data/square_surface_60mm.mat');
[Zfit, c] = poly_fit(X, Y, Z, 1: 15, 'legendre');
Zfit_c = load_matrix('../data/Zfit.bin', 'double');
Zdif = Zfit - Zfit_c;
Zdif = Zdif - nanmean(Zdif);

figure;
subplot(2, 2, 1);
show_surface(X, Y, Z, 1e9, 'nm', 'Original');
subplot(2, 2, 2);
show_surface(X, Y, Zfit, 1e9, 'nm', 'MATLAB');
subplot(2, 2, 3);
show_surface(X, Y, Zfit_c, 1e9, 'nm', 'C++');
subplot(2, 2, 4);
show_surface(X, Y, Zfit - Zfit_c, 1e9, 'nm', 'Difference');

figure;
subplot(1, 3, 1);
show_surface(X, Y, Z, 1e9, 'nm', 'Input surface');
subplot(1, 3, 2);
show_surface(X, Y, Zfit_c, 1e9, 'nm', 'Fit with Q1 ~ Q15');
subplot(1, 3, 3);
show_surface(X, Y, Zfit - Z, 1e9, 'nm', 'Difference');
