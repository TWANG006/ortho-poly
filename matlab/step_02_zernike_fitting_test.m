clear;
close all;
clc;

%% test fit
load('../data/zernike_test_data.mat');
[Zfit, c] = poly_fit(X, Y, Z, 1: 16, 'zernike');
Zfit_c = load_matrix('../data/Zfit_zernike.bin', 'double');
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
show_surface(X, Y, Zfit_c, 1e9, 'nm', 'Fit with Q1 ~ Q16');
subplot(1, 3, 3);
show_surface(X, Y, Zfit - Z, 1e9, 'nm', 'Difference');