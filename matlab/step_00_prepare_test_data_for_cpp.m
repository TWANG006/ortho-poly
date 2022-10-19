clear;
close all;
clc;

addpath('functions/');
data_dir = '../data/';
legendre_file = 'square_surface_60mm.mat';


%% work with data for legendre
load([data_dir legendre_file]); % load the data

% downsample the data
n = 3;
X = X(1: n: end, 1: n: end);
Y = Y(1: n: end, 1: n: end);
Z = Z(1: n: end, 1: n: end);

% write the data
save([data_dir 'legendre_test_data.mat'], 'X', 'Y', 'Z');
write_matrix([data_dir 'X_legendre.bin'], X, 'double');
write_matrix([data_dir 'Y_legendre.bin'], Y, 'double');
write_matrix([data_dir 'Z_legendre.bin'], Z, 'double');

% plot the data
figure;
show_surface(X, Y, Z - nanmin(Z(:)), 1e9, 'nm', 'Test surface for Legendre');