clear;
close all;
clc;

data_dir = '../data/';
surf_file = 'legendre_test_data.mat';
load([data_dir surf_file]);


%% normalization test
% matlab side
X_lg_m = -1 + 2.*(X - min(X(:)))./(max(X(:)) - min(X(:)));
Y_lg_m = -1 + 2.*(Y - min(Y(:)))./(max(Y(:)) - min(Y(:)));

% c++ side
X_lg_c = load_matrix([data_dir 'X_lg.bin'], 'double');
Y_lg_c = load_matrix([data_dir 'Y_lg.bin'], 'double');

figure;
subplot(1, 2, 1);
scatter(X_lg_m(:), Y_lg_m(:));
subplot(1, 2, 2);
scatter(X_lg_c(:), Y_lg_c(:));