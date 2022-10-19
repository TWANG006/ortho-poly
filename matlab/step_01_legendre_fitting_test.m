clear;
close all;
clc;

data_dir = '../data/';
surf_file = 'legendre_test_data.mat';
load([data_dir surf_file]);


%% normalization test
% matlab side
X_nor_m = -1 + 2.*(X - min(X(:)))./(max(X(:)) - min(X(:)));
Y_nor_m = -1 + 2.*(Y - min(Y(:)))./(max(Y(:)) - min(Y(:)));

% c++ side
