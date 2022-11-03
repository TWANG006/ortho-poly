clear;
close all;
clc;

addpath('functions/');

[X, Y] = meshgrid(-1: 0.01: 1, -1: 0.01: 1);
J = 1: 16;

[~, ~, ~, Z3, ~, ~] = legendre_xy_jc(X, Y, J, ones(size(J)));

b = ceil(sqrt(J));
a = b.^2-J+1;

nsm = -a/2.*(~mod(a,2))+(a-1)/2.*(mod(a,2));
nam = 2*b-abs(nsm)-2;


%% plotting
figure;
rows = max(a);
cols = max(b);
colormap jet;
for i = 1: 15
    subplot(cols, rows, rows * (a(i) - 1) + b(i));
    imagesc(Z3(:, :, i)); axis image xy; set(gca, 'xtick', [], 'ytick', []);
    title(['Q' num2str(i)]);
end