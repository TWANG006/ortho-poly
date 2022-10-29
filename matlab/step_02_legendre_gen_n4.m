clear;
close all;
clc;

addpath('functions/');

[X, Y] = meshgrid(-1: 0.01: 1, -1: 0.01: 1);
J = 1: 15;

[~, ~, ~, Z3, ~, ~] = legendre_xy_jc(X, Y, [1: 15], ones(15, 1));
a = floor((1 + sqrt(1 + 8*(J-1)))*0.5);
b = (J-1) - 0.5*(a.*(a-1)) + 1;
figure;
rows = max(a);
cols = max(b);
colormap jet;
for i = 1: 15
    subplot(cols, rows, rows * (a(i) - 1) + b(i));
    imagesc(Z3(:, :, i)); axis image; set(gca, 'xtick', [], 'ytick', []);
    title(['Q' num2str(i)]);
end