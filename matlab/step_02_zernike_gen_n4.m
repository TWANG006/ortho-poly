clear;
close all;
clc;

addpath('functions/');

[X, Y] = meshgrid(-1: 0.01: 1, -1: 0.01: 1);
J = 1: 16;

[~, ~, ~, Z3, ~, ~] = zernike_xy_jc(X, Y, J, ones(length(J), 1));


%% plots
figure;
colormap jet;
rows = sqrt(length(J));
cols = 2 * rows - 1;

id = sqrt(X.^2 + Y.^2) <= 1;

for i = 1: sqrt(length(J))
    curr_col = 2*i - 1;
    cen = i^2;
    cen_id = (i - 1) * cols + (cols + 1) / 2;
    
    subplot(rows, cols, cen_id);
    Z = Z3(:, :, cen);
    Z(~id) = NaN;
    imagesc(Z); axis image xy; set(gca, 'xtick', [], 'ytick', []);
    title(['Q' num2str(cen)]);
    
    for j = 1: (curr_col - 1)/2
        subplot(rows, cols, cen_id + j);
        Z = Z3(:, :, cen - 2 * j);
        Z(~id) = NaN;
        imagesc(Z); axis image xy; set(gca, 'xtick', [], 'ytick', []);
        title(['Q' num2str(cen - 2 * j)]);
        
        subplot(rows, cols, cen_id - j);
        Z = Z3(:, :, cen - 2 * j + 1);
        Z(~id) = NaN;
        imagesc(Z); axis image xy; set(gca, 'xtick', [], 'ytick', []);
        title(['Q' num2str(cen - 2 * j + 1)]);
    end    
end