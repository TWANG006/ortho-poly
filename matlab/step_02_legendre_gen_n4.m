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

% figure;
% colormap jet;
% subplot(5, 9, 5);
% imagesc(Z3(:, :, 1)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q1: Piston');
% 
% subplot(5, 9, 13);
% imagesc(Z3(:, :, 2)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q2: x-tilt');
% subplot(5, 9, 15);
% imagesc(Z3(:, :, 3)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q3: y-tilt');
% 
% subplot(5, 9, 21);
% imagesc(Z3(:, :, 4)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q4: x-defocus');
% subplot(5, 9, 23);
% imagesc(Z3(:, :, 5)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q5');
% subplot(5, 9, 25);
% imagesc(Z3(:, :, 6)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q6: y-defocus');
% 
% subplot(5, 9, 29);
% imagesc(Z3(:, :, 7)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q7: x-primary coma');
% subplot(5, 9, 31);
% imagesc(Z3(:, :, 8)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q8');
% subplot(5, 9, 33);
% imagesc(Z3(:, :, 9)); axis image; set(gca, 'xtick', [], 'ytick', []);title('Q9');
% subplot(5, 9, 35);
% imagesc(Z3(:, :, 10)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q10: y-primary coma');
% 
% 
% subplot(5, 9, 37);
% imagesc(Z3(:, :, 11)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q11: x-primary spherical');
% subplot(5, 9, 39);
% imagesc(Z3(:, :, 12)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q12');
% subplot(5, 9, 41);
% imagesc(Z3(:, :, 13)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q13');
% subplot(5, 9, 43);
% imagesc(Z3(:, :, 14)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q14');
% subplot(5, 9, 45);
% imagesc(Z3(:, :, 15)); axis image; set(gca, 'xtick', [], 'ytick', []); title('Q15: y-primary spherical');