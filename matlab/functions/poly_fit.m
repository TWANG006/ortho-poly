function [Zfit, c] = poly_fit(X, Y, Z, j, type)

% normalization
X_nor = -1 + 2.*(X - min(X(:)))./(max(X(:)) - min(X(:)));
Y_nor = -1 + 2.*(Y - min(Y(:)))./(max(Y(:)) - min(Y(:)));


if(strcmp(type,'zernike'))
    [~, ~, ~, Z3, ~, ~] = zernike_xy_jc(X_nor, Y_nor, j, ones(size(j)));
elseif(strcmp(type,'legendre'))
    [~, ~, ~, Z3, ~, ~] = legendre_xy_jc(X_nor, Y_nor, j, ones(size(j)));
else
    error('Unkown polynomial type.');
end

z3_res = reshape(Z3, [],size(Z3,3));

A = z3_res(~isnan(Z(:)),:);
b = Z(~isnan(Z(:)));

c = A\b;

for i = 1:length(c)
    Z3(:,:,i) = Z3(:,:,i)*c(i);
end

Zfit = sum(Z3,3);


end