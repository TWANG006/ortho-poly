function [Z, Zx, Zy, z3, zx3, zy3, m, n] = legendre_xy_jc(X, Y, J, C)
%LEGENDREXYJC Legendre polynominals in normalized X and Y coordinates 
%   with order vector J and coefficient vector C.

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2020-01-17 Original Version
% Modified by Tianyi on 10/28/2022 for the Q ordering

% Release the order vector J as a column vector............................
J = J(:);

% Calculate the order vectors (n, m).......................................
% b = ceil(sqrt(J));
% a = b.^2-J+1;
% 
% nsm = -a/2.*(~mod(a,2))+(a-1)/2.*(mod(a,2));
% nam = 2*b-abs(nsm)-2;
% 
% n = (nam+nsm)/2;
% m = (nam-nsm)/2;

% Qj is in: row a & col b
a = floor((1 + sqrt(1 + 8*(J-1)))*0.5);
b = (J-1) - 0.5*(a.*(a-1)) + 1;

% Get n, m from a & b
n= a - b;
m = b - 1;

% Calculate the Legendre polynominals with thier coefficients, as well as
% their 1st derivatives....................................................
[z3, zx3, zy3] = legendre_xy_nm(X, Y, n, m);

for k = 1:length(C)
    z3 (:,:,k) = z3 (:,:,k).*C(k);
    zx3(:,:,k) = zx3(:,:,k).*C(k);
    zy3(:,:,k) = zy3(:,:,k).*C(k);
end

% Sum up the modes.........................................................
Z  = sum(z3,3);
Zx = sum(zx3,3);
Zy = sum(zy3,3);





