function [Z, Zx, Zy, z3, zx3, zy3] = zernike_xy_jc(X, Y, J, C, order)
%ZERNIKEXYJC Zernike polynominals in normalized X and Y coordinates
%   with order vector J and coefficient vector C.

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-11-02 Original Version

% Release order vector J and coefficient vector C..........................
C = C(J~=0);
J = J(J~=0);

% Get the order corresonding to the function number........................
% Check the input arguments to see if user defined the order type.
if(nargin < 5)
    order = 'Carl Zeiss'; % by default.
else
    if(isempty(order))
        order = 'Carl Zeiss'; % if empty defination.
    end
end

% Check which zernike order is used.
if(strcmpi(order,'Carl Zeiss'))
    J = J(:);
    b = ceil(sqrt(J));
    a = b.^2-J+1;
    m = -a/2.*(~mod(a,2))+(a-1)/2.*(mod(a,2));
    n = 2*(b-1)-abs(m);
elseif(strcmpi(order,'Original'))
    
    if any(J>35)
        error(['Only the first 36 ''Original'' ' ...
            'Zernike functions are computed (Ind = 0~35).'])
    end
    p = J(:)-1;
    n = ceil((-3+sqrt(9+8*p))/2);
    m = 2*p - n.*(n+2);
else
    error('CreateAberrWithZernike:order',['Unknown order = ',order]);
end

% Zernikes in X and Y......................................................
[z3,zx3,zy3] = zernike_xy_nm(X,Y,n,m);
for k = 1:length(C)
    z3 (:,:,k) = z3 (:,:,k).*C(k);
    zx3(:,:,k) = zx3(:,:,k).*C(k);
    zy3(:,:,k) = zy3(:,:,k).*C(k);
end
Z  = sum(z3,3);
Zx = sum(zx3,3);
Zy = sum(zy3,3);

end % END of Function
