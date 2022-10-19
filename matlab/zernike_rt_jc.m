function Z = zernike_rt_jc(R,Theta,J,C,order)
%ZERNIKERTJC Zernike polynominals in normalized R and T coordinates
%   with order vector J and coefficient vector C.

%   Copyright since 2015 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2015-02-08 Original Version

% Check if the data is inside the normalized aperture......................
InAp = ~isnan(R) & R<=1;

% Get the order corresonding to the function number........................
% Check the input arguments.
if(nargin < 5)
    order = 'Carl Zeiss';
else
    if(isempty(order))
        order = 'Carl Zeiss';
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

% Pass the inputs to the function..........................................
zpf = zernike_rt_nm(R(InAp),Theta(InAp),n,m);

% Weight each aberration
Z = NaN(size(R));
for k = 1:length(C)
    zpf(:,k) = zpf(:,k).*C(k);
end
z = sum(zpf,2);
Z(InAp) = z(:);
end
