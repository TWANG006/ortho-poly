function [z3, zx3, zy3] = zernike_xy_nm(X, Y, n, m, NormFlag)
%ZERNIKEXYNM Zernike polynominals in normalized X and Y coordinates
%   with order vector (n, m).

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-01-01 Original Version

% Check normalization......................................................
if nargin==5 && ischar(NormFlag)
    if ~strcmpi(NormFlag,'normPeak') ...
            && ~strcmpi(NormFlag,'Common') ...
            && ~strcmpi(NormFlag,'normIntegral')
        error('ZernikeXYnm:NormFlag','Unrecognized normalization flag.')
    end
else
    % By default, use 'normPeak' ...
    NormFlag = 'normPeak';
end

% Compute the Zernike functions............................................
Cj = n+1j*m;
sz = size(X);
ln = length(n);
[Theta, R] = cart2pol(X,Y);
z = zernike_rt_nm(R,Theta,n,m,NormFlag);
z3 = reshape(z,sz(1),sz(2),ln);

% Compute the Zernike Cartesian derivatives................................
% [1] J. Ruoff, and M. Totzeck, "Orientation Zernike polynomials: a useful
% way to describe the polarization effects of optical imaging systems,"
% JOSA (2009).

% Calculate the parameters.
am = (m>=0) - (m<0);

Nnm = sqrt((2-(m==0)).*(n+1));
Nnm1 = sqrt((2-((m-1)==0)).*((n-1)+1));
Nnm2 = sqrt((2-((m+1)==0)).*((n-1)+1));
Nnm0 = sqrt((2-(m==0)).*((n-2)+1));

b1 = Nnm./Nnm1;
b2 = Nnm./Nnm2;
b0 = Nnm./Nnm0;
b1(isinf(b1))=0;
b2(isinf(b2))=0;
b0(isinf(b0))=0;
b0(imag(b0)~=0)=0;

m1 = am.*abs(m-1);
m2 = am.*abs(m+1);
m3 = -am.*abs(m-1);
m4 = -am.*abs(m+1);

% Calculate the derivatives.
sz = size(z3);
zx3 = zeros(sz);
zy3 = zeros(sz);
Zx3 = zeros(sz);
Zy3 = zeros(sz);
for k = 1:ln
    % Deal with the initial conditions.....................................
    j1 = find( Cj == (n(k)-1)+m1(k)*1j );
    j2 = find( Cj == (n(k)-1)+m2(k)*1j );
    j3 = find( Cj == (n(k)-1)+m3(k)*1j );
    j4 = find( Cj == (n(k)-1)+m4(k)*1j );
    j0 = find( Cj == (n(k)-2)+ m(k)*1j  );
    % Check j1.
    if isempty(j1)
        Z1 = zeros(sz(1),sz(2));
    else
        if strcmpi(NormFlag,'normPeak')
            Z1 = z3(:,:,j1).*Nnm(j1);
        elseif strcmpi(NormFlag,'Common')
            Z1 = z3(:,:,j1);
        elseif strcmpi(NormFlag,'normIntegral')
            Z1 = z3(:,:,j1)./sqrt(1/pi);
        end
    end
    % Check j2.
    if isempty(j2)
        Z2 = zeros(sz(1),sz(2));
    else
        if strcmpi(NormFlag,'normPeak')
            Z2 = z3(:,:,j2).*Nnm(j2);
        elseif strcmpi(NormFlag,'Common')
            Z2 = z3(:,:,j2);
        elseif strcmpi(NormFlag,'normIntegral')
            Z2 = z3(:,:,j2)./sqrt(1/pi);
        end
    end
    % Check j3.
    if isempty(j3)
        Z3 = zeros(sz(1),sz(2));
    else
        if strcmpi(NormFlag,'normPeak')
            Z3 = z3(:,:,j3).*Nnm(j3);
        elseif strcmpi(NormFlag,'Common')
            Z3 = z3(:,:,j3);
        elseif strcmpi(NormFlag,'normIntegral')
            Z3 = z3(:,:,j3)./sqrt(1/pi);
        end
    end
    % Check j4.
    if isempty(j4)
        Z4 = zeros(sz(1),sz(2));
    else
        if strcmpi(NormFlag,'normPeak')
            Z4 = z3(:,:,j4).*Nnm(j4);
        elseif strcmpi(NormFlag,'Common')
            Z4 = z3(:,:,j4);
        elseif strcmpi(NormFlag,'normIntegral')
            Z4 = z3(:,:,j4)./sqrt(1/pi);
        end
    end
    % Check j0.
    if isempty(j0)
        Zx = zeros(sz(1),sz(2));
        Zy = zeros(sz(1),sz(2));
    else
        Zx = Zx3(:,:,j0);
        Zy = Zy3(:,:,j0);
    end
    % Calculate the derivatives............................................
    Zx3(:,:,k) = ...
        n(k).*( b1(k).*Z1 + am(k).*sign(m(k)+1).*b2(k).*Z2 ) ...
        + b0(k).*Zx;
    
    Zy3(:,:,k) = ...
        n(k).*(-am(k).*sign(m(k)-1).*b1(k).*Z3 + b2(k).*Z4 ) ...
        + b0(k).*Zy;
    
    % Give output according to the NormFlag.
    if strcmpi(NormFlag,'normPeak')
        zx3(:,:,k) = Zx3(:,:,k)./Nnm(k);
        zy3(:,:,k) = Zy3(:,:,k)./Nnm(k);
    elseif strcmpi(NormFlag,'Common')
        zx3(:,:,k) = Zx3(:,:,k);
        zy3(:,:,k) = Zy3(:,:,k);
    elseif strcmpi(NormFlag,'normIntegral')
        zx3(:,:,k) = Zx3(:,:,k).*sqrt(1/pi);
        zy3(:,:,k) = Zy3(:,:,k).*sqrt(1/pi);
    end
    
end


end
