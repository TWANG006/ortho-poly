function z = zernike_rt_nm(R,Theta,n,m,NormFlag)
%ZERNIKERTnm Zernike polynominals in normalized R and T coordinates 
%   with order vector (n, m).

%   Copyright since 2015 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2015-02-08 Original Version

% This function is changed from Paul Fricker(2/28/2012)'s version.

% Check and prepare the inputs.............................................
if ( ~any(size(n)==1) ) || ( ~any(size(m)==1) )
    error('zernfun:NMvectors','N and M must be vectors.')
end

if length(n)~=length(m)
    error('zernfun:NMlength','N and M must be the same length.')
end

n = n(:);
m = m(:);
if any(mod(n-m,2))
    error('zernfun:NMmultiplesof2', ...
          'All N and M must differ by multiples of 2 (including 0).')
end

if any(m>n)
    error('zernfun:MlessthanN', ...
          'Each M must be less than or equal to its corresponding N.')
end

if any( R>1 | R<0 )
    error('zernfun:Rlessthan1','All R must be between 0 and 1.')
end


% vectorize R and Theta.
R = R(:);
Theta = Theta(:);

if ( ~any(size(R)==1) ) || ( ~any(size(Theta)==1) )
    error('zernfun:RTHvector','R and THETA must be vectors.')
end

length_r = length(R);
if length_r~=length(Theta)
    error('zernfun:RTHlength', ...
          'The number of R- and THETA-values must be equal.')
end

% Check normalization......................................................
if nargin==5 && ischar(NormFlag)
    if ~strcmpi(NormFlag,'normPeak') ...
            && ~strcmpi(NormFlag,'Common') ...
            && ~strcmpi(NormFlag,'normIntegral')
        error('ZernikeRTnm:NormFlag','Unrecognized normalization flag.')
    end
else
    % By default, use 'normPeak' ...
    NormFlag = 'normPeak';
end

% Compute the Zernike Polynomials..........................................

% Determine the required powers of r.......................................
m_abs = abs(m);
rpowers = [];
for j = 1:length(n)
    rpowers_update = [rpowers,m_abs(j):2:n(j)];
    rpowers = rpowers_update;
end
rpowers = unique(rpowers);

% Pre-compute the values of r raised to the required powers,
% and compile them in a matrix.............................................
if rpowers(1)==0
    rpowern = arrayfun(@(p)R.^p,rpowers(2:end),'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
    rpowern = [ones(length_r,1) rpowern];
else
    rpowern = arrayfun(@(p)R.^p,rpowers,'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
end

% Compute the values of the polynomials....................................
z = zeros(length_r,length(n));
for j = 1:length(n)
    s = 0:(n(j)-m_abs(j))/2;
    pows = n(j):-2:m_abs(j);
    for k = length(s):-1:1
        p = (1-2*mod(s(k),2))* ...
                   prod(2:(n(j)-s(k)))/              ...
                   prod(2:s(k))/                     ...
                   prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                   prod(2:((n(j)+m_abs(j))/2-s(k)));
        idx = (pows(k)==rpowers);
        z(:,j) = z(:,j) + p*rpowern(:,idx);
    end
    
    if strcmpi(NormFlag,'normPeak')
        % do nothing.
    elseif strcmpi(NormFlag,'Common')
        z(:,j) = z(:,j)*sqrt((1+(m(j)~=0))*(n(j)+1));
    elseif strcmpi(NormFlag,'normIntegral')
        z(:,j) = z(:,j)*sqrt((1+(m(j)~=0))*(n(j)+1)/pi);
    end
end % END: Compute the Zernike Polynomials

% Compute the Zernike functions............................................
idx_pos = m>0;
idx_neg = m<0;

if any(idx_pos)
    z(:,idx_pos) = z(:,idx_pos).*cos(Theta*m_abs(idx_pos)');
end
if any(idx_neg)
    z(:,idx_neg) = z(:,idx_neg).*sin(Theta*m_abs(idx_neg)');
end

end % EOF
