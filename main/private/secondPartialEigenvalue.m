function dxdy = secondPartialEigenvalue(V,d,k,Adx,Ady,Adxdy)
%   secondPartialEigenvalue:
%       Let A be an n by n Hermitian matrix with eigenvalues d (vector) and
%       eigenvectors V (square matrix), both presorted in some desired
%       order.  This computes the second partial derivative of the kth
%       eigenvalue with respect to two choices of variables x and y (which
%       can be the same).
%
%   INPUT:
%       V                           [matrix]
%           Matrix of eigenvectors 
% 
%       d                           [vector]
%           Corresponding eigenvalues (in some desired order)
%
%       k                           [positive integer]
%           Compute the second partial of eigenvalue d(k)
%
%       Adx                         [matrix]
%           The matrix first derivative with respect to x
%
%       Ady                         [matrix or []]
%           The matrix first derivative with respect to x
%           When this is empty, it is assumed that Adxdy is really Adx2,
%           i.e. the matrix second derivative with respect to x
% 
%       Adxdy                       [matrix]  
%           The matrix second derivative with respect to x and y
% 
%       V,Adx,Ady,Adxdy are all assumed to be n by n with d being n by 1.
%       However, if n is even, Adx,Ady,Adxdy can be given as n/2 by n/2
%       matrices, in which case the second partial derivative of eigenvalue
%       d(k) of [0 A; A' 0] is computed.
%   
%   OUTPUT:
%       dxdy                        [real scalar]
%           The second partial derivative of eigenvalue d(k)
% 
%   See also secondPartialSingularValue and svdToEigs.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   secondPartialEigenvalue.m introduced in ROSTAPACK Version 2.0
%
% =========================================================================
% |  ROSTAPACK: RObust STAbility PACKage                                  |
% |  Copyright (C) 2014-2019 Tim Mitchell                                 |
% |                                                                       |
% |  This file is part of ROSTAPACK                                       |
% |                                                                       |
% |  ROSTAPACK is free software: you can redistribute it and/or modify    |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  ROSTAPACK is distributed in the hope that it will be useful,         |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/agpl.html>.                             |
% =========================================================================
    
    [r,c]       = size(Adx);
    if r == length(d)
        multFn  = @multiplyEig; 
    else
        idx1    = 1:r;
        idx2    = r+1:r+c;
        multFn  = @multiplySvd;
    end
    
    dk          = d(k);
    vk          = V(:,k);
    vkhAdxV     = multFn(vk,Adx)*V;
    
    if isempty(Ady)
        terms   = (vkhAdxV.*conj(vkhAdxV))./(dk - d.');
    else
        terms   = (vkhAdxV.*conj(multFn(vk,Ady)*V))./(dk - d.');
    end
    terms(k)    = [];   % get rid of the kth term (which is inf/nan)
    
    dxdy        = multFn(vk,Adxdy)*vk + 2*sum(terms);
       
    % Private functions
    
    function vhA = multiplyEig(v,A)
        vhA     = v'*A;
    end
    
    function vhA = multiplySvd(v,A)
        vhA     = [ (A*v(idx2))', v(idx1)'*A ];
    end
  
end

