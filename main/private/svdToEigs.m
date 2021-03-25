function [X,d] = svdToEigs(U,S,V)
%   svdToEigs(U,S,V):
%       Given a singular value decomposition USV' of some matrix A, form
%       the eigenvectors and eigenvalues of [0 A; A' 0].  
%
%   INPUT:
%       U,S,V                       [matrices]
%           Singular value decomposition USV' of some matrix A
%   
%   OUTPUT:
%       X                           [matrix]
%           Matrix of eigenvectors for [0 A; A' 0]
% 
%       d                           [vector]
%           Corresponding eigenvalues for [0 A; A' 0], ordered from largest
%           to smallest
% 
%   See also secondPartialSingularValue and secondPartialEigenvalue.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   svdToEigs.m introduced in ROSTAPACK Version 2.0
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

    r   = size(U,1);
    c   = size(V,1);
    k   = min(r,c);
    nz  = max(r,c) - k;
    
    if k == 1
        s = S(1);
    else
        s = diag(S);
    end
    
    d   = [s; zeros(nz,1); -flip(s)];
     
    if r == c
        N       = [];
    elseif r < c
        [V,Vn]  = splitAtK(V,k);
        N       = [zeros(r,size(Vn,2)); Vn];
    else
        [U,Un]  = splitAtK(U,k);
        N       = [Un; zeros(c,size(Un,2))];
    end
    % Ensure X will have unit norm columns (faster to rescale now)
    % The columns of N are already unit norm
    U   = U/sqrt(2);
    V   = V/sqrt(2);
    X   = [ [ U; V ] N [fliplr(U); -fliplr(V)] ];  % LR flips necessary!!!
    
end

function [V1,V2] = splitAtK(V,k)
    V2  = V(:,k+1:end);
    V1  = V(:,1:k);
end