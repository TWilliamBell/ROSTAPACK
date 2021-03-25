function partialFn = secondPartialSingularValue(A,V,d)
%   secondPartialSingularValue:
%       Let A be an r by c matrix where V and d are respectively the
%       eigenvectors and eigenvalues of [0 A; A' 0].  This returns a
%       function handle for computing the second partial derivative of any
%       singular value of A with respect to two choices of variables x and
%       y (which can be the same).
%
%   INPUT:
%       A                           [matrix]
%           Matrix for which a second partial derivative of one its
%           singular values is desired
%
%       V                           [matrix]
%           Matrix of eigenvectors for [0 A; A' 0]
% 
%       d                           [vector]
%           Corresponding eigenvalues for [0 A; A' 0], ordered from largest
%           to smallest
%   
%   OUTPUT:
%       partialFn                   [function handle]
%           where 
%               dxdy = partialFn(k,Adx,Adx2)    
%               - dxdy second partial of kth singular value w.r.t. x
%               - Adx is the matrix first derivative with respect to x
%               - Adx2 is the matrix second derivative with respect to x
%               dxdy = partialFn(k,Adx,Ady,Adxdy) 
%               - second partial of kth singular value w.r.t. x and y
%               - Adx is the matrix first derivative w.r.t. x
%               - Ady is the matrix first derivative w.r.t. y
%               - Adxdy is the matrix second derivative w.r.t. x and y
% 
%   See also svdToEigs and secondPartialEigenvalue.
%    
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   secondPartialSingularValue.m introduced in ROSTAPACK Version 2.0
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

    [r,c]       = size(A);
    sv_vals     = min(r,c);    
    partialFn   = @partial; 
    
    function dxdy = partial(k,Adx,Ad1,Ad2)
        
        if k > sv_vals
            error('There is no %dth singular value!',k);
        end
        
        if nargin < 4 
            Ady     = [];
            Adxdy   = Ad1;
        else
            Ady     = Ad1; 
            Adxdy   = Ad2;
        end
        dxdy        = secondPartialEigenvalue(V,d,k,Adx,Ady,Adxdy);
    end
end