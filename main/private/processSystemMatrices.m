function [A,B,C,D,E,D_zero] = processSystemMatrices(A,B,C,D,E)
%   processSystemMatrices(A,B,C,D,E):
%       processSystemMatrices ensures that the given matrices A,B,C,D,E are
%       finite valued and compatibly dimensioned for a linear dynamical
%       system with input and output.  If not, errors are thrown.  When
%       the provided B,C,E are identity matrices, they will be set to [] in
%       the output as the explicit forms are not needed.  When D is zero,
%       D_zero is true and D will either be set to [] (when n=m=p) or
%       sparse(p,m). The values of m and p are inferred from the input
%       matrices.
% 
%   See also systemDimensions.
%
%   
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   processSystemMatrices.m introduced in ROSTAPACK Version 2.0
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

    [n,m,p] = systemDimensions(A,B,C,D);
    
    assert( isSize(A,n,n) && ~isempty(A) && isFiniteValued(A),  ...
            'A must be a square matrix with no infs/nans'       );
    
    if ~isempty(B)
        assert(isSize(B,n,m) && isFiniteValued(B),              ...
            'B must be n by m with no infs/nans'                );
        if isEye(B)
            B = [];
        end
    end
    
    if ~isempty(C)
        assert(isSize(C,p,n) && isFiniteValued(C),              ...
            'C must be p by n with no infs/nans'                );
        if isEye(C)
            C = [];
        end
    end
    
    if ~isempty(D)
        assert(isSize(D,p,m) && isFiniteValued(D),              ...
            'D must be p by m with no infs/nans'                );
        D_zero = isZero(D);
        if D_zero
            if p == n && m == n
                D = [];
            elseif ~issparse(D)
                D = sparse(p,m);
            end
        end
    else
        D_zero = true;
    end
        
    if ~isempty(E)
        assert(isSize(E,n,n) && isFiniteValued(E),              ...
            'E must be n by n with no infs/nans'                );
        if isEye(E)
            E = [];
        end
    end
    
end