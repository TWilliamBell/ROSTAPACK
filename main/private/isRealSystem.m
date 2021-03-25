function [is_real,A,B,C,D,E] = isRealSystem(A,B,C,D,E)
%   isRealSystem:
%       Returns true if all matrices A,B,C,D,E are real valued or [].  Any
%       matrix which has a zero but nevertheless allocated imaginary part,
%       is returned with this removed. 
% 
%   INPUT:
%       A,B,C,D,E                   [matrices]
%   
%   OUTPUT:
%       is_real                     [logical]
%           True if all matrices A,B,C,D,E have zero imaginary parts or are
%           [] matrices
% 
%       A,B,C,D,E                   [matrices]
%           The given A,B,C,D,E matrices but with zero imaginary parts 
%           removed if they were allocated.
%           
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   isRealSystem.m introduced in ROSTAPACK Version 2.0
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

    is_real = true;
    A       = checkMatrix(A);
    B       = checkMatrix(B);
    C       = checkMatrix(C);
    D       = checkMatrix(D);
    E       = checkMatrix(E);
    
    function M = checkMatrix(M)
        if isRealValued(M)
            M = real(M);
        else
            is_real = false;
        end
    end
end