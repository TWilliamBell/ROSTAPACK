function [n,m,p] = systemDimensions(A,B,C,D)
%   systemDimensions(A,B,C,D):
%       Given input-output system A,B,C,D, returns:
%
%           n: the state-space dimension
%           m: the number of inputs
%           p: the number of outputs
% 
%       A must be provided as an explicit matrix.  B,C,D can optionally be
%       given as [], which are used as shortcuts to denote:
% 
%           appropriately dimensioned (sp)eye matrices: B and/or C
%           an appropriately dimensioned zero matrix:   D 
%
%       - If B and D are both empty, m is set to n
%       - If C and D are both empty, p is set to n
%       - If B and D are inconsistently dimensioned, m is set to the larger
%         of the two incompatible dimensions
%       - If C and D are inconsistently dimensioned, p is set to the larger 
%         of the two incompatible dimensions
%
%   See also processSystemMatrices. 
%
% 
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   systemDimensions.m introduced in ROSTAPACK Version 2.0
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
    
    n       = length(A);
    [p,m]   = size(D);
    m       = max(m,size(B,2));
    p       = max(p,size(C,1));
    if m == 0
        m   = n;
    end
    if p == 0
        p   = n;
    end
         
end