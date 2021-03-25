function [u,v] = normalizeRank1Extrapolation(row,col,row_index)
%   normalizeRank1Extrapolation:
%       This is a subroutine of extrapolateRank1Matrix and
%       extrapolateRank2Matrix.
%       
%       Given a row and column vector, and an index in that row, this
%       routine produces vectors u and v, both with unit norm, such that a
%       row and column of u*v' will hopefully approximate the input row and
%       col.  The purpose of this routine is to take the two vector
%       extrapolations that have been done for a chosen row and column of
%       rank-one outer product sequence of matrices, and produce a
%       normalized outer product u*v' version of the implicit rank-one
%       extrapolation.
%     
%   See also extrapolateRank1Matrix and extrapolateRank2Matrix.  
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   normalizeRank1Extrapolation.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  ROSTAPACK: RObust STAbility PACKage                                  |
% |  Copyright (C) 2014-2018 Tim Mitchell                                 |
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

    % The following reconstruction is a good guess but the system 
    % we are trying to solve is actually overdetermined.
    
    % If col and row are purely real, we don't want to unnecessarily
    % introduce imaginary numbers (since this will double the number of 
    % variables for the subsequent optimization step)
    if isreal(col(row_index)) && col(row_index) < 0 
        vj  = sqrt(-col(row_index));
    else
        vj  = sqrt(col(row_index));
    end
    u   = col / conj(vj);
    u   = u / norm(u);
    
    v   = conj( row / u(row_index) );
    v   = v / norm(v);
    
    % u and v will only be complex if row and/or col were.
    
    % We could optionally solve for u and v by specifying it as a
    % constrained optimization problem but we have not observed that doing
    % so has any fungible benefit and it is dramatically more expensive.
    
end
