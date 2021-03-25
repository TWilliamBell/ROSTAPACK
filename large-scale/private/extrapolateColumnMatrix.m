function [u,v] = extrapolateColumnMatrix(U,v_scalars) 
%   extrapolateColumnMatrix:
%       This is a subroutine of extrapolateUpToRank2Matrix.
% 
%       This subroutine performs vector extrapolation on a sequence of
%       columns vectors u_k*v_k, where u_k is the kth column of matrix U
%       and scalar v_k is the kth entry in v_scalars.  U and v_scalars may
%       contain complex entries.
%
%       This returns the extrapolated vector u, normalized to have unit
%       norm, and v = 1;
%
%       Besides the normalized output, this method is basically a straight
%       forward version of MPE vector extrapolation.  
%
%   See also extrapolateVectorMPE.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   extrapolateColumnMatrix.m introduced in ROSTAPACK Version 1.0 
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
    
    % Since m == 1, we can form u*v' matrix explicitly and extrapolate it
    
    UVh_samples = U .* repmat(conj(v_scalars),size(U,1),1);
    
    % we can assume that v_scalars are all identically 1, meaning that
    % extrapolating u*v' is equivalent to u directly
    u           = extrapolateVectorMPE(UVh_samples);  
    u           = u / norm(u);  % ensure extrapolated u has unit norm
    v           = 1;
end
