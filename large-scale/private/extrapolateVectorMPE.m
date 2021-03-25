 function s = extrapolateVectorMPE(V)
%   extrapolateVectorMPE(V):
%       Computes vector extrapolation s of vector sequence v_k, where v_k
%       is a column of matrix V.  
%
%       This code implements the minimal polynomal extrapolation method for 
%       vector sequences of:
%
%       S. Cabay, L.W. Jackson, A polynomial extrapolation method for
%       finding limits and antilimits of vector sequences, SIAM J. Numer.
%       Anal., 13 (1976), pp. 734-752.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   extrapolateVectorMPE.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  extrapolateVectorMPE.m                                               |
% |  Copyright (C) 2016-2019 Tim Mitchell                                 |
% |                                                                       |
% |  This file is originally from URTM.                                   |
% |                                                                       |
% |  URTM is free software: you can redistribute it and/or modify         |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  URTM is distributed in the hope that it will be useful,              |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
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

    % calculate v1-v2, v2-v3, and so forth
    U   = V(:,1:end-1) - V(:,2:end);

    % take off last column uk
    uk  = U(:,end);
    U   = U(:,1:end-1);

    % Compute linear least squares of coefficients c, k-1 of them
    % no need for warnings here so temporarily disable and reenable back to
    % the same warning state
    state = warning();
    warning('off');
    c   = U\(-uk);
    warning(state);

    % set ck to 1
    c   = [c;1];

    g   = c/sum(c);

    % compute vector extrapolate
    s = V(:,1:end-1)*g;

end
