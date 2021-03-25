function z = polar2cplx(r,theta)
%   polar2cplx(r,theta):
%       z = polar2cplx(r,theta)
%       Form a complex number z from its polar coordinates: radius r and
%       angle theta.  If theta is equal to 0 or pi, then z will be as
%       plus/minus r, with exactly zero imaginary part, whereas just using
%       z = r*exp(1i*theta) for theta == pi would result in z having a
%       slightly nonzero imaginary part.
%
%       NOTE: Unlike MATLAB's builtin pol2cart routine, here the radius is 
%       the first input and the angle is the second.
%      
%   See also cplx2polar, cart2cplx, cplx2cart, cart2polar, and polar2cart.
%
% 
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   polar2cplx.m introduced in ROSTAPACK Version 2.0
%
% =========================================================================
% |  polar2cplx.m                                                         |
% |  Copyright (C) 2019 Tim Mitchell                                      |
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

    if theta == 0
        z = r;
    elseif theta == pi
        z = -r;
    else
        assert( isscalar(r) || isscalar(theta) ||       ...
                max(size(r) - size(theta)) == 0,        ...
                'Input arrays must have the same size.' );
        z = r.*exp(1i*theta);
    end
   
end
