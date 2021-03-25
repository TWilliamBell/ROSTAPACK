function [z,z_dr,z_dt,z_dr2,z_dt2,z_drdt] = partialsPolar(r,t)
%   partialsPolar(r,t):
%       Given a number in polar coordinates, this returns it as a complex
%       number z, with its first and second partial derivatives with
%       respect to r and t: z_dr, z_dt, z_dr2, z_dt2, z_drdt.
%
%   See also partialsCartesian.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   partialsPolar.m introduced in ROSTAPACK Version 2.0
%
% =========================================================================
% |  partialsPolar.m                                                      |
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

    % x = radius and y = angle
    eitheta = polar2cplx(1,t);      % ensures exactly -1 if t == pi
    z       = r*eitheta;
    z_dr    = eitheta;
    z_dt    = 1i*z;
    z_dr2   = 0;
    z_dt2   = -z;
    z_drdt  = 1i*eitheta;
    
end