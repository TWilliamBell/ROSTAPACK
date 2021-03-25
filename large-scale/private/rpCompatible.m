function [x,y,yhx] = rpCompatible(x,y,varargin)
%   rpCompatible:
%       Scale y so that the pair of vectors x,y is RP-compatible, that is
%       both x and y have unit norm and yhx = y'*x is a real, positive.  
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   rpCompatible.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  rpCompatible.m                                                       |
% |  Copyright (C) 2016 Nicola Guglielmi, Michael Overton, Tim Mitchell   |
% |                                                                       |
% |  This routine (this single file) is taken from the PSAPSR software    |
% |  package, which is licensed under the GPL v3.  As such, the contents  |
% |  of the file (code and comments) are also licensed under the GPL v3.  |
% |  However, note that this is an exceptional case; PSARNOT and most of  |
% |  its subroutines are licensed under the AGPL v3.                      |
% |                                                                       |
% |  This routine is free software: you can redistribute it and/or modify |
% |  it under the terms of the GNU General Public License as published by |
% |  the Free Software Foundation, either version 3 of the License, or    |
% |  (at your option) any later version.                                  |
% |                                                                       |
% |  This routine is distributed in the hope that it will be useful,      |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    |
% |  General Public License for more details.                             |
% |                                                                       |
% |  You should have received a copy of the GNU General Public License    |
% |  along with this program.  If not, see                                |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
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

    x   = x/norm(x);
    y   = y/norm(y);
    yhx = y'*x; 

    if isreal(yhx)
        % avoid imaginary rounding error
        if yhx < 0
            y   = -y;
            yhx = -yhx;
        end
    else
        scalar = exp(1i*angle(yhx));
        y   = scalar*y; 
        yhx = real(conj(scalar)*yhx); % ensure yhx is strictly real
    end
end
