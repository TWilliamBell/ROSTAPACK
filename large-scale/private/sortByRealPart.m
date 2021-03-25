function [s,ind] = sortByRealPart(lambda)
%   sortByRealPart:
%       Sort complex numbers by real part and break ties by angles. It
%       firsts sort by angle so that ties have consistent ordering by
%       angle, which puts conjugates with positive real part first.
%   
%       Compared to the PSAPSR v1.2 code, this routine has been updated to 
%       improve code conciseness and efficiency.
%
%
%   Start of ROSTAPACK comments:
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   sortByRealPart.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  sortByRealPart.m                                                     |
% |  Copyright (C) 2016 Nicola Guglielmi, Michael Overton, Tim Mitchell   |
% |                                                                       |
% |  This routine (this single file) is taken from the version 1.4 of the |
% |  PSAPSR software package, which is licensed under the GPL v3.  As     |
% |  such, the contents of this individual file (code and comments) are   |
% |  also licensed under the GPL v3.  Note however that this is an        |
% |  exceptional case; ROSTAPACK and most of its subroutines are licensed |
% |  under the AGPL v3.                                                   |
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

    % first sort by angle so ties have consistent ordering by angle, as
    % will exact conjugates
    [~,angle_indices]   = sort(angle(lambda),'descend');
    lambda              = lambda(angle_indices);
    
    % then put desired values first (sort() is stable)
    [~,indices]         = sort(real(lambda),'descend');  % largest modulus
       
    ind = angle_indices(indices);
    s   = lambda(indices);
end
