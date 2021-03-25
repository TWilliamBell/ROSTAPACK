function lines = copyrightNotice(filename)
%   copyrightNotice:
%       Returns a cell array of strings, where each entry is a line of text
%       to print as the copyright notice for this software package and a
%       particular file, given by the caller-supplied string filename.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   copyrightNotice.m introduced in ROSTAPACK Version 1.0
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

lines = {
    'ROSTAPACK: RObust STAbility PACKage --- Version 1.0',              ...
    'Licensed under the AGPLv3, Copyright (C) 2014-2018 Tim Mitchell',  ...
    sprintf('%s()',filename),                                           };
end
