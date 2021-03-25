function intervals = splitAnInterval(intervals,w,width)
%   splitAnInterval:
%       Given a cell array of intervals, i.e. 2 element vectors, this
%       function splits any interval into two seperate ones at w, whenever
%       w is considered sufficiently within the middle portion of said
%       interval.
%
%   INPUT:
%       intervals                   [cell array of 2 element vectors]
%           A list of real-valued intervals
% 
%       w                           [real scalar]
%           The point at which to split any interval
% 
%       width                       [real value in [0,1]]
%           Width of middle portion of interval w must be inside in order
%           to split an interval into two at w.  When width is 0, no
%           interval is ever split; when width is 1, an interval is split
%           at w whenever w is anywhere inside the entire open version of
%           the interval (to ensure a zero length interval is never
%           created).  When width is 0.5, w must be in the middle 50% of
%           an interval in order for it be split.
%   
%   OUTPUT:
%        intervals                   [cell array of 2 element vectors]
%           The list of the real-valued intervals, with possibly one or
%           more intervals split.
%
%   See all formAllIntervals.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   splitAnInterval.m introduced in ROSTAPACK Version 2.0
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

    if width <= 0
        return
    end
    
    tf  = cellfun(@shouldSplit,intervals);
    k   = find(tf); 
    
    if ~isempty(k)
        before      = intervals(1:k-1);
        intk        = intervals{k};
        after       = intervals(k+1:end);   
        s           = {[intk(1) w]; [w intk(2)]};
        intervals   = [ before; s; after ];
    end

    function tf = shouldSplit(interval)
        low     = interval(1);
        high    = interval(2);
        d       = 0.5*(1 - width)*(high - low);
        tf      = (low + d) < w && w < (high - d);
    end
end

