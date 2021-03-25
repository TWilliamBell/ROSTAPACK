function [out_fn, get_count_fn] = functionCounter(fn)
%   functionCounter(fn):
%       Allows one to easily count the number of times a function has been
%       called by replacing the function handle fn with out_fn, a wrapped
%       version that has an internal counter.
%       
%   USAGE:
%       [out_fn, getCounts] = functionCounter(your_fn_handle);
%       out_fn();
%       put_fn();
%       counts = getCounts();  % returns 2 since out_fn() called twice
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   functionCounter.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  functionCounter.m                                                    |
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

    count = 0;
    
    out_fn = @fnWithCounter;
    get_count_fn = @getCount;
    
    function varargout = fnWithCounter(varargin)
        [varargout{1:nargout}] = fn(varargin{:}); 
        count = count + 1;
    end

    function c = getCount()
        c = count;
    end
end
