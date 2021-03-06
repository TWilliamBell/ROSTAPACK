function obj = totalsCounter(get_totals_fn)
%   totalsCounter:
%       This is a helper object for managing totals, counts, and other
%       increasing numerical metadata generated by some process X. 
%
%   INPUT:
%       get_totals_fn 
%           A function with no outputs that, every time it is called,
%           returns a struct with fields that are all numeric types.  The
%           struct must always has the same fields present.  Typically,
%           this function returns a struct containing the total costs
%           accrued so far by some process X. 
%   
%   OUTPUT:
%       obj
%           An "object", a struct with the following methods attached:
%
%           totals = obj.getCumulative()
%               This returns a struct of the current total accrued costs
%               that process X has incurred since the creation of this
%               totalsCounter object.
%
%           totals = obj.getLatest()
%               This returns a struct of the total accrued costs that
%               process X has incurred since obj.getLatest() was last
%               called.  If obj.getLatest() has not been called previously,
%               this will return the total accrued costs incurred since
%               this totalsCounter object was created.
%           
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   totalCounters.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  totalsCounter.m                                                      |
% |  Copyright (C) 2016 Tim Mitchell                                      |
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

    totals_initial  = get_totals_fn();
    totals_latest   = totals_initial;
    fields          = fieldnames(totals_initial);
    
    obj = struct(   'getCumulative',    @getCumulative, ...
                    'getLatest',        @getLatest      );
                
    function data = getLatest()  
        [data, totals_current]  = getTotalsSince(totals_latest);
        totals_latest           = totals_current;
    end

    function data = getCumulative()       
        data = getTotalsSince(totals_initial);
    end

    % private helper function
    
    function [data, totals_current] = getTotalsSince(totals_start)    
        totals_current  = get_totals_fn();
        data            = totals_current;
        
        for j = 1:length(fields)
            field = fields{j};
            data.(field) = totals_current.(field) - totals_start.(field);
        end
    end
end
