function recorder = epsilonExpandRecord(sys_c,opts)
%   epsilonExpandRecord:
%       This is an "object" for recordings stats and the iterate history of
%       epsilonExpand.  
%
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, already initialized at some perturbation
%           epsilon*U*V'.  Note that this is only used to get the data for
%           the initial point and the total costs.  Iterate info must be
%           recorded by passing in the explicit details, as sys_c is never
%           updated internally with the info for rejected steps.
%
%       opts            [required: struct of parameters]
%           A required struct settable parameters necessary to run this
%           routine. 
%
%           .maxit                      [positive integer]
%               The maximum number of iterations allowed for the
%               epsilonExpand procedure.  This is used so that a cell
%               array can usually be preallocated of sufficient size, so
%               that it doesn't have to grow when recording.  
%
%           .line_search_opts.maxit     [positive integer]
%               The maximum number of evaluations possible in each line
%               search, including t = 0.  Also used to determine the size
%               of celll array to preallocate.
%
%           .record_level               [value in {0,1,2}]
%               Determines how much metadata is recorded:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the accepted iterates 
%               2 - additionally includes the rejected line search points.
%
%           .record_UV                  [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be recorded.
%
%   OUTPUT:
%       recorder
%           The recording "object", a struct with the following methods:
%       
%       .addEpsilon(accepted,f,z,epsilon,t,df_eps_t0,newton,damped)    
%           Add the info about the full step (t = 1) to the record.  
%           [This method is only present if opts.record_level > 0]
%
%       .addEpsilonLS(accepted,f,z,epsilon,t,bisected)
%           Add the info about this line search evaluation to the record.  
%           [This method is only present if opts.record_level > 0]
%
%       rec = recorder.getRecord()
%           A cell array containing the recorded info, where the odd
%           entries are the labels and each even entry is the data
%           associated with the preceding label in the cell array.
%
%   See also epsilonExpand and epsilonExpandOptions.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   epsilonExpandRecord.m introduced in ROSTAPACK Version 1.0  
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

    cost_counter            = totalsCounter(@sys_c.getTotals);
    
    record_level            = opts.record_level;
    switch record_level
        case 0
            recorder        = struct('getRecord',@getRecord);
            return
        case 1
            points_per_iter = 1;
            getStructFn     = @infoAccept;
        case 2
            ls_opts         = opts.line_search_opts;
            points_per_iter = ls_opts.maxit;
            getStructFn     = @infoAny;
        otherwise
            error('Unrecognized value of opts.record_level!');
    end
    record_UV               = opts.record_UV;
    if record_UV 
        [U,V]               = sys_c.getUV();
    end
    
    % preallocate data structure for storing all recorded iterates
    rec                     = cell(1,points_per_iter*opts.maxit + 1);
    count                   = 0;
    ls_info_fn              = [];
  
    % Add the starting point  
    [f,z]                   = sys_c.getf();
    epsilon                 = sys_c.getEpsilon();
    addPoint('epsilon_initial',f,z,epsilon,[],true);
  
    recorder    = struct(   'addEpsilon',           @addEpsilon,        ...
                            'addEpsilonLS',         @addEpsilonLS,      ...
                            'getRecord',            @getRecord          );
 
    function addEpsilon(accepted,f,z,epsilon,t,df_eps_t0,newton,damped)
        
        % f == nan means df_eps_t0 was below tolerance
        
        % need to set this whether or not point is accepted
        ls_info_fn = @(t,bisected) lsInfoAccept(t,bisected,df_eps_t0,   ...
                                                newton,damped           );
        
        if record_level == 1 && ~accepted
            return
        end
        
        if accepted
            info        = ls_info_fn(t,false);
            ls_info_fn  = [];
        else
            info        = lsInfo(t,false);
        end
        addPoint('epsilon',f,z,epsilon,info,accepted);
    end

    function addEpsilonLS(accepted,f,z,epsilon,t,bisected)
                 
        if isempty(ls_info_fn) || (record_level == 1 && ~accepted)
            return
        end
           
        if accepted 
            info = ls_info_fn(t,bisected);
        else
            info = lsInfo(t,bisected);
        end
        addPoint('epsilon_ls',f,z,epsilon,info,accepted);
    end

    function addPoint(type,f,z,epsilon,info,accepted)
        cost        = cost_counter.getLatest();
        s           = getStructFn(type,f,z,epsilon,info,cost,accepted);
        count       = count + 1;
        rec{count}  = s;
    end

    function r = getRecord()
        r = { 'cost', cost_counter.getCumulative() };
        if record_level < 1
            return
        end 
        if record_UV 
            r = [ {'U',U,'V',V} r ];
        end
        r = [ r {'iterates', cell2mat(rec(1:count))} ];
    end
end

function s = infoAccept(type,f,z,epsilon,info,cost,~)
    s = struct( 'type',     type,       ...  
                'f',        f,          ...
                'z',        z,          ...
                'epsilon',  epsilon,    ...
                'info',     info,       ...
                'cost',     cost        );
end

function s = infoAny(type,f,z,epsilon,info,cost,accepted)
    s = struct( 'type',     type,       ...  
                'accepted', accepted,   ...
                'f',        f,          ...
                'z',        z,          ...
                'epsilon',  epsilon,    ...
                'info',     info,       ...
                'cost',     cost        );
end

function s = lsInfoAccept(t,bisection,df_epsilon_t0,newton,damped)
    s = struct( 't',                t,              ...
                'bisection',        bisection,      ...
                'df_epsilon_t0',    df_epsilon_t0,  ...
                'newton_type',      newton,         ...
                'damped',           damped          );
end

function s = lsInfo(t,bisection)
    s = struct( 't',            t,          ...
                'bisection',    bisection   );
end
