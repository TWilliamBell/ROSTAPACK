function obj = uvExpandRecord(sys_c,opts)
%   uvExpandRecord:
%       This is an "object" for recordings stats and the iterate history of
%       uvExpand.  
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
%           The following uvExpand parameters are required:
%
%           .maxit  
%           .interp_threshold
%           .extrap_size
%           .line_search_opts.maxit  
%           .line_search_opts.UV_dual_mode  
%                  
%           For more details on these, see uvExpandOptions.
%
%           The following parameters are additionally required: 
%
%           .record_level               [value in {0,1,2}]
%               Determines how much metadata is recorded:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the accepted iterates 
%               2 - additionally includes the rejected points.
%
%           .record_UV                  [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be recorded.
%
%   OUTPUT:
%       recorder
%           The recording "object", a struct with the following methods:
%       
%       .addUV(accepted,f,z,U,V,t,df_UV_t0,df_eps_t0)
%           Add the info about the full step (t = 1) to the record. Note
%           that we require t because it may actually be zero, which
%           happens when the full step perturbation is either identical to
%           the current perturbation or if it wasn't evaluated (if the
%           magnitude of the line search was below its tolerance). 
%           [This method is only present if opts.record_level > 0]
%
%       .addUVLS(accepted,f,z,U,V,t,bisected,flipped_UV)
%           Add the info about this line search evaluation to the record.  
%           [This method is only present if opts.record_level > 0]
%
%       .addInterp(accepted,f,z,U,V,t)
%           Add the info about the interpolation attempt.
%           [This method is only present if opts.record_level > 0]
%
%       .addExtrap(accepted,f,z,U,V,vecs_used)
%           Add the info about the extrapolation step attempt.
%           [This method is only present if opts.record_level > 0]
%
%       rec = recorder.getRecord()
%           A cell array containing the recorded info, where the odd
%           entries are the labels and each even entry is the data
%           associated with the preceding label in the cell array.
%
%   See also uvExpand and uvExpandOptions.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   uvExpandRecord.m introduced in ROSTAPACK Version 1.0
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
    record_UV               = opts.record_UV;
    switch record_level
        case 0
            obj             = struct('getRecord',@getRecord);
            return
        case 1              
            max_ls_evals    = 1;
            getStructFn     = ternOp(record_UV,@infoAcceptUV,@infoAccept);
        case 2
            max_ls_evals    = getMaxLineSearchEvals(opts.line_search_opts);
            getStructFn     = ternOp(record_UV,@infoAnyUV,@infoAny);
        otherwise
            error('Unrecognized value of opts.record_level!');
    end
    
    interp_on               = opts.interp_threshold < inf;
    extrap_on               = opts.extrap_size > 1;
    points_per_iter         = max_ls_evals + interp_on + extrap_on;
    
    % preallocate data structure for storing all recorded iterates
    rec                     = cell(1,points_per_iter*opts.maxit + 1);
    count                   = 0;
    ls_info_fn              = [];
  
    % Add the starting point  
    [f,z]                   = sys_c.getf();
    epsilon                 = sys_c.getEpsilon();
    [U,V]                   = sys_c.getUV();
    addPoint('UV_initial',f,z,U,V,[],true);
    
    obj     = struct(   'addUV',                @addUV,                 ...
                        'addUVLS',              @addUVLS,               ...
                        'addInterp',            @addInterp,             ...
                        'addExtrap',            @addExtrap,             ...
                        'getRecord',            @getRecord              );

    function addUV(accepted,f,z,U,V,t,df_UV_t0,df_eps_t0)
           
        % t == 0 means U and V weren't used because either: 
        % - U,V were identical to previous U,V (f will equal previous f)
        % - df_UV_t0 was below tolerance (f will be nan)
        
        % need to set this whether or not point is accepted
        ls_info_fn      = @(t,tf) lsInfoAccepted(t,tf,df_UV_t0,df_eps_t0);
        
        if record_level == 1 && ~accepted
            return
        end
       
        if accepted 
            info        = ls_info_fn(t,false);
            ls_info_fn  = [];
        else
            info        = lsInfo(t,false);
        end
        addPoint('UV',f,z,U,V,info,accepted);
    end

    function addUVLS(accepted,f,z,U,V,t,bisected,flipped_UV)
        
        if isempty(ls_info_fn) || (record_level == 1 && ~accepted)
            return
        end
            
        if flipped_UV
            % Negative t means the signs of U and V were flipped
            t = -t;
        end
        
        if accepted 
            info = ls_info_fn(t,bisected);
        else
            info = lsInfo(t,bisected);
        end
        addPoint('UV_ls',f,z,U,V,info,accepted); 
    end

    function addInterp(accepted,f,z,U,V,t)
        if isempty(ls_info_fn) || (record_level == 1 && ~accepted)
            return
        end
        
        info.t = t;
        addPoint('interp',f,z,U,V,info,accepted);
    end

    function addExtrap(accepted,f,z,U,V,vecs_used)
        if isempty(ls_info_fn) || (record_level == 1 && ~accepted)
            return
        end
        
        info.vecs_used = vecs_used;
        addPoint('extrap',f,z,U,V,info,accepted);
    end

    function addPoint(type,f,z,U,V,info,accepted)
        cost        = cost_counter.getLatest();
        s           = getStructFn(type,f,z,U,V,info,cost,accepted);
        count       = count + 1;
        rec{count}  = s;
    end

    function r = getRecord()
        r = { 'cost', cost_counter.getCumulative() };
        if record_level < 1
           return
        end
        iterates = cell2mat(rec(1:count));
        r = [ {'epsilon',epsilon} r {'iterates',iterates} ];       
    end
end

function s = infoAccept(type,f,z,~,~,info,cost,~)
    s = struct( 'type',     type,       ...  
                'f',        f,          ...
                'z',        z,          ...
                'info',     info,       ...
                'cost',     cost        );
end

function s = infoAcceptUV(type,f,z,U,V,info,cost,~)
    s = struct( 'type',     type,       ...  
                'f',        f,          ...
                'z',        z,          ...
                'U',        U,          ...
                'V',        V,          ...
                'info',     info,       ...
                'cost',     cost        );
end

function s = infoAny(type,f,z,~,~,info,cost,accepted)
    s = struct( 'type',     type,       ...  
                'accepted', accepted,   ...
                'f',        f,          ...
                'z',        z,          ...
                'info',     info,       ...
                'cost',     cost        );
end

function s = infoAnyUV(type,f,z,U,V,info,cost,accepted)
    s = struct( 'type',     type,       ...  
                'accepted', accepted,   ...
                'f',        f,          ...
                'z',        z,          ...
                'U',        U,          ...
                'V',        V,          ...
                'info',     info,       ...
                'cost',     cost        );
end

function s = lsInfo(t,bisection)
    s = struct( 't',            t,          ...
                'bisection',    bisection   );
end

function s = lsInfoAccepted(t,bisection,df_UV_t0,df_epsilon_t0)
    s = struct( 't',                t,              ...
                'bisection',        bisection,      ...
                'df_UV_t0',         df_UV_t0,       ...
                'df_epsilon_t0',    df_epsilon_t0   );
end

function max_evals = getMaxLineSearchEvals(ls_opts)
    max_evals       = ls_opts.maxit;
    dual_mode       = ls_opts.UV_dual_mode;
    if isnumeric(dual_mode) || dual_mode
        max_evals   = 2*max_evals;
    end
end
