function recorder = recordIterates(sys_c,epsilon_fn,uv_fn,opts)
%   recordIterates:
%       This is an "object" for recordings stats and the iterate history of
%       the findUpperBound and hybridExpansionContraction routines.
%
%   INPUT:
%       sys_c           [required]
%           A system[Type] object.  Note that this is only used to get the
%           total costs incurred.
%
%       epsilon_fn      [required, function handle]
%           This is a function handle to the routine for updating epsilon.
%           Its last output argument should return the routine's metadata
%           regarding its computation.
%      
%       uv_fn           [required, function handle]
%           This is a function handle to the routine for updating the
%           unit-norm perturbation U*V'.  Its last output argument should
%           return the routine's metadata regarding its computation.
%      
%       opts            [required, struct]
%           The following fields are required:
%
%           .maxit                      [positive integer]
%               The maximum number of iterations allowed for the
%               findUpperBound or hybridExpansionContraction routines, as
%               appropriate.  
%
%           .record_level               [value in {0,1,2}]
%               Determines how much metadata is recorded:
%               0 - only basic metadata regarding the total incurred costs 
%               1 - same as level 0 (as other levels will handle adding
%                   data when .record_level is 1) 
%               2 - adds a history of the metadata returned by epsilon_fn
%                   and uv_fn subroutines whenever they are called.
%
%           Note that the level of metadata returned by the epsilon_fn and
%           uv_fn subroutines is governed by their own values of
%           opts.record_level.  We don't modify their opts structs using
%           recordIterates's opts.record_level since we may be passed
%           anonymous function handles of the epsilon_fn and uv_fn, which
%           already have their opts struct embedded.
%
%   OUTPUT:
%       recorder
%           The recording "object", a struct with the following methods:
%       
%       .epsilon_fn
%           A wrapped version of epsilon_fn.  This wrapped version must be
%           called to support recording.
%
%       .uv_fn 
%           A wrapped version of uv_fn.  This wrapped version must be
%           called to support recording.
%
%       rec = recorder.getRecord(halt_status,iters_incurred)
%           Returns a struct of metadata collected with the following
%           fields:
%               .halt_status     
%                   The halt_status of the routine, as supplied by the
%                   input argument.
%               .iters
%                   The number of iterations incurred by the routine, as
%                   supplied by the input argument.
%               .epsilon_phases
%                   How many times epsilon_fn was called.
%               .UV_phases
%                   How many times uv_fn was called.
%               .stats
%                   A struct containing subfields .epsilon and .UV, which
%                   each contain the csum of all of counted stats,
%                   respectively returned by all the calls to epsilon_fn
%                   and uv_fn.
%               .cost
%                   A struct containing subfields .epsilon and .UV, which
%                   each contain the sum of all of costs, respectively
%                   returned by all the calls to epsilon_fn and uv_fn
%               .phases             
%                   Only present if opts.record_level > 1.  This is a cell
%                   array of the metadata/iterate history structs returned
%                   by epsilon_fn and uv_fn, in the order they were called.
%
%   See also findUpperBound and hybridExpansionContraction.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   recordIterates.m introduced in ROSTAPACK Version 1.0
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

    cost_counter        = totalsCounter(@sys_c.getTotals);
    
    % The last argument contains stats and history data
    epsilon_cell_args   = cell(1,nargout(epsilon_fn));
    UV_cell_args        = cell(1,nargout(uv_fn));
    record_level        = opts.record_level;
    
    if record_level > 1
        % plus one more for the upperbound routine
        recs            = cell(1,2*opts.maxit + 1);
    end
    
    count               = 0;
    epsilon_phases      = 0;
    UV_phases           = 0;
    stats               = struct('epsilon',struct(),'UV',struct());
    cost                = struct('epsilon',struct(),'UV',struct());

    recorder            = struct(   'epsilon_fn',   @epsilonPhase,      ...
                                    'uv_fn',        @uvPhase,           ...
                                    'getRecord',    @getRecord          );
                    
    function varargout = epsilonPhase(varargin)
        epsilon_phases          = epsilon_phases + 1;
        varargout               = epsilon_cell_args;
        [varargout{:}]          = epsilon_fn(varargin{:});
        recordPhase(varargout{end},'epsilon');                                              
    end

    function varargout = uvPhase(varargin)
        UV_phases               = UV_phases + 1;
        varargout               = UV_cell_args;
        [varargout{:}]          = uv_fn(varargin{:});
        recordPhase(varargout{end},'UV');
    end

    function recordPhase(phase_recs,type)
        if record_level > 1
            count               = count + 1;
            recs{count}         = addTypeField(phase_recs,type);
        end
        sumData(phase_recs,type);
    end

    function sumData(phase_recs,type)
        stats.(type)    = addStructFields(stats.(type),phase_recs.stats);
        cost.(type)     = addStructFields(cost.(type),phase_recs.cost);
    end

    function rec = getRecord(status,iters)
        cost.total      = cost_counter.getCumulative();
        rec = struct(                                               ...
            'halt_status',      status,                             ...
            'iters',            iters,                              ...
            'epsilon_phases',   epsilon_phases,                     ...
            'UV_phases',        UV_phases,                          ...
            'stats',            stats,                              ...
            'cost',             cost                                );
        if record_level > 1
            rec.phases  = recs(1:count);
        end
    end   
end

function s = addTypeField(s,type)
    s.type  = type;
    n       = numel(fieldnames(s));
    s       = orderfields(s,[n,1:n-1]);
end

function s_result = addStructFields(s_total,s_to_add)
    fields_to_add   = fieldnames(s_to_add);
    
    for j = 1:length(fields_to_add)
        name = fields_to_add{j};
        if isnumeric(s_to_add.(name))
            if isfield(s_total,name)
                s_total.(name) = s_total.(name) + s_to_add.(name);
            else
                s_total.(name) = s_to_add.(name);
            end
        end
    end        
    s_result = s_total;
end
