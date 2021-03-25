function obj = initializePerturbationRecord(sys_c,opts)
%   initializePerturbationRecord:
%       This is an "object" for recordings stats and the iterate history of
%       initializePerturbation.
% 
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, with no perturbation ever set.  Note
%           that this is only used to get the data for the total costs and
%           the correctly-dimensioned zero perturbation U*V'.  Iterate info
%           must be recorded by passing in the explicit details.
%
%       opts            [required: struct of parameters]
%           A required struct settable parameters necessary to run this
%           routine. 
%
%           .record_level       [value in {0,1}]
%               Determines how much metadata is recorded:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the iterates. 
%
%           .record_UV          [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be recorded.
%
%   OUTPUT:
%       recorder
%           The recording "object", a struct with the following methods:
%       
%       .addAStability(f,z)
%           Add the info about the rightmost/outermost eigenvalue of A to
%           the record, where f is the value of the root function for z.
%           [This method is only present if opts.record_level > 0]
% 
%       .addAEigenvalue(f,z,kth_eigenvalue)
%           Add the info about the kth rightmost/outermost eigenvalue of A
%           to the record, where f is the value of the root function for z.
%           [This method is only present if opts.record_level > 0]
%
%       .addPerturbation(f,z,epsilon,U,V,kth_eigenvalue)  
%           Add the info about the kth eigenvalue of the initial
%           perturbation epsilon*U*V', where f is the value of the root
%           function for z. 
%           [This method is only present if opts.record_level > 0]
%
%       rec = recorder.getRecord()
%           A cell array containing the recorded info, where the odd
%           entries are the labels and each even entry is the data
%           associated with the preceding label in the cell array.
%
%   See also initializePerturbation.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   initializePerturbationRecord.m introduced in ROSTAPACK Version 1.0
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
    
    record_level        = opts.record_level;
    switch record_level
        case 0
            obj         = struct('getRecord',@getRecord);
            return
        case 1
            getStructFn = ternOp(opts.record_UV,@infoUV,@infoMinimal);
        otherwise
            error('Unrecognized value of opts.record_level!');
    end

    % data structure for storing all recorded iterates (only ever need at
    % most 3 entries for initializePerturbation)
    rec             = cell(3);
    count           = 0;
    
    [U_zero,V_zero] = sys_c.getUV();
    
    obj             = struct(   'addAStability',    @addAStability,     ...
                                'addAEigenvalue',   @addAEigenvalue,    ...
                                'addPerturbation',  @addPerturbation,   ...
                                'getRecord',        @getRecord          );
                            
    function addAStability(f,z)
        addPoint('A stability',f,z,0,U_zero,V_zero,1);
    end

    function addAEigenvalue(f,z,kth_eigenvalue)
        addPoint('A eigenvalue',f,z,0,U_zero,V_zero,kth_eigenvalue); 
    end

    function addPerturbation(f,z,epsilon,U,V,kth_eigenvalue)
        addPoint('Initial Perturbation',f,z,epsilon,U,V,kth_eigenvalue);
    end

    function addPoint(type,f,z,epsilon,U,V,kth_eigval)
        cost        = cost_counter.getLatest();
        s           = getStructFn(type,f,z,epsilon,U,V,kth_eigval,cost);
        count       = count + 1;
        rec{count}  = s;
    end

    function r = getRecord()
        r           = {'cost',cost_counter.getCumulative()};
        if record_level > 0
            r       = [r {'iterates',cell2mat(rec(1:count))}];
        end
    end
end

function s = infoMinimal(type,f,z,epsilon,~,~,kth_eigenvalue,cost)
    s = struct( 'type',             type,               ...
                'f',                f,                  ...
                'z',                z,                  ...
                'epsilon',          epsilon,            ...  
                'kth_eigenvalue',   kth_eigenvalue,     ...
                'cost',             cost                );
end

function s = infoUV(type,f,z,epsilon,U,V,kth_eigenvalue,cost)
    s = struct( 'type',             type,               ...
                'f',                f,                  ...
                'z',                z,                  ...
                'epsilon',          epsilon,            ...  
                'U',                U,                  ...
                'V',                V,                  ...
                'kth_eigenvalue',   kth_eigenvalue,     ...
                'cost',             cost                );
end
