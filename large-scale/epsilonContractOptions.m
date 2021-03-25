function opts = epsilonContractOptions(user_opts)
%   epsilonContractOptions:
%       Validate a struct of parameters pertaining to the contraction phase
%       of the Hybrid Expansion-Contraction (HEC) algorithm.  If user_opts
%       is [] or not provided, returned opts will be the default
%       parameters.
%
%   USAGE:
%       opts = epsilonContractOptions();
%       opts = epsilonContractOptions(user_opts);
%
%   INPUT:
%       user_opts   Struct of settable algorithm parameters.  No fields are 
%                   required, irrelevant fields are ignored, and user_opts 
%                   may be given as [].
%   
%   OUTPUT:
%       opts        Struct of all user-tunable parameters.  If a field is
%                   provided in user_opts, then the user's value is checked
%                   to whether or not it is a valid value, and if so, it is
%                   set in opts. Otherwise, an error is thrown.  If a field
%                   is not provided in user_opts, opts will contain the
%                   field with the default value.
%
%   PARAMETERS 
%
%   .maxit                      [positive integer | {10}]
%       If the contraction phase has at least made some acceptable
%       reduction in the value of epsilon, such that the root function is
%       still positive, then the contraction phase will terminate after 
%       at most opts.maxit iterations have occurred.  If no progress had
%       been made, the contraction phase will continue iterating beyond
%       opts.maxit iterations, until one of its other termination
%       conditions occurs:
%           - converged
%           - relative reduction target achieved
%           - exhausted floating point precision
% 
%   .tol                        [positive finite real | {1e-12}]
%       The algorithm has converged when the value of the root function f
%       is in (0,opts.tol).  
% 
%   .rel_reduction_tol          [value in [0,1) | {0.01}]
%       The algorithm will halt early if the value of the root function f
%       is less than opts.rel_reduction_tol * f_initial (the initial value
%       of the root function before starting the contraction phase).
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   epsilonContractOptions.m introduced in ROSTAPACK Version 1.0  
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

    persistent validator;
    
    if isempty(validator)
        validator = optionValidator(                                    ...
            'ROSTAPACK', getDefaults(), [], 'contraction_opts'          );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
   
    validator.setUserOpts(user_opts);
    
    try
        validator.setIntegerPositive('maxit');
        validator.setRealInIntervalOO('tol',0,inf);
        validator.setRealInIntervalCO('rel_reduction_tol',0,1);
    catch err
        err.throwAsCaller();
    end
    
    opts = validator.getValidatedOpts();
end

function default_opts = getDefaults()
    default_opts = struct(  'maxit',                10,         ...
                            'tol',                  1e-12,      ...
                            'rel_reduction_tol',    0.01        );          
end
