function opts = upperBoundOptions(user_opts)
%   upperBoundOptions:
%       Validate user options struct for parameters specific to finding an
%       initial destabilizing perturbation using Algorithm Fast Upper Bound
%       [MO16, Section 4.4].  If user_opts is [] or not provided, returned
%       opts will be the default parameters.
%
%   USAGE:
%       opts = upperBoundOptions();
%       opts = upperBoundOptions(user_opts);
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
%   .maxit                      [positive integer | {100}]
%       Maximum number of iterations that are allowed.  Each iteration
%       consists of two phases: expanding via updating epsilon, and
%       expanding via updating the perturbation U*V' (which happens first
%       is configurable).
% 
%   .rel_diff_tol               [positive finite real | {1e-12}]
%       This is the convergence tolerance for the two phases.  When a phase
%       satisfies this tolerance, it is an indication that it can no longer
%       progress and the other phase should be attempted to make further
%       progress.
% 
%   .rel_step_size_tol          [value in [0,1) | {0.1}]
%       Either phase will halt if its current step size is smaller than
%       opts.rel_step_size_tol times the largest step length encountered
%       during this phase.  Once a phase seems to make slow progress
%       relative to some previous large step it has taken, the current
%       phase will halted and the other phase will take over.  Setting this
%       to higher values encourages more frequent switches, with zero
%       completely disabling this feature.  This parameter has no effect
%       when both opts.uv_steps_per_iter and opts.epsilon_steps_per_iter
%       are set to one.
%       
%   .upperbound_quality         [value in [0,1] | {0.99}]
%       When this is zero, the hybrid expansion-contraction will commence
%       immediately once the first destabilizing perturbation has been
%       found.  However, by setting this value to higher numbers, the
%       algorithm will first spend more effort to be better locate the a
%       locally maximizing frequency, which often lead to better
%       approximation quality.
%
%   .update_UV_first            [logical | {false}]
%       This chooses whether or not expansion towards finding a
%       destabilizing perturbation should begin with updating the U*V'
%       perturbation or by increasing epsilon.
%
%   .UV_steps_per_iter          [positive integer | {1}]
%       The maximum number of U*V' perturbation updates allowed per U*V'
%       expansion phase.
%
%   .epsilon_steps_per_iter     [positive integer | {1}]
%       The maximum number of epsilon updates per epsilon expansion phase.
%
%   .epsilon_step_multiplier    [value in [1,inf) | {2}]
%       When increasing epsilon, the Newton step towards the stability
%       boundary will be multiplied by this factor. 
%
%   .epsilon_limit_factor       [value in [0,1] | {0.99}]
%       Once epsilon has exceeded this fraction of its upperbounded, i.e.
%       this value times 1/norm(D), remaining increases in epsilon will 
%       only happen by bisection steps between its current value and
%       1/norm(D).  This parameter has no effect when D is zero.
%
%   .epsilon_line_search_opts   [struct of line search parameters]
%       
%       .maxit                  [positive integer | {10}]
%           Maximum number of line search evaluations (including t = 0)
%           that are allowed.
% 
%       .df_tol                 [nonnegative finite real | {0}]
%           If the derivative of the line search function with respect to t
%           at t = 0 is below this value, don't bother to do the line
%           search, including evaluating the full step at t = 1,
% 
%       .t_tol                  [nonnegative finite real | {1e-14}]
%           Terminate the line search when:
%               abs(t_k - t_{k-1}) <= opts.t_tol 
%           holds.
% 
%       .model                  [value in {1,2,3} | {1}]
%           1 - bisection line search
%           2 - quadratic interpolation line search
%           3 - cubic interpolation line search
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   upperBoundOptions.m introduced in ROSTAPACK Version 1.0   
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
        [defaults,sub_validators] = getDefaults();
        validator = optionValidator(                                    ...
            'ROSTAPACK',defaults,sub_validators,'upperbound_opts'       );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
   
    validator.setUserOpts(user_opts);
    
    try
        validator.setIntegerPositive('maxit');
        validator.setRealInIntervalOO('rel_diff_tol',0,inf);
        validator.setRealInIntervalCO('rel_step_size_tol',0,1);
        validator.setRealInIntervalCC('upperbound_quality',0,1);
        validator.setLogical('update_UV_first');
        validator.setIntegerPositive('UV_steps_per_iter');
        validator.setIntegerPositive('epsilon_steps_per_iter');
        validator.setRealInIntervalCO('epsilon_step_multiplier',1,inf);
        validator.setRealInIntervalCC('epsilon_limit_fraction',0,1);
            
        validator.setStructWithValidation('epsilon_line_search_opts');

    catch err
        err.throwAsCaller();
    end
   
    opts = validator.getValidatedOpts(); 
end

function [default_opts,sub_validators] = getDefaults()
    default_opts = struct(                                          ...
        'maxit',                        100,                        ...
        'rel_diff_tol',                 1e-12,                      ...
        'rel_step_size_tol',            0.1,                        ...
        'upperbound_quality',           0.99,                       ...
        'update_UV_first',              false,                      ...
        'UV_steps_per_iter',            1,                          ...
        'epsilon_steps_per_iter',       1,                          ...
        'epsilon_step_multiplier',      2,                          ...
        'epsilon_limit_fraction',       0.99,                       ...
        'epsilon_line_search_opts',     []                          );
    sub_validators = struct(                                        ...
         'epsilon_line_search_opts',    @epsilonLineSearchOptions   );
end

function opts = epsilonLineSearchOptions(user_opts)

    persistent validator;
    
    if isempty(validator)
        validator = optionValidator(                                    ...
            'ROSTAPACK', getLineSearchDefaults(), [],                   ...
            'upperbound_opts.epsilon_line_search_opts'                  );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
  
    validator.setUserOpts(user_opts);
        
    try
        validator.setIntegerPositive('maxit');
        validator.setRealInIntervalCO('df_tol',0,inf);
        validator.setRealInIntervalCO('t_tol',0,inf);
        validator.setIntegerInRange('model',1,3);
    catch err
        err.throwAsCaller();
    end
    
    opts = validator.getValidatedOpts();
end

function default_opts = getLineSearchDefaults()
    default_opts = struct(  'maxit',                10,                 ...
                            'df_tol',               0,                  ...
                            't_tol',                1e-14,              ...
                            'model',                1                   );
end
