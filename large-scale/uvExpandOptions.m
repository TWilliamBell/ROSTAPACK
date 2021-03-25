function opts = uvExpandOptions(user_opts,varargin)
%   uvExpandOptions:
%       Validate user options struct for parameters specific to finding
%       approximations of the spectral value set abscissa or radius, which
%       is the expansion phase of the Hybrid Expansion-Contaction (HEC)
%       algorithm, implemented by getStabRadBound.  If user_opts is [] or
%       not provided, returned opts will be the default parameters.
%       
%       Note that the interpolation, extrapolation, and line search
%       parameters are also used by the upper bound procedure of
%       getStabRadBound.
%
%   USAGE:
%       opts = uvExpandOptions();
%       opts = uvExpandOptions(user_opts);
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
%   .maxit                      [positive integer | {1000}]
%       Maximum number of iterations that are allowed.  
% 
%   .rel_diff_tol               [positive finite real | {1e-12}]
%       This is the convergence tolerance.  
% 
%   .rel_step_size_tol          [value in [0,1) | {0.1}]
%       The method will halt expanding if its current step size becomes
%       smaller than opts.rel_step_size_tol times the largest step length
%       encountered so far.  In the context of the HEC algorithm, this
%       allows the expansion phases to save cost by only doing complete
%       expansions once HEC has nearly converged.  In theory, this means
%       that HEC's quadratic convergence rate is downgraded from quadratic
%       to superlinear but it is generally a net win due to the dramatic
%       reduction of iterations incurred inside the expansion phases.
%      
%   .interp_threshold           [value in [1,inf] | {1.5}]
%       A predicted interpolation step has to be at least this times better
%       than the full step in order to attempt taking the interpolated
%       step.  When this is set to inf, interpolation is disabled. 
%
%   .extrap_size                [positive integer | {5}]
%       Use the k most recent perturbations to attempt an extrapolation.
%       When this is set to 1, extrapolation is disabled. 
%
%   .extrap_min_size            [integer in {2,...,.extrap_size} | {5}]
%       Allow an extrapolation to be attempted if the method hits its
%       maximum iteration count and there are at least this many previous
%       perturbations available to use for an extrapolation, even if there
%       are less than opts.extrap_size.  
%
%   .extrap_rollover            [nonnegative integer | {0}]
%       The number of the most recent perturbations to be reused for the
%       next extrapolation attempt.  By default, when an extrapolation is
%       attempted, the cache of previous perturbations is completely
%       cleared.  
%   
%   .extrap_initial_skip        [nonnegative integer | {2}]
%       The number of initial perturbations to skip before accruing
%       perturbations for extrapolation purposes.  The initial steps are
%       typically quite large and can thus be poor choices to include in
%       the sequence to extrapolated.
% 
%   .line_search_opts           [struct of line search parameters]
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
%       .UV_delta               [value in (0,1) | {0.999999}]
%           If the current perturbation and the full update step are
%           identical, except for their signs, then interpolating them with
%           t = 0.5 will be undefined (since it will a zero matrix that
%           must have unit norm).  Whenever a value of t would cause a
%           complete annihilation, t is instead multiplied by this number
%           in order to have the line search produce a nonzero
%           perturbation matrix.
%
%       .UV_dual_mode           [logical or value in (0,inf) | {1e-12}]
%           When this is true, dual line searches are done in parallel in
%           both "directions", since we cannot always trust the sign of the
%           derivative of the line search function at t = 0, if its
%           magnitude is small.  When this is false, the sign of the
%           derivative is always trusted and a regular line search is done.
%           When this is a positive value, the dual line search is
%           adaptively invoked whenever the magnitude of the line search
%           derivative falls below this value. 
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   uvExpandOptions.m introduced in ROSTAPACK Version 1.0
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
            'ROSTAPACK',defaults,sub_validators,'expansion_opts'        );
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
        validator.setRealInIntervalCC('interp_threshold',1,inf);
        
        maxit       = validator.getValue('maxit');
        validator.setIntegerInRange('extrap_size',1,maxit);
        extrap_size = validator.getValue('extrap_size');
        
        if extrap_size > 1
            validator.setInteger('extrap_min_size');
            validator.setInteger('extrap_initial_skip');
            validator.setInteger('extrap_rollover');
             
            min_size    = validator.getValue('extrap_min_size');
            initial     = validator.getValue('extrap_initial_skip');
            rollover    = validator.getValue('extrap_rollover');
            
            validator.assert(min_size >= 2 && min_size <= extrap_size,  ...
                'extrap_min_size must be in {2,...,extrap_size}.'       );
            
            validator.assert(initial >= 0 && initial < max(3,maxit-1),  ...
                'extrap_initial_skip must be in {0,...,max(3,maxit-2)}.');
           
            validator.assert(rollover >= 0 && rollover < extrap_size,   ...
                'extrap_rollover must be in {0,...,extrap_size-1}.'     );
        end
        
        validator.setStructWithValidation('line_search_opts');
    catch err
        err.throwAsCaller();
    end
   
    opts = validator.getValidatedOpts(); 
end

function [default_opts,sub_validators] = getDefaults()
    default_opts = struct(                              ...
        'maxit',                1000,                   ...
        'rel_diff_tol',         1e-12,                  ...
        'rel_step_size_tol',    0,                      ...
        'interp_threshold',     1.5,                    ...
        'extrap_size',          5,                      ...
        'extrap_min_size',      5,                      ...
        'extrap_rollover',      0,                      ...
        'extrap_initial_skip',  2,                      ...
        'line_search_opts',     []                      );
    sub_validators = struct(                            ...
        'line_search_opts',     @uvLineSearchOptions    );
end

function opts = uvLineSearchOptions(user_opts)

    persistent validator;
    
    if isempty(validator)
        validator = optionValidator(                                    ...
            'ROSTAPACK', getLineSearchDefaults(), [],                   ...
            'expansion_opts.line_search_opts'                           );
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
        validator.setRealInIntervalOO('UV_delta',0,1);
        
        if validator.isSpecified('UV_dual_mode')
            dm = user_opts.UV_dual_mode;
            validator.assert(   islogical(dm) || (dm > 0 && ~isinf(dm)),...
                                dualModeMsg()                           );
        end
    catch err
        err.throwAsCaller();
    end
    
    opts = validator.getValidatedOpts();
end

function default_opts = getLineSearchDefaults()
    default_opts = struct(  'maxit',                10,                 ...
                            'df_tol',               0,                  ...
                            't_tol',                1e-14,              ...
                            'model',                1,                  ...
                            'UV_delta',             0.999999,           ...
                            'UV_dual_mode',         false               );
end

function m = dualModeMsg()
m = sprintf([                                                           ...
'UV_dual_mode must either be:\n'                                        ...
'- a logical specifying whether the dual line search is on/off\n'       ...
'- a positive finite real value for the adaptive dual line search.'     ]);
end
