function opts = specValSetOptions(user_opts,varargin)
 %   specValSetOptions:
%       Validate user options struct for specValSet.m.  If user_opts is []
%       or not provided, returned opts will be specValSet's default
%       parameters.
%
%   USAGE:
%       opts = specValSetOptions();
%       opts = specValSetOptions(user_opts);
%
%   INPUT:
%       user_opts   Struct of settable algorithm parameters.  No fields are 
%                   required, irrelevant fields are ignored, and user_opts 
%                   may be given as [].
%   
%   OUTPUT:
%       opts        Struct of all user-tunable specValSet parameters.  If a
%                   field is provided in user_opts, then the user's value
%                   is checked whether or not it is a valid value, and if
%                   so, it is set in opts.  Otherwise, an error is thrown.
%                   If a field is not provided in user_opts, opts will
%                   contain the field with its default value.
%
%   BASIC PARAMETERS 
%
%   .discrete_time              [logical | {false}]
%       By default, specValSet computes the spectral value set
%       abscissa.  If you wish to instead compute the spectral value
%       set radius, set .discrete_time to true.
%
%   .ignore_unperturbable       [logical | {true}]                      
%       By default, specValSet computes the epsilon spectral value set
%       abscissa|radius by initializing at a globably rightmost|outermost
%       eigenvalue of (A,E) of those that are also both controllable and
%       observable.  By setting this option to false, the method will be
%       initialized at a globally rightmost|outermost eigenvalue of (A,E),
%       regardless of whether it is controllable and/or observable.  If
%       there exists a non-isolated point in the spectral value set that is
%       to the right (or further out) than this starting eigenvalue, the
%       computed result should be unchanged.  Otherwise, the method may
%       terminate at the initial eigenvalue if it is uncontrollable or
%       unobservable.  Note that an eigenvalue is considered controllable
%       and observable if norm(B'*y) and norm(C*x) are respectively greater
%       than zero (by any amount, no tolerance is used), where x and y are
%       respectively the right and left eigenvectors.
%
%   .warm_start                 [complex scalar | {[]}]                      
%       The method can be warm started at a point other than than a
%       rightmost|outermost eigenvalue of (A,E).  In order for a warm start
%       to help accelerate the computation, the provided point must be in
%       the spectral value set and it should provide a better initial
%       estimate of the spectral value set abscissa|radius, i.e. better
%       than the spectral abscissa|radius.  Warm starts can be particularly
%       useful when computing the spectral value set abscissa|radius for an
%       increasing sequence of epsilon values; each computed point that
%       attains the abscissa|radius for one value of epsilon can be used
%       to warm start the computation for the next larger epsilon value.
%
%   .maxit                      [positive integer | {100}]
%       Maximum number of allowed iterations.  As the method has local
%       quadratic convergence, most problems only require a handful of
%       iterations.
%
%   .tol                        [positive finite real | {1e-14}]
%       The desired relative accuracy of the computed the spectral value
%       set abscissa|radius.  Note if opts.fast_search is disabled, the
%       computation may have significantly worse accuracy than indicated by
%       this tolerance.  When opts.fast_search, there is little benefit to
%       setting this significantly higher than machine precision.
%
%   .suppress_warnings          [logical | {true}] 
%       By default, specValSet turns off all warnings while it is running
%       and restores the warning state when it terminates.  The reason is
%       that when opts.fast_search is true, the LUs and backsolves used in
%       specValSet may often throw warnings.  However, since specValSet is
%       already designed to handle the causes of these warnings, there is
%       no need to ever see them on the console.
%
%   .plotting                   [logical | {false}] 
%       Enabling this will cause specValSet to open a new figure depicting:
%           The eigenvalues of (A,E)
%           - controllable and observable:      black dots
%           - uncontrollable or unobservable:   green dots
%           Initial eigenvalue:                 inside large red square
%           Warm start point (if interior):     blue square 
%           Horizontal/radial searches:         black dashed segments
%           Horizontal/radial start points:     black o's
%           Vertical/circular searches:         red vertical segments/arcs
%           Computed boundary points:           red o's
%           Cross section midpoints:            red x's 
%           Point attaining abscissa|radius:    blue o
%       Note that the spectral value set boundary is not plotted.
%
%   ADVANCED PARAMETERS: in general, these should be left at their defaults
%
%   .fast_search                [logical | {true}]   
%       Use the faster and more accurate horizontal|radial searches based
%       on local root-finding.  Generally speaking, this should not be
%       disabled.  
%
%   .vertical_search_first      [logical | {false}]   
%       For the abscissa case only.  By default, the method starts with a a
%       horizontal search instead of a vertical one.  Starting with a
%       vertical search may cause additional searches to be incurred before
%       convergence, which can significantly increase the run time.  This
%       should not be changed from its default value.
%
%   .bracket_order              [1 or 2 | {2}]
%       Only relevant when opts.fast_search is true.  By default,
%       second-order (Halley) steps are used to first bracket any root
%       finding problem.  Optionally, one may instead set this to 1 to only
%       use first-order (Newton) steps for bracketing, which can be a bit
%       cheaper to compute but may incur more steps until a bracket is
%       found.  For most systems, this should not be changed but some
%       problems with relatively large m,p, using the first-order
%       bracketing could be faster.  
%
%   .root_order                 [any value in [1,2] | {2}]
%       Only relevant when opts.fast_search is true.  By default, a
%       second-order Halley method with bracketing is used to find roots,
%       which should have a cubic rate of convergence.  By setting this to
%       1, a first-order Newton method with bracketing will be used, which
%       should have a quadratic rate of convergence but each iteration is a
%       bit cheaper to compute.  If this option is set to any number in
%       (1,2), an interpolation-based root-finding method, also involving
%       only first-order terms, is used; the exact choice in (1,2) has no
%       effect on the interpolation method.  For most systems, this should
%       not be changed but some problems with relatively large m,p, using
%       one of these two first-order methods could be faster.
%
%   .solve_type                 [any value in [0,1,2] | {0}]                      
%       Sets how inverses of zE-A should be applied:
%           0: via upper triangular Hessenberg factorization 
%           1: via LU
%           2: via LU with permutations 
%       0 is generally recommended for small-scale systems.
%
%   .force_two_norm             [logical | {false}]                      
%       Only relevant when opts.fast_search is true.  If B=C=I, D=0, and
%       n=m=p, by default the fast searches will be done using an alternate
%       more efficient formulation that only needs zE-A, and not its
%       inverse, thus avoiding the expense of LUs and backsolves.  As a
%       result, in this case there is also little to no reason not use the
%       second-order bracketing and root finding options.  However, as this
%       alterate form involves computing the smallest singular value,
%       numerical accuracy can be adversely affected for some very poorly
%       scaled matrices (due to current limitations of the GESDD routine in
%       LAPACK, which will be called by svd in MATLAB). In this case, one
%       can set this option to true to force the computation involving
%       (zE-A)^{-1} (where then it could be computationally beneficial to
%       consider the first-order bracketing abd root finding options).
%   
%   .random_directions          [positive integer | {3}]
%       Only relevant when opts.discrete_time is true.  This sets the
%       number of random radial search directions to be attempted in order
%       to overcome a possibly (nearly) singular pencil or problematic 
%       interior search.  Three seems to be a good minimum choice. 
%
%   .safeguard_width            [real value in [0,1] | {0.75}]
%       Due to rounding errors when calculating double eigenvalues, a
%       safeguard is needed, which overcomes this difficulty.  The problem
%       is when a previous horizontal/radial search corresponds to a double
%       eigenvalue but at a nonlocally rightmost/outermost boundary point.
%       In this case, the next vertical/circular should find two cross
%       sections that are adjacent at this double eigenvalue, but the
%       rounding errors may cause only the union of the two cross sections
%       to be detected.  In this case, the algorithm would stagnate but the
%       safeguard corrects this problem by splitting any cross section into
%       two sections whenever its midpoint is too close to the previous
%       horizontal/radial search point.  The downside is that the splitting
%       means an additional horizontal/radial search may be incurred, even
%       when the safeguard is not actually needed.  When opts.fast_search
%       is enabled (the default and recommended setting), this additional
%       cost is negligible so the default is to be liberal in invoking the
%       safeguard.  When this is set to zero, the safeguard is completely
%       disabled (not recommended).  When set to a number in (0,1], the
%       safeguard is invoked whenever the previous horizontal/radial search
%       is within the middle opts.safeguard_width fraction of said interval,
%       with 0.5 meaning the middle half of the interval, 1 being the
%       entire interval (except its endpoints), etc.
% 
%   See also specValSet.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   specValSetBoundOptions.m introduced in ROSTAPACK Version 2.0.
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

    persistent validator;
    
    if isempty(validator)
        validator   = optionValidator('specValSet',getDefaults());
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
 
    validator.setUserOpts(user_opts);
    
    try
        % System options
        validator.setLogical('discrete_time');
        
        % Starting options
        validator.setLogical('ignore_unperturbable');
        if validator.isSpecified('warm_start')
            validator.setDimensioned('warm_start',1,1);
            validator.setFiniteValued('warm_start');
        end
               
        % Termination options
        validator.setIntegerPositive('maxit');
        validator.setRealInIntervalOO('tol',0,inf);
        
        % Performance options
        validator.setLogical('fast_search');
        validator.setLogical('vertical_search_first');
        validator.setIntegerInRange('bracket_order',1,2);
        validator.setRealInIntervalCC('root_order',1,2);
        validator.setIntegerInRange('solve_type',0,2);
        validator.setLogical('force_two_norm');
        validator.setIntegerInRange('random_directions',1,inf);
        validator.setRealInIntervalCC('safeguard_width',0,1);
        
        % Plotting / output options
        validator.setLogical('suppress_warnings');
        validator.setLogical('plotting');   
    catch err
        err.throwAsCaller();
    end
   
    opts = validator.getValidatedOpts(); 
end

function default_opts = getDefaults()
    default_opts = struct(                              ...
        'discrete_time',            false,              ...
        'ignore_unperturbable',     true,               ...
        'warm_start',               [],                 ...
        'maxit',                    100,                ...
        'tol',                      1e-14,              ...
        'fast_search',              true,               ...
        'vertical_search_first',    false,              ...
        'bracket_order',            2,                  ...
        'root_order',               2,                  ...
        'solve_type',               0,                  ...
        'force_two_norm',           false,              ...
        'random_directions',        3,                  ...
        'safeguard_width',          0.75,               ...
        'suppress_warnings',        true,               ...
        'plotting',                 false               );
end