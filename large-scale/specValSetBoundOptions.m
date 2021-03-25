function opts = specValSetBoundOptions(user_opts)
%   specValSetBoundOptions:
%       Validate user options struct for specValSetBound.m.  If user_opts
%       is [] or not provided, returned opts will be specValSetBound's
%       default parameters.
%
%   USAGE:
%       opts = specValSetBoundOptions();
%       opts = specValSetBoundOptions(user_opts);
%
%   INPUT:
%       user_opts   Struct of settable algorithm parameters.  No fields are 
%                   required, irrelevant fields are ignored, and user_opts 
%                   may be given as [].
%   
%   OUTPUT:
%       opts        Struct of all user-tunable specValSetBound parameters.  
%                   If a field is provided in user_opts, then the user's
%                   value is checked to whether or not it is a valid value,
%                   and if so, it is set in opts.  Otherwise, an error is
%                   thrown.  If a field is not provided in user_opts, opts
%                   will contain the field with its default value.
%
%   PARAMETERS 
%
%   .discrete_time              [logical | {false}]
%       By default, specValSetBound approximates the spectral value set
%       abscissa.  If you wish to instead approximate the spectral value
%       set radius, set .discrete_time to true.
%
%   .real_frobenius_norm        [logical | {false}]                      
%       By default, specValSetBound approximates the complex-valued
%       spectral value set abscissa|radius, that is, for spectral value
%       sets allowing complex perturbations, bounded by the spectral norm.
%       However, when this parameter is set to true, specValSetBound
%       instead approximates the real-valued Frobenius-norm bounded
%       spectral value set abscissa|radius, where only real-valued
%       perturbations are allowed and which are bounded by the Frobenius
%       norm.
%
%   .complex_ode_iteration      [logical | {false}]                      
%       For approximating the complex-valued spectral value set abscissa or
%       radius (i.e. when .real_frobenius_norm is false), specValSetBound
%       has two possible iteration types that can be used.  When this
%       parameter is true, the iteration is based on the SVSAR method of
%       [GGO13].  If system matrix D is also nonzero, this iteration type
%       requires solving a pair of special linear systems on each
%       iteration, which usually has a negligible cost (see
%       .phi_solver_opts for more details).  When this parameter is set to
%       true, an ODE-based iteration, described in [GGMO17, Appendix B in
%       Supplementary Materials], is used instead.  Note that this
%       parameter has no effect when .real_frobenius_norm is set to true.
%
%   .sparse_mode                [integer in {0,1,2} | {2}]
%       Controls whether specValSetBound will use dense or sparse
%       eigenvalue solvers.
%           0 - force dense eigenvalue solves
%           1 - force sparse eigenvalue solves
%           2 - automatically select dense or sparse eigenvalue solves
%               based on how the system matrix A is provided.  If A is
%               given as a sparse matrix, an outer product, or as a
%               function handle, sparse eigenvalues solves will be used.
% 
%   .left_eigenvector           [logical | {false}]     
%       Setting this to true will cause the left eigenvector, corresponding
%       to the eigenvalue that the method has converged to, to also be
%       returned, along with value abs(y'*x).  Note that enabling this may
%       cause an additional eigenvalue solve to be incurred.  If sparse
%       eigensolves are being used, there is a possibility that this
%       additional required solve can fail.  In this case, the left
%       eigenvector and abs(y'*x) will each be returned as [].
%
%   .initial_perturbation       [struct | {initialPerturbationOptions()}]
%       By default, specValSetBound will compute an initial perturbation
%       matrix U*V' based off of the right and left eigenvectors of the
%       rightmost (continuous-time systems) or outermost (discrete-time
%       systems) eigenvalue of system matrix A.  However, the user may
%       optionally choose the kth rightmost/outermost eigenvalue and/or
%       specify their own initial perturbation matrix.  When a perturbation
%       is supplied, the code chooses the kth eigenvalue of the perturbed
%       system matrix, for the user-supplied perturbation U*V'.  For more
%       details, see parameters .kth_eigenvalue, .u, .v from
%       initialPerturbationOptions; note that the .epsilon is ignored, as
%       epsilon is provided as its own bonafide input argument to
%       specValSetBound.
%   
%   .expansion_opts             [struct | {uvExpandOptions()}]
%       This is a struct of the main parameters controlling
%       specValSetBound's algorithm.  For more details, see
%       uvExpandOptions.
%
%   .eig_solver_opts            [struct | {eigSolverOptions()}]
%       This is a struct of the parameters for sparse eigenvalue solves.
%       If dense eigensolves are being used, this struct of parameters is
%       ignored.  For more details, see eigenSolverOptions.
%
%   .phi_solver_opts            [struct | {phiSystemsSolverOptions()}]   
%       If .real_frobenius_norm and .complex_ode_iteration are both false
%       (their defaults), each iteration requires solving a pair of
%       specialized linear systems, which are typically inexpensively
%       solved using various direct methods.  By default, specValSetBound
%       will select the most efficient method automatically.  However, one
%       can also force which one of the available methods is employed via
%       this set of options.  Furthermore, this set of options allows one
%       to set pcg's options, when it is being used to iteratively solve
%       the systems.  For more details, see phiSystemsSolverOptions.
%
%   .count_multiplies           [logical | {false}]
%       When this is set to true, specValSetBound will additionally return
%       the total number of matrix multiplications for A,A',B,B',C,C',D,D'.
%       Note enabling this adds some overhead.
%      
%   .record_level               [integer in {0,1,2,3} | {0}]
%       Controls how much metadata is included getStabRadBound's second
%       output argument.  Increasing this value increases memory and
%       computation overhead.
%           0 - only basic metadata regarding the total incurred costs
%           1 - breaks down the metadata on per phase basis, that is, for  
%               initialization and the expansion phase
%           2 - adds a history of the accepted iterates 
%           3 - additionally includes the rejected iterates as well 
%               (from line searches, rejected extrapolations, etc).
%       Note that this parameter has no effect if specValSetBound is called
%       with only one output argument.
%
%   .record_UV                  [logical | {false}]
%       If this is true, the pair of vectors/matrices for each perturbation
%       U*V' will be additionally included in the recorded history that is
%       returned in specValSetBound's second output argument.  Setting this
%       value to true can significantly increase memory consumption.  Note
%       that this parameter has no effect when specValSetBound is called
%       with only one output argument or if opts.record_level == 0.
%       
%   .print_level                [integer in {0,1,2,3} | {1}]
%       Level of detail printed to console regarding optimization progress:
%           0 - no printing whatsoever
%           1 - prints info for each accepted step 
%           2 - additionally prints rejected steps, which may be either 
%               extrapolations, interpolations, or if the line search fails
%               (the last of which can only happens on the very last step)
%           3 - additionally prints all line search evaluations.
%
%   .print_info_msg             [logical | {true}]
%       If .discrete_time is true, specValSetBound will print out an
%       message box in orange text, warning the user that the iteration
%       printing shows the progression of the root function for determininfg
%       stability, namely the spectral value set radius MINUS 1.  To
%       disable this message, set this parameter to false.
%
%   .print_use_orange           [logical | {true}]
%       By default, specValSetBound makes selected use of an undocumented
%       MATLAB feature to enable printing certain text on the console in
%       orange (without issusing a warning).  However, the user is the
%       given option to disable it, since support cannot be guaranteed
%       (since it is an undocumented feature).
%
%   .print_ascii                [logical | {false}]
%       By default, specValSetBound's printed output uses the extended
%       character map, so nice looking tables can be made.  However,
%       diary() does not capture these symbols.  So, if you need to record
%       the output, you can restrict the printed output to only use the
%       basic ASCII character map, which may look better when captured by
%       diary().
%       Note: on Windows, opts.print_ascii is forced to true as the
%       extended ASCII used by MATLAB on WIndows are not monospaced.
%
%   See also initialPerturbationOptions, uvExpandOptions,
%   eigenSolverOptions, and phiSystemsSolverOptions.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   specValSetBoundOptions.m introduced in ROSTAPACK Version 1.0
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
        validator = optionValidator(    'specValSetBound',              ...
                                        defaults,sub_validators         );
    end
    
    if nargin < 1 || isempty(user_opts)
        % We must remove the unused epsilon field here, and not above, as
        % the default parameters are populated into each validated 
        % substruct of options by the validator object.   
        opts = removeUnusedFields(validator.getDefaultOpts());
        return
    end
   
    % Also remove the epsilon field here, otherwise the subvalidator may
    % assert it value (instead of ignoring it) 
    validator.setUserOpts(removeUnusedFields(user_opts));
    
    try 
        validator.setLogical('discrete_time');
        validator.setLogical('real_frobenius_norm');
        validator.setLogical('complex_ode_iteration');
        validator.setIntegerInRange('sparse_mode',0,2);
        validator.setLogical('left_eigenvector');
        
        validator.setStructWithValidation('initial_perturbation');
        validator.setStructWithValidation('expansion_opts');
        validator.setStructWithValidation('eig_solver_opts');
        validator.setStructWithValidation('phi_solver_opts');
                
        validator.setLogical('count_multiplies');
        validator.setIntegerInRange('record_level',0,3);
        validator.setLogical('record_UV');
        validator.setIntegerInRange('print_level',0,3);
        validator.setLogical('print_info_msg');
        validator.setLogical('print_use_orange');
        validator.setLogical('print_ascii');
       
        opts = removeUnusedFields(validator.getValidatedOpts());
        
        % Extended ASCII chars in MATLAB on Windowsare not monospaced so
        % don't support them.
        if ~opts.print_ascii
            validator.assert(~ispc(),                                   ...
                'only opts.print_ascii == true is supported on Windows.');
        end
  
    catch err
        err.throwAsCaller();
    end
    
    function opts = removeUnusedFields(opts)
        % Since opts.initial_perturbation.epsilon is ignored (as epsilon is
        % explicitly provided to specValSetBound as its own bonafide input
        % argument), let's completely remove all traces of this unnecessary
        % parameter field, so users aren't tempted to use it.
        if isfield(opts,'initial_perturbation')
            if isfield(opts.initial_perturbation,'epsilon')
                init_pert = rmfield(opts.initial_perturbation,'epsilon');
                opts.initial_perturbation = init_pert;
            end
        end
    end
end



function [default_opts,sub_validators] = getDefaults()
    default_opts = struct(                                              ...
        'discrete_time',        false,                                  ...
        'real_frobenius_norm',  false,                                  ...
        'complex_ode_iteration',false,                                  ...
        'sparse_mode',          2,                                      ...
        'left_eigenvector',     false,                                  ...
        'initial_perturbation', [],                                     ...
        'expansion_opts',       [],                                     ...
        'eig_solver_opts',      [],                                     ...
        'phi_solver_opts',      [],                                     ...
        'count_multiplies',     false,                                  ...
        'record_level',         0,                                      ...
        'record_UV',            false,                                  ...
        'print_level',          1,                                      ...
        'print_info_msg',       true,                                   ...
        'print_use_orange',     true,                                   ...
        'print_ascii',          false || ispc()                         );
    
    sub_validators = struct(                                            ...
        'initial_perturbation', @initialPerturbationOptions,            ...
        'expansion_opts',       @uvExpandOptions,                       ...
        'eig_solver_opts',      @eigenSolverOptions,                    ...
        'phi_solver_opts',      @phiSystemsSolverOptions                );                       
end
