function opts = getStabRadBoundOptions(user_opts)
%   getStabRadBoundOptions:
%       Validate user options struct for getStabRadBound.m.  If user_opts
%       is [] or not provided, returned opts will be getStabRadBound's
%       default parameters.
%
%   USAGE:
%       opts = getStabRadBoundOptions();
%       opts = getStabRadBoundOptions(user_opts);
%
%   INPUT:
%       user_opts   Struct of settable algorithm parameters.  No fields are 
%                   required, irrelevant fields are ignored, and user_opts 
%                   may be given as [].
%   
%   OUTPUT:
%       opts        Struct of all user-tunable getStabRadBound parameters.
%                   If a field is provided in user_opts, then the user's
%                   value is checked to whether or not it is a valid value,
%                   and if so, it is set in opts.  Otherwise, an error is
%                   thrown.  If a field is not provided in user_opts, opts
%                   will contain the field with its default value.
%
%   PARAMETERS 
%
%   .maxit                      [positive integer | 100]
%       Maximum number of allowed iterations in the convergent phase, which
%       is the Hybrid Expansion-Contraction (HEC) Algorithm [MO16].
%
%   .discrete_time              [logical | {false}]
%       By default, getStabRadBound approximates the complex or real
%       stability radius for continous-time systems.  If the system is a 
%       discrete-time system, set this option to true.
%
%   .real_frobenius_norm        [logical | {false}]                      
%       By default, getStabRadBound approximates the complex stability
%       radius, the stability radius under complex-valued perturbations,
%       and whose reciprocal is the H-infinity norm.  However, by setting
%       this option to true, getStabRadBound will instead approximate the
%       real stability radius, where only real-valued perturbations are
%       allowed and the perturbations are bounded by the Frobenius norm.
%
%   .complex_ode_iteration      [logical | {false}]                      
%       When approximating the complex stability radius, (i.e. when
%       .real_frobenius_norm is false), getStabRadBound has two possible
%       iteration types that can be used for the expansion subphase of the
%       Hybrid Expansion-Contraction Algorithm [MO16].  When this option is
%       true, the iteration is based on the SVSA/SVSR method of [GGO13]. If
%       system matrix D is also nonzero, this iteration type requires
%       solving a pair of special linear systems on each iteration, but the
%       cost to solve them is usually negligible (for more details, see
%       .phi_solver_opts).  When this parameter is set to true, an
%       ODE-based iteration, described in [GGMO17, Appendix B in
%       Supplementary Materials], is used instead.  Note that this
%       parameter has no effect when .real_frobenius_norm is set to true.
%
%   .check_stability            [logical | {true}]
%       Compute the right/outermost eigenvalue of system matrix A in order
%       to assert that the system is stable.  Note that if no initial
%       perturbation is provided by the user, this computation will happen
%       regardless, in order to initialize the algorithm.  However, if a
%       complete initial perturbation is provided by the user, that is
%       epsilon*U*V', then the user may optionally disable asserting A's
%       stability by setting this option to false, as this eigenvalue solve
%       is typically the most difficult and expensive.  Note that this
%       should only ever be done if the user knows a priori that A is
%       stable, as that is a requirement for getStabRadBound's convergent
%       HEC phase.
%
%   .sparse_mode                [integer in {0,1,2} | {2}]
%       Controls whether getStabRadBound will use dense or sparse 
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
%       By default, getStabRadBound will compute an initial unit-norm
%       perturbation matrix U*V' based off of the right and left
%       eigenvectors of the rightmost (continuous-time systems) or
%       outermost (discrete-time systems) eigenvalue of system matrix A.
%       However, the user may optionally choose the kth rightmost/outermost
%       eigenvalue and/or specify their own initial perturbation matrix.
%       When a perturbation is supplied, the code chooses the kth
%       eigenvalue of the perturbed system matrix, for the user-supplied
%       perturbation U*V'. For more details, see
%       initialPerturbationOptions.
%   
%   .upperbound_opts            [struct | {upperBoundOptions()}]
%       This is a struct of parameters that are specific to finding an
%       initial destabilizing perturbation, which is required as a starting
%       condition for the convergent HEC phase of getStabRadBound.  For
%       more details, see upperBoundOptions.
%   
%   .expansion_opts             [struct | {uvExpandOptions()}]
%       This is a struct of the expansion subroutine's parameters, which is
%       used in finding an initial destabilizing perturbation as well as a
%       subphase of the convergent portion of getStabRadBound's algorithm.
%       For more details, see uvExpandOptions.
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
%       solved using various direct methods.  By default, getStabRadBound
%       will select the most efficient method automatically.  However, one
%       can also force which one of the available methods is employed via
%       this set of options.  Furthermore, this set of options allows one
%       to set pcg's options, when it is being used to iteratively solve
%       this pair of systems.  For more details, see
%       phiSystemsSolverOptions.
%
%   .count_multiplies           [logical | {false}]
%       When this is set to true, getStabRadBound will additionally return
%       the total number of matrix multiplications for A,A',B,B',C,C',D,D'.
%       Note enabling this adds some overhead.
%      
%   .record_level               [integer in {0,1,2,3} | {0}]
%       Controls how much metadata is included getStabRadBound's second
%       output argument.  Increasing this value increases memory and
%       computation overhead.
%           0 - only basic metadata regarding the total incurred costs
%           1 - breaks down the metadata on per phase basis, for  
%               initialization, finding an upper bound, and the HEC phase
%           2 - adds a history of the accepted iterates 
%           3 - additionally includes the rejected iterates as well 
%               (from line searches, rejected extrapolations, etc).
%       Note that this parameter has no effect if getStabRadBound is called
%       with only one output argument.
%
%   .record_UV                  [logical | {false}]
%       If this is true, the pair of vectors/matrices for each perturbation
%       U*V' will be additionally included in the recorded history that is
%       returned in getStabRadBound's second output argument.  Setting this
%       value to true can significantly increase memory consumption.  Note
%       that this parameter has no effect when getStabRadBound is called
%       with only one output argument or if opts.record_level = 0.
%       
%   .print_level                [integer in {0,1,2,3,4} | {1}]
%       Level of detail printed to console regarding progress:
%           0 - no printing whatsoever
%           1 - main printing level, overview of each phase:
%               initialization, finding an upper bound, and then the
%               convergent phase using HEC
%           2 - additionally prints out the individual iteration histories
%               of the subphases of finding an upper bound and HEC:
%                   - finding an initial upper bound
%                   - additional optimization (expansion) to improve the
%                     the initial maximizing frequency guess
%                   - HEC contractions 
%                   - HEC expansions
%           3 - additionally prints rejected steps in the above four
%               subphase histories, which may be either extrapolations,
%               interpolations, or if the line search fails (the last of
%               which can only happens on a subphase's very last step)
%           4 - additionally prints all line search evaluations in the
%               above four subphases.
%
%   .print_ascii                [logical | {false}]
%       By default, getStabRadBound's printed output uses the extended
%       character map, so nice looking tables can be made.  However,
%       diary() does not capture these symbols.  So, if you need to record
%       the output, you can restrict the printed output to only use the
%       basic ASCII character map, which may look better when captured by
%       diary().
%       Note: on Windows, opts.print_ascii is forced to true as the
%       extended ASCII used by MATLAB on WIndows are not monospaced.
%
%   See also getStabRadBound, initialPerturbationOptions,
%   upperBoundOptions, uvExpandOptions, epsilonContractOptions,
%   eigenSolverOptions, and phiSystemsSolverOptions.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   getStabRadBoundOptions.m introduced in ROSTAPACK Version 1.0
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
        validator = optionValidator('ROSTAPACK',defaults,sub_validators);
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
   
    validator.setUserOpts(user_opts);
    
    try 
        validator.setIntegerPositive('maxit');
        validator.setLogical('discrete_time');
        validator.setLogical('real_frobenius_norm');
        validator.setLogical('complex_ode_iteration');
        validator.setLogical('check_stability')
        validator.setIntegerInRange('sparse_mode',0,2);
        validator.setLogical('left_eigenvector');
        
        validator.setStructWithValidation('initial_perturbation');
        validator.setStructWithValidation('upperbound_opts');
        validator.setStructWithValidation('expansion_opts');
        validator.setStructWithValidation('contraction_opts');
        validator.setStructWithValidation('eig_solver_opts');
        validator.setStructWithValidation('phi_solver_opts');
                
        validator.setLogical('count_multiplies');
        validator.setIntegerInRange('record_level',0,3);
        validator.setLogical('record_UV');
        validator.setIntegerInRange('print_level',0,4);
        validator.setLogical('print_ascii');

        opts = validator.getValidatedOpts();
        
        % Extended ASCII chars in MATLAB on Windows are not monospaced so
        % don't support them.
        if ~opts.print_ascii
            validator.assert(~ispc(),                                   ...
                'only opts.print_ascii == true is supported on Windows.');
        end
        
    catch err
        err.throwAsCaller();
    end
end

function [default_opts,sub_validators] = getDefaults()
    default_opts = struct(                                              ...
        'maxit',                        100,                            ...
        'discrete_time',                false,                          ...
        'real_frobenius_norm',          false,                          ...
        'complex_ode_iteration',        false,                          ...
        'check_stability',              true,                           ...
        'sparse_mode',                  2,                              ...
        'left_eigenvector',             false,                          ...
        'initial_perturbation',         initialPerturbationOptions(),   ...
        'upperbound_opts',              upperBoundOptions(),            ...
        'expansion_opts',               uvExpandOptions(),              ...
        'contraction_opts',             epsilonContractOptions(),       ...
        'eig_solver_opts',              eigenSolverOptions(),           ...
        'phi_solver_opts',              phiSystemsSolverOptions(),      ...
        'count_multiplies',             false,                          ...
        'record_level',                 0,                              ...
        'record_UV',                    false,                          ...
        'print_level',                  1,                              ...
        'print_ascii',                  false || ispc()                 );
    
    sub_validators = struct(                                            ...
        'initial_perturbation', @initialPerturbationOptions,            ...
        'upperbound_opts',      @upperBoundOptions,                     ...
        'expansion_opts',       @uvExpandOptions,                       ...
        'contraction_opts',     @epsilonContractOptions,                ...
        'eig_solver_opts',      @eigenSolverOptions,                    ...
        'phi_solver_opts',      @phiSystemsSolverOptions                );  
end
