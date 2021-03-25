function [A_stable,info] = initializePerturbation(sys_c,opts)
%   initializePerturbation:
%       This routine calculates the starting eigenvalue to initialize a
%       spectral-value-set-based method at.  If a user supplies a complete
%       perturbation epsilon*U*V', where U*V' must have unit norm, the
%       initialization will be the kth eigenvalue (with the eigenvalues
%       with largest real part sorted first for continuous-time systems;
%       eigenvalues of largest modulus first for discrete-time systems) of
%       matrix M(epsilon*U*V'), where
%
%           M(Delta) = A + B Delta (I - D Delta)^{-1} C
% 
%       is the perturbed system matrix.  If the user does not specify any
%       perturbation information, or only epsilon or only U*V', this
%       routine will calculate the missing the values.  In this case, the
%       initialization will be done using the kth eigenvalue of system
%       matrix A and a perturbation epsilon*U*V' will either be entirely or
%       partially computed from this kth eigenvalue, as necessary.
% 
%       Regardless of what the user has provided, the final result of this
%       routine is that the sys_c is initialized at some perturbation
%       epsilon*U*V' where U*V' has unit norm.
% 
%       This routine will also check whether or not A is stable, unless the
%       user provides their own complete perturbation and specifically
%       disables this check.  In this case, the eigenvalues of A will not
%       be computed.  
%
%       Note that this method does not do any checks on the user-provided
%       inputs.  It is assumed that all are valid input values.
%
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, with no perturbation ever set.
%
%       opts            [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine.  
%
%           .check_stability            [logical]
%               When false and a complete perturbation epsilon*U*V' is
%               provided, the stability check on matrix A will be skipped.
%               If a complete perturbation is not provided, the stability
%               check on A will be done regardless of this parameter, since
%               the eigenvalues of A will be needed to create/complete the
%               initial perturbation.
%               
%           .initial_perturbation       [struct]
%               This struct is for specifying a partial/complete initial
%               perturbation and selecting which eigenvalue to initialize
%               at.  For more details, see initialPerturbationOptions.
%
%           .epsilon_step_multiplier    [value in [1,inf)]
%               When epsilon is not provided by the user, an initial guess
%               will be calculated using the eigenvalues of A, which is 
%               then multiplied by the value of this parameter. 
%
%           .epsilon_limit_fraction     [value in [0,1]] 
%               If D is nonzero, then epsilon < 1/norm(D) must always hold.
%               If necessary, the initial guess computed will be capped to
%               be at most the value of this parameter times 1/norm(D).  
%               This parameter has no effect if D is zero or if the user
%               supplies their own value of epsilon.
% 
%           .record_level       [value in {0,1}]
%               Determines how much metadata is gathered in the info output
%               argument:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the initialization iterates
%
%           .record_UV          [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be saved in the info output argument.
%
%           .print_level        [value in {0,1}]
%               0 - no printing whatsoever
%               1 - prints info about any or all of the stability of A, the
%                   initial eigenvalue selection, and the initial
%                   perturbation.
% 
%           .print_ascii        [logical]
%               Fallback to standard ASCII character set for printing table
%               borders
%
%           .printer
%               If a printerInit object is already configured, it can be
%               provided here, to prevent instantiating a new one.
% 
%   OUTPUT:
%       Note that that the state of sys_c will be modified by this routine,
%       regardless of its result.  It will be initialized at some
%       perturbation epsilon*U*V', where U*V' has unit norm.
%
%       A_stable
%           True if A if the computed rightmost/outermost eigenvalues are
%           within the stable region (which should indicate that A is
%           stable, unless eigs failed to converge to the desired
%           eigenvalues).  If the user provides their own complete
%           perturbation and opts.check_stability is false, the eigenvalues
%           of A will not be computed and this return value be set to true
%           regardless.  If this value is false, then it has been detected
%           that A is unstable.
%           
%       f0
%           Computed value of the stability root function matrix A.
%
%       z0 
%           The associated computed rightmost/outermost eigenvalue of A.  
%
%       info
%           Struct of data containing metadata about the computation:
%
%       .cost 
%           Total costs in terms of the number of eigenvalue solves and, if
%           applicable, the number of linear system solves by
%           phiSystemsSolver.  The number of eigs and/or pcg iterations
%           will also be included if applicable.
%   
%       .iterates           [only present if opts.record_level > 0]
%           A struct array with the following fields of data for each step
%           in the initialization:
% 
%           .type               Type of the step:
%                               - A stability:
%                                 Right/outermost eigenvalue
%                               - A eigenvalue:
%                                 kth right/outermost eigenvalue
%                               - Initial Perturbation:
%                                 1st or kth right/outermost eigenvalue
%           .f                  Value of the root function 
%           .z                  The corresponding eigenvalue 
%           .U and .V           Perturbation U*V'  [opts.record_UV == true]
%           .kth_eigenvalue     Which eigenvalue of A or M(Delta) this is   
%           .cost               Struct of cost data to obtain this point.
%
%   See also initializePerturbationRecord, initialPerturbationAssert, and
%   initialPerturbationOptions.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   initializePerturbation.m introduced in ROSTAPACK Version 1.0
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

    check_stability         = opts.check_stability;
    init_pert_info          = opts.initial_perturbation;
    
    kth_eigenvalue          = init_pert_info.kth_eigenvalue;
    epsilon                 = init_pert_info.epsilon;
    U                       = init_pert_info.U;
    V                       = init_pert_info.V;
    
    have_epsilon            = ~isempty(epsilon);
    have_UV                 = ~isempty(U) && ~isempty(V);
    user_supplied_Delta0    = have_epsilon && have_UV;
    
    [print_level,printer]   = printerConfig(opts,'init');
    return_info             = nargout > 1;
    return_iterates         = return_info && opts.record_level > 0;
    if return_info
        record              = initializePerturbationRecord(sys_c,opts);
    end
    
    epsilon_limit           = sys_c.getEpsilonUpperBound();
    A_stable                = true;
         
    if print_level 
        printer.msg(getRootFunctionMessage(sys_c.isDiscreteTime()));
    end
    
    f_init = nan;
    if user_supplied_Delta0        
        if check_stability 
            [A_stable,f_init] = computeEigA(1);
            if ~A_stable
                prepareExit('Matrix A is not stable - computation done.');
                return
            end
        end          
    else
        % If epsilon,u,v are not all provided, an initial perturbation
        % can be computed using the eigenvalues of A.  
        [A_stable,f_init,z_init,x_init] = computeEigA(kth_eigenvalue);
        if check_stability && ~A_stable
            prepareExit('Matrix A is not stable - computation done.');
            return
        end
        
        [eps_A,u_A,v_A] = sys_c.computeDeltaFromA(f_init,z_init,x_init);
        sys_c.snapShot();

        if ~have_epsilon
            % make sure epsilon_A isn't too close to 1/norm(D)
            epsilon = opts.epsilon_step_multiplier * eps_A;
            epsilon = min(  epsilon,                                    ...
                            opts.epsilon_limit_fraction * epsilon_limit );
        end
        if ~have_UV
            U = u_A;
            V = v_A;
        end
             
        % We already used the kth_eigenvalue from A to construct/complete
        % the initial perturbation so we should take the right/outermost
        % eigenvalue from now on
        kth_eigenvalue = 1;
    end
    [f,z] = sys_c.computeEigDelta(epsilon,U,V,kth_eigenvalue);  
    sys_c.snapShot();
   
    if print_level
        printer.initPert(epsilon,f,f-f_init,kth_eigenvalue);
    end
    
    if return_iterates
        record.addPerturbation(f,z,epsilon,U,V,kth_eigenvalue);
    end
    
    prepareExit('Initial perturbation computed.');
    
    function [is_stable,fk,zk,xk] = computeEigA(kth_eval)
        try 
            [f1,z1,fk,zk,xk] = sys_c.computeEigA(kth_eval);        
        catch err
            if nargout < 4
                msg     = stabilityCheckFailureMsg();
            else
                msg     = initialPerturbationFailureMsg();
            end
            ME = MException('getStabRadBound:initialPerturbation',msg);
            ME = addCause(ME,err);
            ME.throw();
        end  
        
        sys_c.snapShot();
        is_stable = f1 < 0;
        if print_level
            printer.matrixA(f1,1);
        end
        if return_iterates
            record.addAStability(f1,z1);
        end
        if kth_eval > 1
            if print_level
                printer.matrixA(fk,kth_eval);
            end
            if return_iterates
                record.addAEigenvalue(fk,zk,kth_eval);
            end
        end
    end

    function prepareExit(msg)
        if print_level
            printer.msg(msg);
            printer.close();
        end
        if return_info
            rec     = record.getRecord();
            info    = struct(rec{:});
        end
    end
end

function m = getRootFunctionMessage(discrete_time)
    if discrete_time
        m = 'Initializing - Root Function: Spectral Radius minus 1';
    else
        m = 'Initializing - Root Function: Spectral Abscissa';
    end
end

function m = initialPerturbationFailureMsg()
m = [                                                                   ...
'Failed to compute the rightmost/outermost eigenvalues of A in order \n'...
'obtain an initial perturbation.  If sparse mode is enabled, and the \n'...
'eigensolver failed, adjusting its options via \n\n'                    ...
'    opts.eig_solver_opts\n\n'                                          ...
'may allow the algorithm to proceed successfully.\n\n'                  ... 
'Alternatively, an initial perturbation can be provided explicitly \n'  ...
'by setting each of the following options appropriately: \n\n'          ...
'    opts.initial_perturbation.epsilon\n'                               ...
'    opts.initial_perturbation.U\n'                                     ...
'    opts.initial_perturbation.V\n\n'                                   ...
'Otherwise, make sure eigsPlus is installed and updated for your \n'    ...
'current version of MATLAB, via calling makeEigsPlus().\n\n'            ];
end


function m = stabilityCheckFailureMsg()
m = [                                                                   ...
'Stability check of matrix A computation failed.  If sparse mode is \n' ...
'enabled, and the eigensolver failed, adjusting its options via \n\n'   ...
'    opts.eig_solver_opts\n\n'                                          ...
'may allow the algorithm to successfully check the stability of A. \n\n'... 
'Alternatively, one can bypass the stability check (at your own risk)\n'...
'by specifying an initial perturbation via providing each of \n\n'      ...
'    opts.initial_perturbation.epsilon\n'                               ...
'    opts.initial_perturbation.U\n'                                     ...
'    opts.initial_perturbation.V\n\n'                                   ...
'and setting\n\n'                                                       ...
'    opts.check_stability \n\n'                                         ...
'to false.  Note that this should only be done if you know a priori \n' ...
'that matrix A is stable.'                                              ...
'Otherwise, make sure eigsPlus is installed and updated for your \n'    ...
'current version of MATLAB, via calling makeEigsPlus().\n\n'            ];
end
