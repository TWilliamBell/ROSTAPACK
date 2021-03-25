function [epsilon_new,f_new,restart,info] = epsilonContract(sys_c,opts)
%   epsilonContract:
%       The contraction phase of hybridExpansionContraction.  This routine
%       employs a hybrid Newton-bisection method, loosely based off of the
%       rtsafe method, to reduce the norm of destabilizing perturbation
%       epsilon*U*V' by reducing epsilon so that a rightmost/outermost
%       eigenvalue of the perturbed system matrix M(epsilon*U*V') is
%       "contracted" back to the stability boundary (or just to the right
%       of it by a small tolerance).  In other words, epsilonContract finds
%       a root of the function f which determines whether or not the matrix
%       M(epsilon*U*V') is stable, where f is negative when the matrix is
%       stable and nonnegative otherwise.
%
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, already initialized at some perturbation
%           epsilon*U*V' such that the perturbed system matrix
%           M(epsilon*U*V') is unstable.  That is the root function f is
%           positive.
%
%       opts            [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine.  For the main algorithmic parameters, type:
% 
%           >> help epsilonContractOptions
%           
%           Additional parameters are:
%
%           .restart            [struct of restart]
%               epsilonContract can be restarted from where a previous run
%               of epsilonContract terminated by settings opts.restart to
%               the restart output argument from the previous call.
% 
%           .record_level       [value in {0,1,2}]
%               Determines how much metadata is gathered in the info output
%               argument:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the accepted iterates 
%               2 - additionally includes the rejected iterates as well.  
%
%           .record_UV          [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be saved in the info output argument.
%
%           .print_level        [value in {0,1,2}]
%               0 - no printing whatsoever
%               1 - prints info for each accepted step 
%               2 - additionally prints rejected steps
% 
%           .print_ascii        [logical]
%               Fallback to standard ASCII character set for printing table
%               borders
%
%           .printer
%               If a printerContract object is already configured, it can
%               be provided here, to prevent instantiating a new one.
% 
%   OUTPUT:
%       If epsilonContract makes a valid reduction, that is, epsilon_new is
%       less than its original value and f_new is > 0, then sys_c will be 
%       modified.  This enables the subsequent expansion phase to reuse the
%       current computations without having to pass the data explicitly or
%       recompute things like eigenvectors.
%
%       If epsilonContract cannot make a valid reduction, sys_c is restored
%       back to the state of its last snapshot.
%           
%       epsilon_new
%           The possibly contracted value of epsilon
% 
%       f_new
%           The possibly new value of the root function for epsilon_new.
%           f_new is guaranteed to be in (0,f_initial).  If epsilonContract
%           has converged, f_new is guaranteed to be in (0,opts.tol).
%
%       restart             [note that no assertions are done]
%           Struct of data to restart a subsequent call to epsilonContract
%           from the last iteration, which is useful if the method did not
%           converge before halting.
%
%       info
%           Struct of data containing metadata about the computation:
%
%       .halt_status
%           0:  Maxit reached
%           1:  Converged: f_new is in (0,opts.tol)
%           2:  Sufficient progress made: 
%                   f is in (0,opts.rel_reduction_tol*f_initial)
%           3:  Root is bracketed by two consecutive floating point
%               numbers, cannot make further progress 
% 
%       .stats 
%           .iters              Total number of iterations
%           .newton_steps       Total number of Newton steps
%           .bisection_steps    Total number of bisection steps
%
%       .cost 
%           Total costs in terms of eigenvalue solves (and the number of
%           eigs iterations, if applicable)
%   
%       .U and .V           [only present if opts.record_UV == true]
%           The current perturbation U*V' 
%
%       .iterates           [only present if opts.record_level > 0]
%           A struct array with data for each iterate of the
%           Newton-bisection iteration.  Each entry has the following field
%           fields describing each point: 
% 
%           .type           Type of the step, always 'epsilon'
%           .accepted       True if point is a valid contraction
%           .f              Value of the root function 
%           .z              The corresponding eigenvalue 
%           .epsilon        Value of epsilon
%           .newton_step    True if a Newton step was taken, false is a
%                           bisection step was taken
%           .cost           Struct of cost data to obtain this point
% 
%       .accepted_index     [only present if opts.record_level > 0]
%           The index of the accepted point in .iterates
%
%       COMMON COLUMN HEADERS
%       Iter
%           The current iteration number of the relevant routine/phase
%       Epsilon
%           Value   Current value of epsilon in Delta = epsilon*U*V'
%           Diff    Absolute difference between current and previous values
%       Root Fn   
%           Value   Current value of r(z), where r(z) is either 
%                       r(z) = real(z)      (continuous time)
%                       r(z) = abs(z) - 1   (discrete time)
%                   and z is (typically) a globally rightmost/outermost
%                   eigenvalue of
%                       M(Delta) = A + B*Delta*(I - D*Delta)^{-1}*C
%                   for Delta = epsilon*U*V' with U*V' of unit norm.  Note
%                   that the minimal norm of Delta such that r(z) is zero 
%                   is the complex|real stability radius.
%           Diff    Absolute difference between current and previous values
%       
%       Note that:
%       1)  The value of epsilon (and the difference with its previous 
%           value) will only be printed when Iter is 0 or if epsilon has
%           been changed.  Otherwise, these two columns will remain blank.
%       2)  The values of epsilon and/or r(z) will be printed in ()'s for
%           rejected steps.  Rejected steps are only printed for
%           sufficiently high values of opts.print_level.
%
%       ROUTINE SPECIFIC COLUMN HEADERS
%       Step Type
%            -      initial value
%           N       Newton step 
%           B       bisection step
% 
%   See also epsilonContractOptions, epsilonContractRecord, and
%   hybridExpansionContraction.
% 
%
%   For more details on the algorithm implemented by this routine, see
%   [MO16, Section 4].
% 
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   epsilonContract.m introduced in ROSTAPACK Version 1.0 
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

    maxit                       = opts.maxit;
    tol                         = opts.tol;
    reduction_tol               = opts.rel_reduction_tol;
    
    [print_level,printer]       = printerConfig(opts,'contract');      
    return_info                 = nargout > 3;
    return_iterates             = return_info && opts.record_level > 0;
    if return_info
        record                  = epsilonContractRecord(sys_c,opts);
    end
    
    if isfield(opts,'restart')
        f_initial               = opts.restart.f;
        % sys_c must already be set to this perturbation U,V!
        U                       = opts.restart.U;
        V                       = opts.restart.V;
        epsilon_initial         = opts.restart.epsilon;
        df_epsilon              = opts.restart.df_epsilon;
        epsilon_lb              = opts.restart.epsilon_lb;
        epsilon_ub              = opts.restart.epsilon_ub;
        most_contracted         = opts.restart.most_contracted;
    else
        f_initial               = sys_c.getf();
        [U,V]                   = sys_c.getUV();
        epsilon_initial         = sys_c.getEpsilon();
        df_epsilon              = sys_c.getDfEpsilon();
        % Save the left eigenvector just computed for initial point
        sys_c.snapShot();
        epsilon_lb              = 0;
        epsilon_ub              = epsilon_initial;
        % keep track of iterate of the iterate with the smallest value of 
        % epsilon such that f is nonnegative.
        most_contracted         = epsilon_initial;
    end
    
    rel_reduction_target        = reduction_tol*f_initial;
    
    if print_level
        printer.msg('Contracting epsilon');
        printer.init(epsilon_initial,f_initial);
    end
    
    iter                = 0;
    number_bisect       = 0;
          
    % termination tolerance trick to get f_actual in [0,tol)
    f_actual            = f_initial;
    f                   = shiftF(f_actual);
   
    epsilon             = epsilon_initial;
    delta_old           = epsilon;
    accepted            = false;
    halt_status         = 3; % default code for exhausted precision
    update_accepted     = false;
    
    converged           = f_actual >= 0 && f_actual < tol;
    if converged 
        halt_status     = 1;
        accepted        = true;
        epsilon_new     = epsilon;   
    end
    
    while ~converged 
        % loop will exit when any of the following conditions hold
        % - tol condition is met
        % - rel_reduction_target is met
        % - reduction step is accepted and iter >= maxit
        % - floating point precision is exhausted 
        
        % Compute the Newton step as default for new iterate
        delta       = f/df_epsilon;
        epsilon_new = epsilon - delta;
      
        % Use bisection step instead if:
        % a) Newton step is out of bounds or
        % b) Newton step is not decreasing fast enough (rtsafe)
        if epsilon_new <= epsilon_lb  || epsilon_new >= epsilon_ub ...
                || abs(delta) > abs(0.5*delta_old)
            epsilon_new     = 0.5*(epsilon_lb + epsilon_ub);
            % if the computed mean of epsilon_lb and epsilon_ub fails to be
            % strictly between the two bounds, then this indicates that the
            % precision of the hardware has been exhausted and we must quit
            if epsilon_new <= epsilon_lb || epsilon_new >= epsilon_ub
                break
            end
            number_bisect   = number_bisect + 1;
            delta_old       = epsilon_ub - epsilon_lb;
            newton_step     = false;
        else
            delta_old       = delta;
            newton_step     = true;
        end
        
        [f_actual_new,z]    = sys_c.computeEigEpsilon(epsilon_new);
        f_diff              = f_actual_new - f_actual;
        f_actual            = f_actual_new;
        % only count function evaluations as a completed iteration
        iter                = iter + 1;
        
        % the best answer is the smallest value of epsilon such
        % that f >= 0.
        accepted = f_actual >= 0 && epsilon_new < most_contracted;
        if accepted
            update_accepted     = true;
            most_contracted     = epsilon_new;
            sys_c.updateEigsV0();
            sys_c.snapShot();
        end
        
        if return_iterates
            record.addEpsilon(accepted,f_actual,z,epsilon_new,newton_step);
        end
                
        if print_level
            if newton_step
                type = 'N';
            else
                type = 'B';
            end
            if accepted
                print_fn = @printer.stepEF;
            else
                print_fn = @printer.stepEFRejected;
            end     
            epsilon_diff = epsilon_new - epsilon;
            print_fn(iter,epsilon_new,epsilon_diff,f_actual,f_diff,type);
        end      
        
        % check for termination conditions 
        if f_actual >= 0 && f_actual < tol 
            % Absolute tolerance termination: f_actual in [0,tol)
            halt_status = 1;
            break
        elseif f_actual >= 0 && f_actual < rel_reduction_target
            % Relative reduction termination: requires f_actual >= 0
            halt_status = 2;
            break
        elseif most_contracted < epsilon_initial && iter >= maxit 
            % Some reduction achieved and at least maxit iters incurred
            halt_status = 0;
            break
        end
        
        % termination tolerance trick to get f_actual in [0,tol)
        f           = shiftF(f_actual);
        df_epsilon  = sys_c.getDfEpsilon();
        if accepted 
            % Save the left eigenvector just computed for the point, but
            % only for accepted points!  
            sys_c.snapShot();
        end
               
        if f < 0 
            epsilon_lb = epsilon_new;          
        else
            epsilon_ub = epsilon_new;
        end  
        epsilon = epsilon_new;         
    end
    
    if halt_status == 0 || halt_status == 2
        restart     = struct(   'f',                f_actual,           ...
                                'U',                U,                  ...
                                'V',                V,                  ...
                                'epsilon',          epsilon,            ...
                                'df_epsilon',       df_epsilon,         ...
                                'most_contracted',  most_contracted,    ...
                                'epsilon_lb',       epsilon_lb,         ...
                                'epsilon_ub',       epsilon_ub          );
    else
        restart     = [];
    end
 
    if ~accepted
        sys_c.restoreSnapShot();
        f_actual    = sys_c.getf();
        epsilon_new = most_contracted;
    end
    f_new           = f_actual; 
    
    if print_level
        printer.msg(getTerminationMessage(halt_status,update_accepted));
        printer.close();
    end
    
    if return_info                    
        rec     = record.getRecord();
        info    = struct(   'halt_status',  halt_status,            ...
                            'stats',        getStats(),             ...
                            rec{:}                                  );
    end
    
    function s = getStats()
        s = struct( 'iters',            iter,                       ...
                    'newton_steps',     iter - number_bisect,       ...
                    'bisection_steps',  number_bisect               );
    end

    function f = shiftF(f)
        f = f - 0.5*tol;
    end 
end

% For epsilonContract
function m = getTerminationMessage(code,contracted)
    switch code
        case 0
            m = ternOp(contracted,  'some contraction achieved.',       ...
                                    'no contraction achieved yet.'      );
            m = ['Maximum number of iterations reached - ' m];
        case 1 
           
            m = ternOp(contracted,  'attained.',                        ...
                                    'already attained at initial point.');
            m = ['Converged: absolute tolerance condition ' m];
        case 2
            m = 'Relative reduction condition attained - returning with partial contraction.';
        case 3
            m = ternOp(contracted, 'Some ', 'No ');
            m = { [m  'contraction achieved but root is bracketed by '  ...
                      'two consecutive floating point numbers.'],       ...
               'This may indicate that the tolerances are set too tight.'};
    end
end
