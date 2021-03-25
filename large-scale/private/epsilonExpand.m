function [  epsilon,        ...
            epsilon_diff,   ...
            f,f_diff,       ...
            info            ] = epsilonExpand(sys_c,opts,halt_on_unstable)
%   epsilonExpand:
%       For a fixed given unit-norm perturbation U*V' such that the
%       perturbed system matrix M(epsilon*U*V') is still stable for the
%       initial value of epsilon, this routine attempts to increase epsilon
%       so that M(epsilon*U*V') will no longer be stable.  This routine is
%       a subroutine of findUpperBound, which aims to find an initial
%       destabilizing perturbation that is necessary in order to run
%       hybridExpansionContraction.
%
%       Given the root function f, which is negative when M(epsilon*U*V')
%       is stable and nonnegative otherwise, we can increase the root
%       function towards being positive by increasing epsilon, provided
%       that the derivative of the root function is positve w.r.t. to
%       epsilon at the current perturbation.  Thus, when f < 0 holds, the
%       method can take a multiple of the Newton step, to hopefully
%       overshoot the root and make the root function positive.  When the
%       attempt does not satisfy monotonicity, a line search is done so
%       that f is always increased on each step.
% 
%       Note that this method can also attempt to increase the root
%       function f, even if it is already positive.
%      
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, already initialized at some perturbation
%           epsilon*U*V'.
%
%       opts            [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine.  For the main algorithmic parameters, type:
% 
%           >> help epsilonExpandOptions
%           
%           Additional parameters are:
%
%           .record_level       [value in {0,1,2}]
%               Determines how much metadata is gathered in the info output
%               argument:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the accepted iterates 
%               2 - additionally includes the rejected line search points.
%
%           .record_UV          [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be saved in the info output argument.
%
%           .print_level        [value in {0,1,2}]
%               0 - no printing whatsoever
%               1 - prints info for each accepted step 
%               2 - additionally prints rejected line search iterates
% 
%           .print_ascii        [logical]
%               Fallback to standard ASCII character set for printing table
%               borders
%
%           .printer
%               If a printerExpand object is already configured, it can
%               be provided here, to prevent instantiating a new one.
%
%       halt_on_unstable        [optional, logical]
%           Provide this as true if you want the routine to halt at the
%           first destabilizing perturbation it encounters. 
% 
%   OUTPUT:
%       If epsilonExpand makes a valid expansion, that is epsilon is
%       greater than its original value and f > f_initial, then sys_c will
%       be modified.  This enables the subsequent phases to reuse the
%       current computations without having to pass the data explicitly or
%       recompute things like eigenvectors.
%
%       If epsilonExpand cannot make a valid expansion, sys_c is restored
%       back to the state of its last snapshot.
%           
%       epsilon
%           The possibly increased value of epsilon
% 
%       epsilon_diff
%           epsilon - epsilon_initial, which is guaranteed to be >= 0
% 
%       f
%           The possibly new value of the root function for epsilon
%
%       f_diff
%           f - f_initial, whcih is guaranteed to be >= 0
%
%       info
%           Struct of data containing metadata about the computation:
%
%       .halt_status
%           0:  Maximum number of iterations reached
%           1:  Converged: relative difference condition attained
%               opts.rel_diff_tol satisfied
%           2:  Relative step size condition attained, halting early
%               opts.rel_step_size_tol satisfied
%           3:  Halted at first destabilizing perturbation encountered
%           4:  Line search failed to produce a monotonic step
%           5:  Magnitude of line search's derivative sufficiently small
%
%           Note that cases 1, 4, and 5 are indications of convergence to a
%           locally optimally point or stationary point (for a fixed
%           perturbation U*V').
% 
%       .stats 
%           .iters              Total number of iterations
%           .steps_accepted     Total number of expansion steps accepted 
%
%       .cost 
%           Total costs in terms of eigenvalue solves (and the number of 
%           eigs iterations, if applicable)
%   
%       .U and .V           [only present if opts.record_UV == true]
%           The current perturbation U*V' 
%
%       .iterates           [only present if opts.record_level > 0]
%           A struct array with data for each iterate of epsilonExpand.
%           Each entry has the following fields describing each point:
% 
%           .type           Type of the step: 'epsilon_initial', 'epsilon' 
%                           or 'epsilon_ls', the last of which is for
%                           points evaluated in the line search
%
%           .accepted       True if point is a valid contraction
%           .f              Value of the root function 
%           .z              The corresponding eigenvalue 
%           .epsilon        Value of epsilon
%           .info           Struct of line search metadata about the point:
%               .t              Value of t (in [0,1]) for line search
%               .bisection      True if a bisection line search step was
%                               taken (in lieu of using a quadratic or
%                               cubic interpolation step). 
%               .df_epsilon_t0  Value of the derivative of the root
%                               function with respect to epsilon, at t = 0
%               .newton_type    True if a Newton-based step was taken (used
%                               when the root function is negative).  False
%                               if epsilon was just increased by a scalar
%                               (used when f >= 0).
%               .damped         True if the Newton-based step was damped,
%                               to prevent epsilon from exceeding its upper
%                               bound.
%           .cost           Struct of cost data to obtain this point
%
%           Note that:
%           - .info.df_epsilon_t0, .info.newton_type, and .info.damped 
%             fields are only present on points that were accepted
%           - the .info.accepted field is only present when
%             opts.record_level > 1.  
%
%   CONSOLE:
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
%           N       multiple of Newton step to increase epsilon (r(z) < 0)
%           S       scale epsilon to increase it (r(z) >= 0)
%           ND      damped version of N, due to N step exceeding 1/norm(D)
%           SD      damped version of S, due to S step exceeding 1/norm(D)
%           LS      rejected line search evaluation 
%       Line Search
%           #       N,S,ND,SD: total number of line search evaluations 
%                   LS: total number of line search evaluations so far
%           t       N,S,ND,SD: final step length in (0,1] 
%                   LS: current trial step length
%
%   See also epsilonExpandOptions, epsilonExpandRecord, and findUpperBound.
%
%
%   For more details on the algorithm implemented by this routine, see
%   [MO16, Section 4.4].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   epsilonExpand.m introduced in ROSTAPACK Version 1.0
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
    
    maxit                   = opts.maxit;
    rel_diff_tol            = opts.rel_diff_tol;
    rel_step_size_tol       = opts.rel_step_size_tol;
    if nargin < 3 
        halt_on_unstable    = false;
    end
    
    step_multiplier         = opts.step_multiplier;
    limit_fraction          = opts.limit_fraction;
    
    ls_opts                 = opts.line_search_opts;
    ls_maxit                = ls_opts.maxit;
    ls_df_tol               = ls_opts.df_tol;
    ls_t_tol                = ls_opts.t_tol;
    ls_model_fn             = lineSearchInterpolation(ls_opts.model);
  
    [print_level,printer]   = printerConfig(opts,'expand');
    return_info             = nargout > 4;
    return_iterates         = return_info && opts.record_level > 0;
    if return_info
        record              = epsilonExpandRecord(sys_c,opts);
    end
    
    % if rel_step_size_tol is positive, quit when step size falls below 
    % some factor of largest step.  If rel_step_size_size is 0, max_f_diff
    % will be inf, and thus never updated to a finite value, which in turn
    % means that rel_step_size_halt will always be -inf and thus will never
    % trigger the relative step size termination condition
    rel_step_size_halt      = -inf;
    max_f_diff              = ternOp(rel_step_size_tol > 0, -inf, inf);
       
    f_initial               = sys_c.getf();
    epsilon_initial         = sys_c.getEpsilon();
    epsilon_upperbound      = sys_c.getEpsilonUpperBound();
     
    if print_level > 0 
        printer.msg('Expanding via increasing epsilon');
        printer.init(epsilon_initial,f_initial);
    end
    
    iter                = 0;
    steps_accepted      = 0;
    
    if (halt_on_unstable && f_initial > 0)
        if print_level
            printer.msg('Destabiling perturbation already attained.');
            printer.init(epsilon_initial,f_initial);
        end
        [epsilon,epsilon_diff,f,f_diff] = prepareTermination(3);
        return
    elseif print_level
        printer.msg('Expanding via epsilon updates');
        printer.init(epsilon_initial,f_initial);
    end
    
    halt_status = 0;
    for iter = 1:maxit
        
        f0 = sys_c.getf();          
        [ft,ft_diff,epsil_t,epsil_t_diff,t,ls_iters,type] = computeStep();
        
        if ft_diff > 0
            steps_accepted  = steps_accepted + 1;
            sys_c.updateEigsV0();
            sys_c.snapShot();
          
            if print_level > 0 
                printer.stepEF( iter, epsil_t, epsil_t_diff,            ...
                                ft, ft_diff, type, ls_iters, t          );
            end
          
            % Unlike uvExpand, don't bother trying an interpolation step.
            % The nonsmooth jitter observed in f as only epsilon varies
            % causes the derivative information to be highly localized.
            % Interpolatory quadratic models do not seem to be useful.
                
            if haltUpdating()
                break
            end
        else
            sys_c.restoreSnapShot();
            
            if print_level > 1 && t > 0
                printer.stepEFRejected( iter, epsil_t, epsil_t_diff,    ...
                                        ft, ft_diff, type, ls_iters, t  );
            end
            
            if ls_iters > 0     % failed to satisfy monotonicity
                halt_status = 4;
            else                % df_epsilon was too small
                halt_status = 5;
            end
            break
        end   
    end
    
    [epsilon,epsilon_diff,f,f_diff] = prepareTermination(halt_status);
    
    function [epsilon,epsilon_diff,f,f_diff] = prepareTermination(code)
        if print_level
            printer.msg(getTerminationMessage(code));
            printer.close();
        end 
        epsilon         = sys_c.getEpsilon();
        epsilon_diff    = epsilon - epsilon_initial;
        f               = sys_c.getf();
        f_diff          = f - f_initial;
        
        if return_info
            rec     = record.getRecord();
            info    = struct(   'halt_status',  code,               ...
                                'stats',        getStats(),         ...
                                rec{:}                              );
        end
    end
    
    function s = getStats()
        s = struct(     'iters',            iter,               ...
                        'steps_accepted',   steps_accepted      );
    end

    function halt = haltUpdating()
        % Compute step lengths (absolute and relative) in the change in f
        % values, in contrast to uvExpand where the step lengths 
        % are measured with respect to changes in the eigenvalues in the 
        % complex plane.  For uvExpand, the eigenvalues may oscillate and 
        % oscillation can indicate that a right/outer-most point has not 
        % yet been obtained.  Here, only epsilon is changed and thus
        % oscillation seems less likely to occur and furthermore, we 
        % neither expect nor intend to use this routine to converge to 
        % locally optimal points.  It thus suffices to have its termination 
        % criteria only be based upon how the f values are evolving.  
        f_new           = sys_c.getf();
        f_diff_abs      = abs(f_new - f0);
        if f0 == 0
            f_diff_rel  = inf;
        else
            f_diff_rel  = f_diff_abs / abs(f0);
        end
       
        % check whether any of the convergence critera are true
        halt = true;
        if f_diff_rel < rel_diff_tol
            halt_status = 1;
            return
        elseif f_diff_abs < rel_step_size_halt
            halt_status = 2;
            return
        elseif halt_on_unstable && f_new > 0
            halt_status = 3;
            return
        end
        halt = false;
        
        % update the longest step in f if necessary
        if f_diff_abs > max_f_diff
            max_f_diff          = f_diff_abs;
            rel_step_size_halt  = rel_step_size_tol*max_f_diff;
        end
    end
       
    function [  epsilon,    ...
                newton,     ...
                damped,     ...
                type        ] = computeNextEpsilon(f0,df_epsilon0,epsilon0)
        if f0 < 0   % take a multiple of the Newton step towards 0
            epsilon = epsilon0 - step_multiplier*(f0/df_epsilon0);
            newton  = true;
        else        % else just multiply the current epsilon
            epsilon = epsilon0*step_multiplier;
            newton  = false;
        end
        
        % make sure epsilon doesn't get too close to upperbound bound
        if epsilon > limit_fraction*epsilon_upperbound
            epsilon = 0.5*(epsilon0 + epsilon_upperbound);
            damped  = true;
        else
            damped  = false;
        end
        
        if print_level > 0
            type = getTypeStr(newton,damped);
        else
            type = [];
        end
    end

    function [ft,ft_diff,et,et_diff,t,ls_iters,type] = computeStep()
        
        epsilon0    = sys_c.getEpsilon();
        df_epsilon  = sys_c.getDfEpsilon();
        % Save the left eigenvector just computed for current point
        sys_c.snapShot();   
        
        [epsilon1,newton,damped,type] = computeNextEpsilon( f0,         ...
                                                            df_epsilon, ...
                                                            epsilon0    );
            
        [ft,ft_diff,et,et_diff,t,ls_iters] = deal(nan,0,epsilon1,0,0,0);
           
        if df_epsilon < ls_df_tol          
            if return_iterates
                 record.addEpsilon( false,ft,nan,et,t,                  ...
                                    df_epsilon,newton,damped            );
            end
            return
        end
        
        t           = 1;
        ls_iters    = 1; 
        [ft,z1]     = sys_c.computeEigEpsilon(epsilon1);
        ft_diff     = ft - f0;
        accepted    = ft_diff > 0;
        et_diff     = epsilon1 - epsilon0;
        
        if return_iterates
            record.addEpsilon(  accepted,ft,z1,epsilon1,t,              ...
                                df_epsilon,newton,damped                );
        end
    
        if accepted    
           return
        elseif print_level > 2
            printer.lineSearchEF(epsilon1,et_diff,ft,ft_diff,1,1);           
        end
        
        % otherwise, attempt a line search for monotonicity
        last_eval_empty = false;
        ls_eval_fn = getLineSearchEvaluator(f0,ft,                      ...
                                            epsilon0,epsilon1,df_epsilon);
        for ls_iters = 2:ls_maxit 
            [   ls_f_diff, ls_t, ls_f,      ...
                ls_epsilon, ls_epsilon_diff ] = ls_eval_fn(ls_iters);
            if isempty(ls_t)    % line search didn't evaluate 
                last_eval_empty = true;
                break
            end
            ft_diff = ls_f_diff;
            t       = ls_t;
            ft      = ls_f;
            et      = ls_epsilon;
            et_diff = ls_epsilon_diff;
            if ft_diff > 0       % line search point accepted
                break
            end
        end 
        if last_eval_empty
            ls_iters = ls_iters - 1;
        end
    end
    
    function eval_fn = getLineSearchEvaluator(  f0, f1,                 ...
                                                epsilon0, epsilon1,     ...
                                                df_epsilon0             )
        t_prev              = 1;
        ft_prev             = f1;
        get_t_from_model_fn = ls_model_fn(f0,f1,df_epsilon0);
    
        eval_fn = @evaluate;
        
        function [diff,t,ft,epsilon,epsilon_diff] = evaluate(iter)
            
            if isempty(t_prev)
                [diff,t,ft,epsilon] = deal(-inf,[],[],[]);
                return
            end
            
            [t,~,bisected] = get_t_from_model_fn(t_prev,ft_prev);
            if abs(t - t_prev) <= ls_t_tol
                [diff,t,ft,epsilon] = deal(-inf,[],[],[],[]);
                return
            end
            
            epsilon         = epsilon1*t + (1-t)*epsilon0;
            [ft,zt]         = sys_c.computeEigEpsilon(epsilon);
            diff            = ft - f0;
            epsilon_diff    = epsilon0 - epsilon;
            accepted        = diff > 0;
            
            if return_iterates
                record.addEpsilonLS(accepted,ft,zt,epsilon,t,bisected);
            end
            
            if accepted
                return;
            end
            
            if print_level > 2
                printer.lineSearchEF(epsilon,epsilon_diff,ft,diff,iter,t);
            end
            
            t_prev  = t;
            ft_prev = ft;
        end
    end
end

% For epsilonExpand
function m = getTerminationMessage(halt_status)
    switch halt_status
        case 0 
            m = 'Maximum number of iterations reached.';
        case 1
            m = 'Converged: relative difference condition attained.';
        case 2
            m = 'Relative step size condition attained.';
        case 3
            m = 'Halted at first destabilizing perturbation encountered.';
        case 4
            m = 'Line search failed to produce a monotonic step.';
        case 5
            m = 'Magnitude of line search''s derivative was sufficiently small.';
    end
    if halt_status > 3
        m = [m ' This is often a sign of a convergence.'];
    end
end

function type = getTypeStr(newton,damped)
    numerical_type = newton + 2*damped;
    switch numerical_type
        case 0
            type = 'S';
        case 1
            type = 'N';
        case 2
            type = 'SD';
        case 3
            type = 'ND';
    end
end   
