function [f,f_diff,info] = uvExpand(sys_c,opts,halt_on_unstable)
%   uvExpand:
%       For a fixed given perturbation level epsilon and an initial
%       unit-norm perturbation U*V', this routine produces a sequence of
%       new unit-norm perturbations that monotonically push an eigenvalue
%       of M(epsilon*U*V') rightward (continuous-time systems) or outward
%       (discrete-time systems), where
%
%           M(Delta) = A + B Delta (I - D Delta)^{-1} C
% 
%       is the perturbed system matrix.  This method almost always
%       converges to a locally rightmost/outermost point of the
%       epsilon-spectral value set, which in turn provides a locally-
%       optimal lower bound to the spectral value set abscissa|radius.  In
%       practice, this routine often converges to a point that is globally
%       rightmost/outermost, in which case, this method computes spectral
%       value set abscissa|radius.
%      
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, already initialized at some valid
%           perturbation epsilon*U*V'.
%
%       opts            [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine.  For the main algorithmic parameters, type:
% 
%           >> help uvExpandOptions
%           
%           Additional parameters are:
%
%           .record_level       [value in {0,1,2}]
%               Determines how much metadata is gathered in the info output
%               argument:
%               0 - only basic metadata regarding the total incurred costs
%               1 - adds a history of the accepted iterates 
%               2 - additionally includes the rejected points.
%
%           .record_UV          [logical]
%               Whether or not the sequence of perturbation
%               vectors/matrices U,V should also be saved in the info
%               output argument.
%
%           .print_level        [value in {0,1,2,3}]
%               0 - no printing whatsoever
%               1 - prints info for each accepted step 
%               2 - prints info for each rejected step
%               3 - additionally prints rejected line search evaluations.
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
%       If uvExpand makes a valid expansion, that is f is greater than its
%       original value (f_diff > 0), then sys_c will be modified and
%       updated to the new perturbation.  This enables the subsequent
%       phases to reuse the current computations without having to pass the
%       data explicitly or recompute things like eigenvectors.
%
%       If uvExpand cannot make a valid expansion, sys_c is restored back
%       to the state of its last snapshot.
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
%           6:  New perturbation U*V' was identical to previous one.
%
%       .converged
%           True if .halt_status was 1, 4, 5, or 6, as these are all
%           indications of convergence to a locally rightmost/outermost
%           point of the spectral value set (or at least a stationary
%           point).
%
%       .stats 
%           .iters              Total number of iterations
%           .steps_accepted     Total number of expansion steps accepted 
%           .interps            Total number of interpolations attempted
%           .interps_accepted   Total number of interpolations accepted
%           .extraps            Total number of extrapolations attempted
%           .extraps_accepted   Total number of extrapolations accepted
%
%       .epsilon
%           The value of epsilon used for the computation.
%
%       .cost 
%           Total costs in terms of the number of eigenvalue solves and, if
%           applicable, the number of linear system solves by
%           phiSystemsSolver.  The number of eigs and/or pcg iterations
%           will also be included if applicable.
%   
%       .iterates           [only present if opts.record_level > 0]
%           A struct array with data for each iterate of uvExpand.
%           Each entry has the following fields describing each point:
% 
%           .type           Type of the step:
%                           - UV_initial:
%                             The initial perturbation U*V'.
%                           - UV: 
%                             A full U*V' update step.
%                           - UV_ls:
%                             An evaluation inside the line search.
%                           - interp:
%                             An interpolation attempt.
%                           - extrap:
%                             An extrapolation attempt.
%           .accepted       True if point was accepted 
%           .f              Value of the root function 
%           .z              The corresponding eigenvalue 
%           .U and .V       The unit-norm perturbation U*V' [only present 
%                           if opts.record_UV == true]
%           .info           Struct of line search metadata about the point:              
%               .t              Value of t (in [0,1]) for line search
%               .bisection      True if a bisection line search step was
%                               taken (in lieu of using a quadratic or
%                               cubic interpolation step). 
%               .df_UV_t0       Value of the derivative of the root
%                               function with respect to t at t = 0 on the
%                               interpolated unit-norm path between the old
%                               perturbation and the full update step
%               .df_epsilon_t0  Value of the derivative of the root
%                               function with respect to epsilon at t = 0
%                               (that is, we take U(0) and V(0) as the 
%                               fixed perturbation, not the (un)accepted
%                               full update or line search evaluation)
%           .cost           Struct of cost data to obtain this point
%
%           Note that:
%           - .info.df_UV_t0 and .info.df_epsilon_t0 fields are only   
%             present on points that were accepted
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
%           UV      regular UV' update step
%           LS      rejected line search evaluation 
%           E       extrapolated step
%           I       interpolated step
%       Line Search
%           #       UV: total number of line search evaluations 
%                   LS: total number of line search evaluations so far
%                   (not applicable for E, I, or the initial value)
%           t       UV: final step length in (0,1] (typically accepted)
%                   LS: current trial step length
%                   I:  optimal length of step determined by interpolation
%                   (not applicable for E or the initial value)
%
%   See also uvExpandOptions and uvExpandRecord.
%
%
%   For more details on basic algorithms implemented by this routine, see
%   [GGO13] and [GGMO17] respectively.  For more details on the various
%   acceleration, efficiency, and reliablity enhancements that have also
%   been implemented, see [MO16, Section 6], [Mit14, Chapter 4 and Section
%   6.3.5], and [GGMO17, Section 7].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   uvExpand.m introduced in ROSTAPACK Version 1.0
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
    
    ls_opts                 = opts.line_search_opts;
    ls_maxit                = ls_opts.maxit;
    ls_df_tol               = ls_opts.df_tol;
    ls_t_tol                = ls_opts.t_tol;
    ls_model_fn             = lineSearchInterpolation(ls_opts.model);
    ls_UV_delta             = ls_opts.UV_delta;
    ls_dual_mode            = ls_opts.UV_dual_mode;
   
    if islogical(ls_dual_mode)
        if ls_dual_mode
            get_ls_eval_fn  = @getLineSearchEvaluatorDual;
        else
            get_ls_eval_fn  = @getLineSearchEvaluator;
        end
    else % ls_mode is a tolerance for choosing between the two modes
        get_ls_eval_fn      = @getLineSearchEvaluatorAdaptive;
    end
    
    interpolate_threshold   = opts.interp_threshold;
    interpolation_on        = interpolate_threshold < inf;
  
    extrap_size             = opts.extrap_size;
    extrapolation_on        = extrap_size > 1;
    extrap_min_size         = opts.extrap_min_size;
    extrap_rollover         = opts.extrap_rollover;
    extrap_initial_skip     = opts.extrap_initial_skip;
    
    [print_level,printer]   = printerConfig(opts,'expand');
    return_info             = nargout > 2;
    return_iterates         = return_info && opts.record_level > 0;
    if return_info
        record              = uvExpandRecord(sys_c,opts);
    end
   
    % if rel_step_size_tol is positive, quit when step size falls below 
    % some factor of largest step.  If rel_step_size_size is 0, max_f_diff
    % will be inf, and thus never updated to a finite value, which in turn
    % means that rel_step_size_halt will always be -inf and thus will never
    % trigger the relative step size termination condition
    rel_step_size_halt      = -inf;
    max_z_diff              = ternOp(rel_step_size_tol > 0, -inf, inf);
  
    epsilon                 = sys_c.getEpsilon();
    f_initial               = sys_c.getf();
          
    % counters
    iter                = 0;
    steps_accepted      = 0;
    interps             = 0;
    interps_accepted    = 0;
    extraps             = 0;
    extraps_accepted    = 0;
   
    if (halt_on_unstable && f_initial > 0) 
        if print_level
            printer.msg('Destabiling perturbation already attained.');
            printer.init(epsilon,f_initial);
        end
        [f,f_diff] = prepareTermination(3);
        return
    elseif print_level
        printer.msg('Expanding via UV updates');
        printer.init(epsilon,f_initial);
    end
    
    if interpolation_on
        quad_model_fn = lineSearchInterpolation(2);
    end
    if extrapolation_on
        extrapolator = extrapolateObject(   extrap_size,        ...
                                            extrap_min_size,    ...
                                            extrap_rollover     );
    end
    
    converged   = false;
    halt_status = 0;
    sys_c.updatePhiSolver();
    
    for iter = 1:maxit
        
        [f0,z0]                             = sys_c.getf();
        [U0,V0]                             = sys_c.getUV();
        [U1,V1,df_UV,df_eps]                = sys_c.computeNextUV();
        % Save the left eigenvector just computed for current point
        sys_c.snapShot();
        
        [f1,ft,ft_diff,Ut,Vt,t,ls_iters]    = computeStep();
     
        if ft_diff > 0
            steps_accepted  = steps_accepted + 1;
            sys_c.snapShot();
            
            % print accepted step info 
            if print_level > 0 
                printer.stepF(iter,ft,ft_diff,'UV',ls_iters,t);
            end
            
            % try an interpolated step
            if interpolation_on
                tryInterpolationStep();
            end
            
            % possibly try an extrapolated step
            if extrapolation_on && iter > extrap_initial_skip
                tryExtrapolationStep();
            end
            sys_c.updateEigsV0();
            sys_c.updatePhiV0(); 
                        
            if haltUpdating() 
                break
            end      
        else         
            sys_c.restoreSnapShot();
            
            if print_level > 1 && t > 0
                printer.stepFRejected(iter,ft,ft_diff,'UV',ls_iters,t);
            end
            
            if ls_iters > 0     % failed to satisfy monotonicity 
                halt_status = 4;
            elseif isnan(ft)    % df_UV was smaller than tolerance
                halt_status = 5;
            else                % U1 V1 was identical to U0 V0
                halt_status = 6;
            end
            converged = true;
            break;
        end
    end
    
    [f,f_diff] = prepareTermination(halt_status);
    
    function [f,f_diff] = prepareTermination(code)
        if print_level
            printer.msg(getTerminationMessage(code));
            printer.close();
        end
        f           = sys_c.getf();
        f_diff      = f - f_initial;
        if return_info
            rec     = record.getRecord();
            info    = struct(   'halt_status',      code,               ...
                                'converged',        converged,          ...
                                'stats',            getStats(),         ...
                                rec{:}                                  );
        end
    end

    function s = getStats()
        s    = struct(  'iters',            iter,               ...
                        'steps_accepted',   steps_accepted,     ...
                        'interps',          interps,            ...
                        'interps_accepted', interps_accepted,   ...
                        'extraps',          extraps,            ...
                        'extraps_accepted', extraps_accepted    );
    end
    
    function halt = haltUpdating()
        % Compute step lengths (absolute and relative) but, crucially,
        % of the eigenvalues in the complex plane, not the step lengths in
        % the f values.  Even though f always monotonically increases, the
        % corresponding eigenvalues can be oscillating (e.g. in the 
        % continuous-time case, the imaginary part may oscillate) and such
        % oscillation can be a sign that a locally right/outer-most point 
        % has not yet be obtained despite that f may be only increasing
        % extremely slowly.
        [f_new,z_new]   = sys_c.getf();
        if sys_c.isReal()
            % make sure we don't measure distance caused to conjugacy
            z_diff_abs  = abs(nonnegConj(z_new) - nonnegConj(z0));
        else
            z_diff_abs  = abs(z_new - z0);
        end
        if z0 == 0 || z_new == 0
            z_diff_rel  = inf;
        else
            z_diff_rel  = z_diff_abs / abs(z0);
        end
           
        % check whether any of the convergence critera are true
        halt = true;
        if z_diff_rel < rel_diff_tol 
            halt_status = 1;
            converged   = true;
            return
        elseif z_diff_abs < rel_step_size_halt 
            halt_status = 2;
            return
        elseif halt_on_unstable && f_new > 0
            halt_status = 3;
            return
        end
        halt = false;
        
        if z_diff_abs > max_z_diff
            max_z_diff          = z_diff_abs;
            rel_step_size_halt  = rel_step_size_tol*max_z_diff;
        end     
    end
    
    function [f_diff,accepted] = checkPoint(f_new,f_old)    
        f_diff      = f_new - f_old;
        accepted    = f_diff > 0;
        if accepted 
            sys_c.snapShot();
        else
            sys_c.restoreSnapShot();
        end
    end

    function [f1,ft,ft_diff,Ut,Vt,t,ls_iters] = computeStep()
    
        % set default values of output arguments
        [f1,ft,ft_diff,Ut,Vt,t,ls_iters] = deal(nan,nan,0,U1,V1,0,0);
       
        % check that the new perturbation isn't identical to old one or
        % if line search derivative is sufficiently small, don't bother
        U1V1_is_U0V0 = areIdentical(U0,U1) && areIdentical(V0,V1);
        if U1V1_is_U0V0 || df_UV < ls_df_tol
            if U1V1_is_U0V0
                [ft,zt] = sys_c.getf();
            else
                zt = nan;
            end
            if return_iterates
                record.addUV(false,ft,zt,Ut,Vt,t,df_UV,df_eps);
            end
            return
        end
       
        % try the full step
        t           = 1;
        ls_iters    = 1;
        [f1,z1]     = sys_c.computeEigUV(U1,V1);
        ft          = f1;
        ft_diff     = f1 - f0; 
        accepted    = ft_diff > 0;
        
        if return_iterates
            record.addUV(accepted,f1,z1,U1,V1,t,df_UV,df_eps);
        end
        
        if accepted  
            return
        elseif print_level > 2
            printer.lineSearchF(f1,ft_diff,1,1);
        end
         
        % otherwise, attempt the line search for monotonicity
        % start from 2 since full step counts as one iteration
        last_eval_empty = false;
        ls_eval_fn = get_ls_eval_fn(f0,f1,U0,V0,U1,V1,df_UV);
        for ls_iters = 2:ls_maxit % do a line search
            [ls_f_diff,ls_t,ls_f,ls_Ut,ls_Vt] = ls_eval_fn(ls_iters);
            if isempty(ls_t)    % line search didn't evaluate 
                last_eval_empty = true;
                break
            end
            ft_diff = ls_f_diff;
            t       = ls_t;
            ft      = ls_f;
            Ut      = ls_Ut;
            Vt      = ls_Vt;
            if ft_diff > 0      % line search point accepted
                break
            end
        end    
        if last_eval_empty
            ls_iters = ls_iters - 1;
        end
    end

    function tryInterpolationStep()
        interp_fn                   = quad_model_fn(f0,f1,df_UV);
        [t_in,f_in,bisected]        = interp_fn(t,ft);
        if ~bisected && (f_in - f0 > interpolate_threshold * ft_diff)
            interps                 = interps + 1;
            [u_in,v_in,t_in]        = uvLineInterpolate(U0,V0,U1,V1,t_in);
            [f_in,z_in]             = sys_c.computeEigUV(u_in,v_in);
            [f_in_diff,accepted]    = checkPoint(f_in,ft);
            
            if accepted 
                [Ut,Vt]             = deal(u_in,v_in);
                interps_accepted    = interps_accepted + 1;
                if print_level > 0
                    printer.interp(f_in,f_in_diff,t_in);                  
                end
            else
                if print_level > 1
                    printer.interpRejected(f_in,f_in_diff,t_in);
                end
            end
            
            if return_iterates
                record.addInterp(accepted,f_in,z_in,u_in,v_in,t_in);
            end
        end
    end

    function tryExtrapolationStep()
        % need to use getf() since we don't whether interpolation step
        % was taken after line search search step and we must getf() before
        % we compute a new eigenvalue
        f_current                   = sys_c.getf();
        try_early                   = iter == maxit;
        [U_ex,V_ex,k_used]          = extrapolator.update(Ut,Vt,try_early);
        if ~isempty(U_ex)
            extraps                 = extraps + 1;
            [f_ex,z_ex]             = sys_c.computeEigUV(U_ex,V_ex);
            [f_ex_diff,accepted]    = checkPoint(f_ex,f_current);
              
            if return_iterates
                record.addExtrap(accepted,f_ex,z_ex,U_ex,V_ex,k_used);
            end
              
            if accepted
                extraps_accepted    = extraps_accepted + 1;
                extrapolator.startNewSequence(U_ex,V_ex);     
                if print_level > 0
                    printer.extrap(f_ex,f_ex_diff);
                end
            else
                if print_level > 1
                    printer.extrapRejected(f_ex,f_ex_diff);
                end
            end
        end
    end

    function [Ut,Vt,t] = uvLineInterpolate(U0,V0,U1,V1,t)

        Ut = t*U1 + (1-t)*U0;
        Vt = t*V1 + (1-t)*V0;

        % checks if either Ut or Vt is identically zero
        % : are necessary in case Ut or Vt are matrices
        while ~any(Ut(:)) || ~any(Vt(:)) 
            % reduce t a little more
            t   = ls_UV_delta*t;  
            Ut  = t*U1 + (1-t)*U0;
            Vt  = t*V1 + (1-t)*V0;
            % must loop and check other vector is now not zero
        end

        [Ut,Vt] = sys_c.normalizeUV(Ut,Vt);
    end

    function eval_fn = getLineSearchEvaluator(f0,f1,U0,V0,U1,V1,df_UV,flip_UV)
        if nargin < 8
            flip_UV             = false;
        elseif flip_UV
            U1                  = -U1;
            V1                  = -V1;
        end
        t_prev                  = 1;
        ft_prev                 = f0;
        get_t_from_model_fn     = ls_model_fn(f0,f1,df_UV);
        uv_interpolate_fn       = @(t) uvLineInterpolate(U0,V0,U1,V1,t);
        
        eval_fn = @evaluate;
        
        function [diff,t,ft,Ut,Vt] = evaluate(iter)
           
            if isempty(t_prev)
                [diff,t,ft,Ut,Vt] = deal(-inf,[],[],[],[]);
                return
            end

            [t,~,bisected] = get_t_from_model_fn(t_prev,ft_prev);
            if abs(t-t_prev) <= ls_t_tol
                [diff,t,ft,Ut,Vt,t_prev] = deal(-inf,[],[],[],[],[]);
                return
            end

            [Ut,Vt,t]   = uv_interpolate_fn(t);
            [ft,zt]     = sys_c.computeEigUV(Ut,Vt);
            diff        = ft - f0;
            accepted    = diff > 0;
            
            if return_iterates
                record.addUVLS(accepted,ft,zt,Ut,Vt,t,bisected,flip_UV);
            end
            
            if accepted    
                return
            end
            
            % here we can print iterates that failed
            
            if print_level > 2
                if flip_UV 
                    t_print = -t;
                else
                    t_print = t;
                end
                printer.lineSearchF(ft,diff,iter,t_print);
            end
            
            t_prev  = t;
            ft_prev = ft;
        end
        
    end

    function eval_fn = getLineSearchEvaluatorDual(f0,f1,U0,V0,U1,V1,df_UV)
        
        eval_fn1    = getLineSearchEvaluator(f0,f1,U0,V0,U1,V1,df_UV,false);
        eval_fn2    = getLineSearchEvaluator(f0,f1,U0,V0,U1,V1,df_UV,true);     
        can_eval1   = true;
        can_eval2   = true;
     
        eval_fn     = @evaluate;
        
        function [diff,t,ft,Ut,Vt] = evaluate(iter)
            args_set = false;
            
            if can_eval1
                [diff,t,ft,Ut,Vt]       = eval_fn1(iter);
                args_set                = true;
                if isempty(t)
                    can_eval1           = false;
                elseif diff > 0
                    return
                end
            end
            
            if can_eval2
                [diff2,t2,ft2,Ut2,Vt2]  = eval_fn2(iter);
                if isempty(t2) 
                    can_eval2           = false;
                else
                    [diff,t,ft,Ut,Vt]   = deal(diff2,t2,ft2,Ut2,Vt2);
                    return
                end
            end   
            
            if ~args_set
                [diff,t,ft,Ut,Vt] = deal(-inf,[],[],[],[]);
            end
        end
    end

    function eval_fn = getLineSearchEvaluatorAdaptive(f0,f1,U0,V0,U1,V1,df_UV)
        if df_UV <= ls_dual_mode 
            eval_fn = getLineSearchEvaluatorDual(f0,f1,U0,V0,U1,V1,df_UV);
        else
            eval_fn = getLineSearchEvaluator(f0,f1,U0,V0,U1,V1,df_UV,false);
        end
    end
end

function z = nonnegConj(z)
    z = real(z) + 1i*sign(imag(z))*imag(z);
end

% For uvExpand
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
        case 6
            m = 'New perturbation was identical to current.';
    end
    if halt_status > 3
        m = [m ' This is often a sign of a convergence.'];
    end
end
