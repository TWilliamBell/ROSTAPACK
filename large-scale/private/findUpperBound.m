function [found_upperbound,iters,info] = findUpperBound(sys_c,opts)
%   findUpperBound:
%       This routine implements Algorithm Fast Upper Bound to find a
%       destabilizing perturbation epsilon*U*V', which in turn gives an
%       epsilon that is larger than the complex|real stability radius.
%       
%       The hybridExpansionContraction routine requires an initial
%       destabilizing perturbation.
%
%       Given a perturbation epsilon*U*V' such that the perturbated system
%       matrix M(epsilon*U*V') is stable, where
% 
%           M(Delta) = A + B Delta (I - D Delta)^{-1} C,
%
%       this routine alternates between:
%           a) increasing epsilon (epsilonExpand) to push a 
%              rightmost/outermost eigenvalue towards and hopefully past
%              the stability boundary
%           b) and then updating the unit norm perturbation U*V' 
%              (uvExpand) to expand rightward/outward.
%
%       If epsilon fixed, uvExpand may begin to take small steps as it
%       approaches a locally rightmost point in the stability region. Thus,
%       it is instead better to periodically increase epsilon for a small
%       number of uvExpand iterations.   
%
%       Once this routine has found a destabilizing perturbation has been
%       found, it can then optionally attempt to improve the quality of
%       destabilizing perturbation, by continuing to update the
%       perturbation matrix U*V' (via uvExpand, keeping epsilon fixed), in
%       order to better locate the maximizing frequency.
%
%   INPUT:
%       sys_c                   [required]
%           A system[Type] object, already initialized at some perturbation
%           epsilon*U*V'.  M(epsilon*U*V') be stable or unstable. 
%
%       opts                    [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine. 
%
%           .upperbound_opts    [struct of main options]
%               These are the main options governing the method. For more
%               details, see upperBoundOptions.
%
%           .expansion_opts     [struct of options for uvExpand]
%               This is the set of interpolation, extrapolation, and line
%               search options to be used for the uvExpand subroutine. Note
%               that the tolerances and .maxit are ignored, as their values
%               are governed via this method and the parameters set in
%               .upperbound_opts.  The only exception is that
%               .expansion_opts.maxit is used as iteration limit in the
%               last uvExpand phase to better locate the maximizing
%               frequency.  For more details, see uvExpandOptions.
% 
%           .record_level       [value in {0,1,2,3}]
%               Determines how much metadata is gathered in the info output
%               argument:
%               0 - no recording
%               1 - only basic metadata, total incurred costs, etc
%               2 - adds histories of the accepted iterates in each phase
%               3 - additionally includes the rejected iterates as well
%
%           .record_UV          [logical]
%               Whether or not the perturbation vectors/matrices U,V should
%               also be saved in the info output argument.
%
%           .print_level        [value in {0,1,2,3,4}]
%               0 - no printing whatsoever
%               1 - main printing level
%               2 - additionally prints out the iteration histories of the 
%                   epsilonExpand and uvExpand subroutines
%               3 - additionally prints rejected steps
%               4 - additionally prints rejected line search evaluations.
% 
%           .print_ascii        [logical]
%               Fallback to standard ASCII character set for printing table
%               borders
%
%           .printer
%               An optional struct of printers can be provided here, so 
%               new printers don't need to be instantiated if they already
%               exist for other routines.
% 
%   OUTPUT:
%       found_upperbound        [logical]
%           This will be true if the routine does finds an epsilon that is
%           larger than the complex|real stability radius.  
%
%       iters                   [nonnegative integer]
%           Number of iterations incurred.  Each iteration consists of an
%           epsilonExpand and uvExpand call.  As such, half iters are
%           possible.
%           
%       info
%           Struct of data containing metadata about the upper bound
%           computation:
%
%       .halt_status
%           0:  Maxit reached
%           1:  Upper bound found
% 
%       .iters
%           .upperbound         Number of steps to find an upper bound
%           .frequency          Number of steps to improve the frequency
% 
%       .epsilon_phases
%           Number of times the algorithm switched to increasing epsilon 
%           (i.e. the number of times epsilonExpand was invoked).
%   
%       .UV_phases
%           Number of times the algorithm switched to updating UV' (i.e.
%           the number of times uvExpand was invoked).  This includes both
%           the UV phases for finding an upper bound and improving the
%           frequency estimate once an upper bound has been found.
%
%       .stats 
%           .epsilon            Accumulated stats from epsilonExpand
%           .UV                 Accumulated stats from uvExpand
%           These two substructs contain the sums of each statistic given 
%           in info.stats returned by epsilonExpand and uvExpand,
%           respectively.  See their help docs for descriptions of the
%           statistics they each provide.
%
%       .cost 
%           .epsilon            Accumulated cost over all epsilon phases
%           .UV                 Accumulated cost over all UV phases
%           .total              Entire cost of findUpperBound
%           Costs are given in terms of the number of eigenvalue solves,
%           eigs iterations, etc.
%   
%       .phases                 [only present if opts.record_level > 1]
%           A cell array of the info output arguments returned by
%           epsilonExpand and uvExpand, in the order they were invoked by
%           findUpperBound.  Each info struct has the additional field
%           'type' to label whether it was an 'epsilon' or 'UV' phase.
%
%   See also epsilonExpand, epsilonExpandOptions, recordIterates, 
%   upperBoundOptions, uvExpand, and uvExpandOptions.
%
%
%   For more details, see [MO16, Section 4.4].
% 
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   findUpperBound.m introduced in ROSTAPACK Version 1.0
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

    [print_level,printers]  = printersConfig(opts);
    if print_level  
        printer             = printers.main;
    end
    
    ub_opts                 = opts.upperbound_opts;
    maxit                   = ub_opts.maxit;
    upperbound_quality      = ub_opts.upperbound_quality;
    update_UV_first         = ub_opts.update_UV_first;
    record_level            = opts.record_level;
    [UV_opts,UV_maxit]      = getOptionsUV();
    epsilon_opts            = getOptionsEpsilon();
                                                 
    f                       = sys_c.getf();
    f_initial               = f;
    epsilon                 = sys_c.getEpsilon();
    found_upperbound        = f > 0;
    
    iters                   = struct('upperbound',0,'frequency',0);
    info                    = [];
    return_info             = record_level > 0 && nargout > 2;
       
    % These must not be anonymous function handles, since we need to be
    % able to ask the number of output arguments requested by the caller.
    epsilon_fn              = @epsilonExpand;
    uv_fn                   = @uvExpand;      
    if return_info
        record              = recordIterates(sys_c,epsilon_fn,uv_fn,opts);
        epsilon_fn          = record.epsilon_fn;
        uv_fn               = record.uv_fn;
    end    
        
    % do alternating strategy to find an initial destabilizing perturbation
    if ~found_upperbound
        if print_level
            printer.msg('Commencing Upper Bound Search');
            printer.init(epsilon,f);
        end

        [update1_fn,update2_fn,print1_fn] = setupUpdateFunctions();

        half_iter       = false;
        second_failed   = false;
        for iter = 1:maxit    
            
            print_iter = iter;
            [f,first_failed]    = update1_fn();
            found_upperbound    = f > 0;
            stagnated           = first_failed && second_failed;
            if found_upperbound || stagnated
                if print_level > 1
                    printHalfStepSummary(print1_fn);
                end
                half_iter       = true;
                break
            end

            print_iter = [];
            [f,second_failed]   = update2_fn();
            found_upperbound    = f > 0;
            stagnated           = first_failed && second_failed ;
            if found_upperbound || stagnated  
                break
            end
        end
        iters.upperbound        = ternOp(half_iter,iter - 0.5,iter);
        epsilon = sys_c.getEpsilon();
        
        stagnated = first_failed && second_failed;
        if print_level
            printer.msg(getTerminatationMessage(found_upperbound,stagnated));
            printer.close();
        end     
    elseif upperbound_quality == 0
        if print_level 
            printer.msg('Upper bound already obtained at initial perturbation.');
            printer.close();
        end
    end
    
    if found_upperbound && upperbound_quality > 0      
        attemptToImproveFrequency();
    end
       
    if return_info
        info = record.getRecord(double(found_upperbound),iters);
    end
    
    function f_new = attemptToImproveFrequency()
        % get most recent f and epsilon since they may have been updated
        [f0,z0] = sys_c.getf();  
        if print_level
            if print_level > 1 && f_initial <= 0
                fprintf('\n\n\n');
            end
            printer.msg(getFrequencyMessage(upperbound_quality,f_initial));
            % print most recent values (since the alternating upper bound
            % search may have updated them)
            printer.init(epsilon,f0);
            if print_level > 1
                printer.close();
            end
        end
        UV_opts = getOptionsFrequency(UV_opts,UV_maxit,upperbound_quality);
        [f_new,f_diff,uv_info]  = uv_fn(sys_c,UV_opts);
        iters.frequency         = uv_info.stats.iters;
        accepted                = f_diff > 0;
        if accepted 
            f = f_new;
        end
        if print_level
            if print_level > 1 
                printer.msg('UPPER BOUND PROCEDURE: Summary of frequency improvement routine.');
            end     
            printer.stepF(  [],     accepted,                           ...
                            f_new,  f_diff,                             ...
                            'UV',   iters.frequency, uv_info.halt_status);
            if ~accepted 
                freq_del = 0;
            else
                [~,z1] = sys_c.getf();
                if sys_c.isDiscreteTime()
                    freq_del = abs(angle(z0) - angle(z1));
                else
                    freq_del = abs(imag(z0) - imag(z1));
                end
            end
            printer.msg(sprintf('Frequency delta: %.6g',freq_del));
            printer.close();
        end  
    end
    
    function [f,failed] = uvStep()
        [f,f_diff]  = uv_fn(sys_c,UV_opts,true);
        failed      = f_diff <= 0;
    end

    function [f,failed,print_fn] = uvStepWithPrintFn()   
        [f,f_diff,uv_info]  = uv_fn(sys_c,UV_opts,true);
        failed              = f_diff <= 0;
        print_fn = @() printer.stepF(                               ...
            print_iter,     ~failed,                                ...
            f,              f_diff,                                 ...
            'UV',           uv_info.stats.iters, uv_info.halt_status);
    end

    function [f,failed] = epsilonStep()   
        % We want epsilonExpand to halt as soon as it finds a destabilizing
        % perturbation, hence we pass true as the third argument 
        [~,eps_diff,f]  = epsilon_fn(sys_c,epsilon_opts,true);
        failed          = eps_diff <= 0;
    end

    function [f,failed,print_fn] = epsilonStepWithPrintFn()    
        % We want epsilonExpand to halt as soon as it finds a destabilizing
        % perturbation, hence we pass true as the third argument 
        [   eps_new, eps_diff,  ...
            f, f_diff, eps_info ] = epsilon_fn(sys_c,epsilon_opts,true);
        failed                           = eps_diff <= 0;
        print_fn = @() printer.stepEF(                                  ...
            print_iter,     ~failed,                                    ...
            eps_new,        eps_diff,                                   ...
            f,              f_diff,                                     ...
            'EPS',          eps_info.stats.iters, eps_info.halt_status  );
    end

    function [f,failed] = stepPart1WithPrinting(update_fn)
        [f,failed,print_fn] = update_fn();
        if mod(iter,20) == 0
            printer.header();
        end
        print_fn();
    end

    function [f,failed] = stepPart2WithPrinting(update_fn)
        [f,failed,print_fn] = update_fn();
        print_fn();
    end

    function [f,failed] = stepPart1WithSubprinting(update_fn)
         [f,failed,print1_fn] = update_fn();      
    end

    function [f,failed] = stepPart2WithSubprinting(update_fn)
        [f,failed,print_second_fn] = update_fn();
        printFullStepSummary(print1_fn,print_second_fn);
    end

    function printHalfStepSummary(print_first_fn)
        printer.msg('Single (half) step of upper bound search summary');
        printer.header();
        print_first_fn();
        printer.close();
        fprintf('\n');
    end

    function printFullStepSummary(print_first_fn,print_second_fn)
        printer.msg('Single step of upper bound search summary');       
        printer.header();
        print_first_fn();
        print_second_fn();
        printer.close();
        fprintf('\n');
    end

    function [update1_fn,update2_fn,print1_fn] = setupUpdateFunctions()
        print1_fn = [];
        if print_level < 1
            if update_UV_first
                update1_fn = @uvStep;
                update2_fn = @epsilonStep;
            else
                update1_fn = @epsilonStep;
                update2_fn = @uvStep;
            end
            return
        else
        
            if update_UV_first
                step1_fn = @uvStepWithPrintFn;
                step2_fn = @epsilonStepWithPrintFn;
            else
                step1_fn = @epsilonStepWithPrintFn;
                step2_fn = @uvStepWithPrintFn;
            end   
            if print_level > 1
                printer.close();
                update1_fn = @() stepPart1WithSubprinting(step1_fn);
                update2_fn = @() stepPart2WithSubprinting(step2_fn);
            else
                update1_fn = @() stepPart1WithPrinting(step1_fn);
                update2_fn = @() stepPart2WithPrinting(step2_fn);        
            end
        end
    end

    function [sub_opts,UV_maxit] = getOptionsUV()
        sub_opts                    = opts.expansion_opts;
        UV_maxit                    = sub_opts.maxit;
        sub_opts.maxit              = ub_opts.UV_steps_per_iter;
        sub_opts.rel_diff_tol       = ub_opts.rel_diff_tol;
        sub_opts.rel_step_size_tol  = ub_opts.rel_step_size_tol;
        sub_opts                    = addRecordAndPrintOptions(sub_opts);  
    end

    function opts = getOptionsEpsilon()
        opts = struct(                                                  ...
            'maxit',                ub_opts.epsilon_steps_per_iter,     ...
            'rel_diff_tol',         ub_opts.rel_diff_tol,               ...
            'rel_step_size_tol',    ub_opts.rel_step_size_tol,          ...
            'step_multiplier',      ub_opts.epsilon_step_multiplier,    ...
            'limit_fraction',       ub_opts.epsilon_limit_fraction,     ...
            'line_search_opts',     ub_opts.epsilon_line_search_opts    );
        opts = addRecordAndPrintOptions(opts);  
    end
    
    function sub_opts = addRecordAndPrintOptions(sub_opts)
        sub_opts.record_level       = max(record_level - 1,0);
        sub_opts.record_UV          = opts.record_UV;
        sub_opts.print_level        = max(print_level - 1,0);
        if sub_opts.print_level > 0 
            sub_opts.printer        = printers.expand;
        end
    end
end

function UV_opts = getOptionsFrequency(UV_opts,UV_maxit,upperbound_quality)
    UV_opts.maxit               = UV_maxit;
    UV_opts.rel_step_size_tol   = 1 - upperbound_quality;
end

function m = getFrequencyMessage(upperbound_quality,f0)
    m = sprintf('Attempting to improve frequency of upper bound with quality level: %.2e',upperbound_quality);
    if f0 > 0
        m = {'Upper bound already obtained at initial perturbation.' m};
    end
end

function m = getTerminatationMessage(found_upperbound,stagnated)
    if found_upperbound
        m = 'UPPER BOUND PROCEDURE: Minimal destabilizing epsilon bounded from above.';
    elseif stagnated
        m = {                                                           ...
            'UPPER BOUND PROCEDURE: STAGNATED',                         ...
            ' - Iteration has converged to a stabilizing perturbation.', ...
            ' - Adjust options and/or try a different starting perturbation.'};
    else
        m = {                                                           ...
            'UPPER BOUND PROCEDURE: MAX ITERS REACHED',                 ...
            ' - Increase opts.upperbound_opts.maxit and try again.'      };
    end
end
