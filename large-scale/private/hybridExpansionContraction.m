function [halt_status,iters,info] = hybridExpansionContraction(sys_c,opts)
%   hybridExpandionContraction:
%       This routine implements the Hybrid Expansion-Contraction (HEC)
%       Algorithm to find a destabilizing perturbation of locally minimal
%       norm, whose norm is an upper bound to the complex|real stability
%       radius.
% 
%       Often this method converges to a perturbation whose norm is
%       globally minimal, in which case this value is the complex|real
%       stability radius.  When this does not happen, the resulting upper
%       bound is usually still a good approximation to the complex|real
%       stability radius.
%
%       Given a perturbation epsilon*U*V' such that U*V' has unit norm and
%       the perturbated system matrix M(epsilon*U*V') is unstable, where
% 
%           M(Delta) = A + B Delta (I - D Delta)^{-1} C,
%
%       this routine alternates between:
%           a) decreasing epsilon (epsilonContract) to bring a
%              rightmost/outermost eigenvalue back to the stability
%              boundary
%           b) and then updating the unit-norm perturbation U*V' 
%              (uvExpand) to expand rightward/outward again.
%
%       HEC is generally quadratically convergent, though this assumes that
%       both expansion and contraction phases are solved accurately.  In
%       practice, it can actually be faster to only solve these phases
%       inaccurately, even though it technically degrades HEC's convergence
%       rate to superlinear.  The reason in particular is that the
%       expansion phase typically only has linear convergence so reducing
%       its number of iterations is often a net gain, even if HEC requires
%       a few more "outer" iterations.
%       
%   INPUT:
%       sys_c           [required]
%           A system[Type] object, already initialized at some perturbation
%           epsilon*U*V' such that M(epsilon*U*V') is unstable and U*V' has
%           unit norm.
%
%       opts            [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine. 
%
%           .maxit              [positive integer]
%               Maximum allowed number of HEC iterations.
%
%           .expansion_opts     [struct of options for uvExpand]
%               This is the set of parameters for the uvExpand subroutine. 
%
%           .contraction_opts   [struct of options for epsilonContract]
%               This is the set of parameters for the epsilonContract
%               subroutine.
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
%                   epsilonContract and uvExpand subroutines
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
%       halt_status             
%           1:  HEC converged to tolerances.
%           2:  HEC stagnated as it can no longer make any contraction or
% `             expansion progress but the tolerances have not been
%               satisfied.  This is generally an indication that the
%               tolerances are too tight for the particular problem.
%           3:  Maximum number of iterations reached.
%
%       iters                  
%           Number of iterations incurred.  Each iteration consists of an
%           epsilonContract and uvExpand call.  As such, half iters are
%           possible.
%           
%       info
%           Struct of data containing metadata about the HEC computation:
%
%       .halt_status
%           Duplicate of first output argument halt_status.
% 
%       .iters
%           Duplicate of second output argument iters.
% 
%       .epsilon_phases
%           Number of contraction phases in HEC (i.e. the number of times
%           epsilonContract was invoked).
%   
%       .UV_phases
%           Number of expansion phases in HEC (i.e. the number of times
%           uvExpand was invoked).
%
%       .stats 
%           .epsilon            Accumulated stats from contraction phases
%           .UV                 Accumulated stats from expansion phases
%           These two substructs contain the sums of each statistic given 
%           in info.stats returned by epsilonContract and uvExpand,
%           respectively.  See their help docs for descriptions of the
%           statistics they each provide.
%
%       .cost 
%           .epsilon            Accumulated cost over all contractions
%           .UV                 Accumulated cost over all expansions
%           .total              Entire cost of hybridExpansionContraction
%           Costs are given in terms of the number of eigenvalue solves,
%           eigs iterations, etc.
%   
%       .phases                 [only present if opts.record_level > 1]
%           A cell array of the info output arguments returned by
%           epsilonContract and uvExpand, in the alternating order they
%           were invoked by hybridExpansionContraction.  As the routine
%           always commences with a contraction phase, contractions are
%           stored in the odd indices while expansions are stored in the
%           even indices.  Each info struct has the additional field 'type'
%           to label whether it was a contraction ('epsilon') or expansion 
%           ('UV') phase.
%
%   See also epsilonContract, epsilonContractOptions, recordIterates, and
%   uvExpand, uvExpandOptions.
%
%
%   For more details on the Hybrid Expansion-Contraction (HEC) Algorithm,
%   see [MO16].  Additional notes on properly implementing HEC in practice
%   are in [GGMO17, Appendix A].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   hybridExpansionContraction.m introduced in ROSTAPACK Version 1.0
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
    maxit                   = opts.maxit;
    record_level            = opts.record_level;
    expand_opts             = getOptions('expand');
    contract_opts           = getOptions('contract');
    hec_tol                 = 2*contract_opts.tol;
         
    if print_level
        printer = printers.main;
        printer.msg('Commencing Hybrid Expansion-Contraction (HEC)');
        printer.init(sys_c.getEpsilon(),sys_c.getf()); 
        [contract_print_fn,expand_print_fn,print_con_fn] = setupPrinters();
    end
    
    return_info     = record_level > 0 && nargout > 2;
       
    % These must not be anonymous function handles, since we need to be
    % able to ask the number of output arguments requested by the caller.
    contract_fn     = @epsilonContract;
    expand_fn       = @uvExpand;     
    
    if return_info
        record      = recordIterates(sys_c,contract_fn,expand_fn,opts);
        contract_fn = record.epsilon_fn;
        expand_fn   = record.uv_fn;
    end

    half_iter           = false;
    expanded            = true;
    stagnated           = false;
    converged           = false;
    expand_converged    = false;
    
    epsilon_opts        = contract_opts;
    
    for iter = 1:maxit
        
        epsilon_pre     = sys_c.getEpsilon();
        f_pre           = sys_c.getf();
              
        [epsilon,f,restart,con_info] = contract_fn(sys_c,epsilon_opts);
        
        % The difference values returned by contract_fn must always be
        % calculated manually since if epsilon_opts.restart is set,
        % contractEpsilon won't know the original values of epsilon and f.
        epsilon_diff    = epsilon - epsilon_pre;
        f_diff          = f - f_pre;
        contracted      = epsilon_diff < 0;
        contract_status = con_info.halt_status;
        
        if print_level
            contract_print_fn(con_info.stats.iters,contract_status);
        end
       
        % If epsilon has been reduced or if it hasn't but the previous
        % expansion phase did not converge (hit maxit or rel_step_size_tol
        % was used), then HEC should try another expansion phase.
        % Otherwise, epsilon is unchanged (meaning the precision has been
        % exhausted) and the previous expansion phase did converge.  This
        % means HEC cannot continue and it has either converged or
        % stagnated, depending on whether f is below hec_tol.
        if epsilon_diff == 0 && expand_converged
            half_iter       = true;
            if f < hec_tol
                converged   = true;
            else
                stagnated   = true;
            end
            if print_level > 1
                printHalfStepSummary();
            end
            break
        end
        
        [f,f_diff,exp_info] = expand_fn(sys_c,expand_opts);
        expanded            = f_diff > 0;
        expand_converged    = exp_info.converged;
       
        if print_level
            expand_print_fn(exp_info.stats.iters,exp_info.halt_status);
        end
        
        if expand_converged && f < hec_tol
            converged = true;
            break
        elseif f_diff > 0
            % otherwise, perturbation updated, so do a new contraction 
            epsilon_opts = contract_opts;
        elseif contract_status == 3
            % Contraction phase exhausted precision and f_diff is zero, 
            % which implies expand_converged is true and thus it must be 
            % that f >= hec_tol. Since no further progress can be made, HEC 
            % has stagnated.
            stagnated = true;
            break
        else
            % As f_diff is zero and f >= hec_tol holds, contract_status 
            % cannot be 1. Thus it is either 0 or 2 which indicates that 
            % the contraction phase made some reduction of epsilon but it
            % didn't converge to tolerance (1) or exhaust precision (3) so 
            % it may be able to make further contraction progress.
            % Therefore, the contraction phase should be restarted from
            % its last iterate and lower/upper bounds (which is
            % possible since the expansion phase did not update the
            % perturbation.
            epsilon_opts            = contract_opts;
            epsilon_opts.restart    = restart;
        end    
    end
    iters = ternOp(half_iter,iter-0.5,iter);
    
    if print_level
        printer.msg(getTerminationMessage(converged,stagnated));
        printer.close();
    end
   
    if converged
        halt_status = 1;
    elseif stagnated
        halt_status = 2;
    else
        halt_status = 3;
    end
    
    if return_info
        info = record.getRecord(halt_status,iters);
    end
    
    function [contract_print,expand_print,print_con_fn] = setupPrinters()
        print_con_fn        = [];
        if print_level > 1
            printer.close();
            contract_print  = @setContractSummaryPrintFn;
            expand_print    = @printFullStepSummary;
        else
            contract_print  = @contractPrint;
            expand_print    = @expandPrint;
        end
    end
    
    function contractPrint(siters,tc)
         if mod(iter,20) == 0
             printer.header();
         end
         printer.stepEF(    iter,       contracted,                 ...
                            epsilon,    epsilon_diff,   f,  f_diff, ...
                            'CON',      siters,         tc          );
    end
    
    function setContractSummaryPrintFn(siters,tc)
        print_con_fn = @() printer.stepEF(                          ...
                            iter,       contracted,                 ...
                            epsilon,    epsilon_diff,   f,  f_diff, ...
                            'CON',      siters,         tc          );
    end

    function expandPrint(siters,tc)
        printer.stepF([],expanded,f,f_diff,'EXP',siters,tc);
    end

    function printFullStepSummary(siters,tc)
        printer.msg('Single step of HEC summary');
        printer.header();
        print_con_fn();
        printer.stepF([],expanded,f,f_diff,'EXP',siters,tc);
        printer.close();
        fprintf('\n');
    end

    function printHalfStepSummary()
        printer.msg('Single (half) step of HEC summary');
        printer.header();
        print_con_fn();
        printer.close();
        fprintf('\n');
    end

    function sub_opts = getOptions(name)
        switch name
            case 'expand'
                sub_opts        = opts.expansion_opts;
            case 'contract'
                sub_opts        = opts.contraction_opts;
            otherwise 
                error('getOptions: field ''%s'' is not recognized!',name);
        end
        
        sub_opts.record_level   = max(record_level-1,0);
        sub_opts.record_UV      = opts.record_UV;
        sub_opts.print_level    = max(print_level-1,0);
        if sub_opts.print_level > 0
            sub_opts.printer    = printers.(name);
        end
    end
end

% For hybridExpansionContraction
function m = getTerminationMessage(converged,stagnated)
    if converged
        m = 'HEC: CONVERGED';
    elseif stagnated
        m = {                                                           ...
            'HEC: STAGNATED',                                           ...
            ' - Failed to obtain a destabilizing perturbation that is locally minimal to tolerances.',...
            ' - This may be an indication that the tolerances are too tight for this problem.'};
    else
        m = {                                                           ...
            'HEC: MAX ITERS REACHED',                                   ...
            ' - Rerunning with opts.maxit set larger may improve accuracy'};
    end
end
