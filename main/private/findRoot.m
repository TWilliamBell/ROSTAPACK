function [x,f,delta,n_ub,n_root,n_full] = findRoot( bracketFn,rootFn,   ...
                                                    x0,f,df,step,       ...
                                                    use_interp,tol      )
%   findRoot:
%       Given a function f and initial point such that f(x0) > 0, findRoot
%       computes a root x of f such that x > x0 and f(x+d) < 0 for all
%       positive d sufficiently small (but only guaranteed in exact
%       arithemetic).  First, findRoot finds an upperbound ub_x on some
%       root such that f(ub_x) < 0 and then proceeds with one of three
%       root-finding methods, with bracketing to ensure convergence.
% 
%       The bracketing phase will take advantage of second-order
%       information (for Halley steps) if bracketFn returns both first and
%       second derivatives (i.e. as finite non-nan values), otherwise it
%       will just use first derivatives (for Newton steps).
%
%       The root finding phase will similarly use second-order information
%       (also Halley steps) if rootFn returns both first and second
%       derivatives and use_interp is set to false.  If rootFn only returns
%       first derivatives and use_interp is set to false, then the root
%       finding phase will use Newton steps.  If use_interp is set to true,
%       a first-order interpolation-based scheme will be used; in this
%       case, rootFn should not bother computing second derivatives since
%       they could be costly to compute and they will be ignored.
%    
%   INPUT:
%       bracketFn           [function of single real scalar input]
%       rootFn              [function of single real scalar input]
%           For a given real scalar x, these functions should have the
%           following signature:
%
%               [f,g,h] = fn(x),
%
%           where f(x), f'(x), and possibly f''(x) are returned as f, g,
%           and h, respectively.  If f''(x) is not computed, it must
%           nonetheless be returned, as nan.  bracketFn and rootFn are
%           provided separately, instead of as one function, to support
%           mixing and matching first- and second-order methods between the
%           bracketing and root finding phases.
%
%       x0                  [real scalar]
%           The initial starting point such that f(x0) > 0 holds.
% 
%       f                   [real scalar]
%           The value of f(x0) with f(x0) > 0
% 
%       df                  [real scalar]
%           The value of f'(x0)
% 
%       step                [real scalar]
%           Typically the value of the Newton or Halley step at x0, i.e. so
%           that x1 would be x0 + step.  Alternatively, this could be any
%           positive number as the bracketing method will consider
%           evaluating x0 + step as the first iteration in trying to find
%           an upper bound.
%
%       use_interp          [logical]
%           False means either Newton or Halley root finding with
%           bracketing will be used (depending on whether rootFn computes
%           second derivatives).  True means an interpolation scheme (with
%           automatically has bracketing) using only first derivatives will
%           be used.  In this case, the provided rootFn should not bother
%           computing second derivatives (as this would add expense
%           needlessly).
% 
%       tol                 [nonnegative real]
%           Desired relative accuracy in the computed root x.
%
%   OUTPUT:
%       x
%           The computed root
% 
%       f   
%           The computed value of f at this root, i.e. f(x)
%   
%       delta
%           The Newton/Halley step at root x, i.e. the next step if it were
%           to be taken. 
% 
%       n_ub
%           Number of steps incurred to find an upper bound on a root
% 
%       n_root
%           Number of steps incurred to find a root, not including the
%           bracketing phase. 
% 
%       n_full
%           Number of full steps in the root phase, i.e. all that were not
%           bisection steps.
%
% 
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   findRoot.m introduced in ROSTAPACK Version 2.0
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
  
    % Fixed parameter to not accept iterates that result in small but
    % nonzero abs(f).  If f is very flat, this could trigger termination
    % before high relative precision in root x is obtained root.  Instead,
    % we will only use that x is no longer changing more than than the
    % relative tolerance given by input tol.
    root_tol    = 0;    
    
    % Counters for the steps incurred for bracketing and root finding,
    % along with the number of full steps taken in the root finding phase.
    n_ub        = 0;  
    n_root      = 0;
    n_full      = 0;
       
    % Bracket root as quickly as possible by taking the larger of:
    % - twice the Newton/Halley step to the right (hence absolute value)
    % - difference of x - x0
    % Make sure initial difference is some positive number that will
    % increase x if it is added to it
    x                   = x0;
    delta               = step;
    diff                = max(1e-6,0.01*abs(x0));
    while f > 0
        lb_x            = x;
        lb_f            = f;
        lb_df           = df;
        if isinf(delta) % prevent delta_taken from ever being inf
            delta       = 0; 
        end
        delta_taken     = max(2*abs(delta),diff);
        x               = x + delta_taken;
        [f,df,delta]    = evaluateFn(bracketFn,x);                  
        diff            = x - x0;  
        n_ub            = n_ub + 1;
    end
    if f == 0
        return
    end
    ub_x            = x;
    ub_f            = f;
    ub_df           = df;
    reduction       = inf;
    altered         = true;
   
    if use_interp
        nextStepFn  = @interpBisection;
    else
        nextStepFn  = @dampedNewtonHalley;
    end
    
    % Find a root between the upper and lower bounds 
    while true

        x_old   = x;      
        x       = nextStepFn(x,delta);
        
        if x_old == x || x <= lb_x || x >= ub_x
            % Reached limits of the hardware
            x   = x_old;
            return
        end
     
        % Compute the Newton/Halley step as default for new iterate
        [f,df,delta]    = evaluateFn(rootFn,x);
        n_root          = n_root + 1;
        
        if abs(f) <= root_tol || abs(delta) <= tol*abs(x)
            % Converged: x is a root
            return
        elseif f < 0
            reduction   = f/ub_f;
            ub_x        = x;
            ub_f        = f;
            ub_df       = df;
        else
            reduction   = f/lb_f;
            lb_x        = x;
            lb_f        = f;
            lb_df       = df;
        end
        
        % If bracket is less than tolerance, no use in continuing
        if ub_x - lb_x <= abs(lb_x)*tol
            return
        end
    end
   
    function [f,df,delta] = evaluateFn(fn,x)
        [f,df,ddf]  = fn(x);
        delta       = rootStep(f,df,ddf);
        if ~isfinite(f)
            f       = inf;
        end
        if ~isfinite(delta)
            delta   = inf;
        end
    end

    function x = dampedNewtonHalley(x,delta) 
        if reduction > 0.9 && ~altered
            % The Newton/Halley step did not reduce the function value at
            % one of the bounds by at least 10% so take a damped step.
            % Since x is one of the bounds, ensureInBracket must return a
            % damped step, unless convergence has been attained
            % numerically.
            [x,altered] = ensureInBracket(x);
        else
            % Otherwise, take the Newton/Halley step, damped if necessary
            [x,altered] = ensureStepIsValid(x,delta);
            if ~altered
                n_full  = n_full + 1;
            end
        end
    end

    function [x,altered] = ensureStepIsValid(x,delta)
        if (delta > 0 && x == ub_x) || (delta < 0 && x == lb_x)
            % If delta is in the wrong direction, just bisect.
            x           = 0.5*(lb_x + ub_x);
            altered     = true;
        else
            % Otherwise, take the Newton/Halley step, possibly damped to
            % ensure the step remains within the lower and upper bounds.
            [x,altered] = ensureInBracket(x+delta); 
        end
    end
    
    function [x,altered] = ensureInBracket(x)
        
        % If the step is within the bounds, return it unaltered. Otherwise,
        % the step is too long so instead take a damped interpolated step
        % that is weighted toward the bound it would violate.  If this is
        % not numerically possible, i.e. the bounds are too close so
        % convergence has been actually already been attained numerically,
        % then return one of the bounds to trigger termination.
        
        altered = true;
        sw      = 0.1*(ub_x - lb_x);
        if x <= lb_x
            x = lb_x + sw;
            if x <= lb_x
                x = min(lb_x*(1 + 1e-15),ub_x);
            end
        elseif x >= ub_x
            x = ub_x - sw;
            if x >= ub_x
                x = max(ub_x*(1 - 1e-15),lb_x);
            end
        else
            altered = false;
        end
    end

    function x = interpBisection(varargin)
        
        % If the reduction from previous interpolation was insufficient,
        % then the interpolation must not have been sufficiently good, so
        % then do a bisection step.  But, make sure an interpolation is
        % at least attempted after every sufficient reduction. 
        
        if reduction > 0.9 && ~altered
            x       = 0.5*(lb_x + ub_x);
            altered = true;
        else
            x       = quadInterpStep();
            altered = false;
        end
    end
        
    function x = quadInterpStep()

        % use the curvature information from the bound that is closer 
        if abs(ub_f) < (lb_f)   % upper bound is closer to root
            x1  = lb_x;
            f1  = lb_f;
            x2  = ub_x;
            f2  = ub_f;
            df2 = ub_df;
        else                    % lower bound is closer to root
            x1  = ub_x;
            f1  = ub_f;
            x2  = lb_x;
            f2  = lb_f;
            df2 = lb_df;
        end
        [r1,r2] = getRootsOfQuadInterp(x1,f1,x2,f2,df2);
        
        if inBracket(r1)
            x       = r1;
            n_full  = n_full + 1;
        elseif inBracket(r2)
            x       = r2;
            n_full  = n_full + 1;
        else
            x       = 0.5*(lb_x + ub_x);
        end
    end

    function tf = inBracket(x)
        tf = ~isempty(x) && lb_x < x && x < ub_x;
    end
end

function [r1,r2,cFn] = getRootsOfQuadInterp(x1,f1,x2,f2,df2)
    A = [   1   x1  x1^2
            1   x2  x2^2
            0   1   2*x2 ];
    b = [ f1; f2; df2 ];
    c = A\b;

    cFn     = @(x) c(3)*x^2 + c(2)*x + c(1);
    [r1,r2] = quadRoots(c(3),c(2),c(1));
end

function [r1,r2] = quadRoots(a,b,c)
    r1      = [];
    r2      = [];
    b2_4ac  = b^2 - 4*a*c;
    if b2_4ac < 0
        return
    end
    discr   = sqrt(b2_4ac);
    r1      = real((-b + discr)/(2*a));
    r2      = real((-b - discr)/(2*a));
end