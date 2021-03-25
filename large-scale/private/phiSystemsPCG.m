function [  solve_fn,           ...
            update_epsilon_fn,  ...
            update_initial_fn,  ...
            get_totals_fn       ] = phiSystemsPCG(D_obj,opts)
%   phiSystemsPCG:
%       This is a subroutine of phiSystemsSolver.
%
%       This routine solves the following pair of linear systems:
%       
%           Phi_m = (I - epsilon^2 D' * D)
%           Phi_p = (I - epsilon^2 D * D'),
%
%       which arise when approximating the complex-valued spectral value
%       set abscissa and/or radius.
%
%       This routine solves each system by iteratively solving both systems
%       using pcg (MATLAB's preconditioned conjugate gradients routine).
%       This method of solving this systems is only appropriate when one of
%       the other direct methods is either:
% 
%           a) not computationally feasible because m and p are both large,
%           where [p,m] = size(D), and D was not provided as a low-rank
%           outer product, or 
%
%           b) D was only given implicitly as a function handle.
%
%       Note that this routine has a facility to recycle information from
%       previous solves in order to warm start subsequent related solves,
%       with the hope that this recycled information will help reduce the
%       number of pcg iterations needed on the next solves.
%       
%   INPUT:
%       D_obj
%           A matrixObject representation.  If D is not an explicit full or
%           sparse matrix, a full representation will be created.
% 
%       opts
%           Optional struct of options.  For more details, see
%           phiSystemsSolverOptions.
%   
%   OUTPUT:
%       [b,c] = solve_fn(By,Cx,y,x)
%           Given right hand sides By and Cx, this returns:
%               b = Phi_m^(-1) By
%               c = Phi_p^(-1) Cx.
%           Input arguments y and x are also required and are used to
%           prepare the initial vector when recycling is enabled.
%
%       update_epsilon_fn(new_epsilon)  
%           Calling this will update the value of epsilon with the user's
%           new value.  Note that this method does NOT check if the new
%           value violates epsilon's upper bound: epsilon < 1/norm(D).
%
%       update_initial_fn()
%           Calling this will update pcg's pair of initial vectors (one for
%           each linear system) with information from the most recently
%           completely pair of solves.  If the previous and next linear
%           systems are similar, doing this recycling will hopefully reduce
%           the number of pcg iterations needed to solve the next pair of
%           linear systems.
% 
%       [phi_m_iters,phi_p_iters] = get_totals_fn()
%           Get the total number of pcg iterations incurred of all of the
%           linear systems solved so far.
%
%   See also phiSystemsBackSlash, phiSystemsFactorizeBoth,
%   phiSystemsFactorizeLowRank, phiSystemsFactorizeSmaller,
%   phiSystemsSolver, and phiSystemsSolverOptions.
%
%
%   For more details, see [GGO13, Section 3], [Mit14, Section 4.1] and 
%   [MO16, Section 6].  
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   phiSystemsPCG.m introduced in ROSTAPACK Version 1.0     
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
        
    % Default options for CG and initial vector
    maxit               = opts.maxit;
    tol                 = opts.tol;
    initial_vector_type = opts.initial_vector_type;
    
    [p,m]               = D_obj.getSize();
    applyD              = @D_obj.apply;
    applyDh             = @D_obj.applyHermitian;
        
    epsilon_sq          = 0; 
    
    % D as a functiom handle - use pcg()
    applyPhi_m          = @(v) v;
    applyPhi_p          = @(v) v;
    
    % Stored solutions b,c and y,x from By and Cx righthand sides
    % to use for warm starting initial vectors
    b_warm              = [];
    y_warm              = [];
    c_warm              = [];
    x_warm              = [];
    % Store most recent b,c and y,x so we can choose when to update 
    % the warm starting initial vectors (only when a u,v step is accepted)
    b_most_recent       = [];
    y_most_recent       = [];
    c_most_recent       = [];
    x_most_recent       = [];
    
    % iter info
    phi_m_iters         = 0;
    phi_p_iters         = 0; 
    
    phi_m_iters_total   = 0;
    phi_p_iters_total   = 0;
    
    switch initial_vector_type
        case 0 
            getCGInitialVector = @zeroVector;
        case 1
            getCGInitialVector = @randomVector;
        case 2
            getCGInitialVector = @passThrough;
        case 3
            getCGInitialVector = @warmStartCGVector;
        otherwise
            error('PCG warm start vector.');
    end    
    
    solve_fn            = @solveSystems;
    update_epsilon_fn   = @updateEpsilon;  
    update_initial_fn   = @updateInitialVectors;
    get_totals_fn       = @getTotals;
    
    function updateEpsilon(epsilon)
        new_epsilon_sq  = epsilon^2;
        if new_epsilon_sq ~= epsilon_sq
            epsilon_sq = new_epsilon_sq;        
            applyPhi_m = @(v) v - epsilon_sq*applyDh(applyD(v));
            applyPhi_p = @(v) v - epsilon_sq*applyD(applyDh(v));
        end
    end

    function [b,c] = solveSystems(By,Cx,y,x)        
        if isempty(b_warm)
            b0 = zeros(m,1);
            c0 = zeros(p,1);
        else
            b0 = getCGInitialVector(b_warm,y_warm,y);
            c0 = getCGInitialVector(c_warm,x_warm,x);
        end
              
        [b, ~, ~, phi_m_iters]  ...
            = pcg(applyPhi_m,By,tol,maxit,[],[],b0);
        
        [c, ~, ~, phi_p_iters]  ...
            = pcg(applyPhi_p,Cx,tol,maxit,[],[],c0);
        
        phi_m_iters_total = phi_m_iters_total + phi_m_iters;
        phi_p_iters_total = phi_p_iters_total + phi_p_iters;
                        
        b_most_recent = b;
        y_most_recent = y;
        
        c_most_recent = c;
        x_most_recent = x; 
    end

    function [phi_m_iters,phi_p_iters] = getTotals()
        phi_m_iters = phi_m_iters_total;
        phi_p_iters = phi_p_iters_total;
    end

    function updateInitialVectors()
        b_warm = b_most_recent;
        y_warm = y_most_recent;
        
        c_warm = c_most_recent;
        x_warm = x_most_recent;     
    end
end

function x0 = warmStartCGVector(soln_prev,x_prev,x)
    xdot    = x_prev'*x;
    theta   = real((0.5i)*log( (xdot^2) / (abs(xdot)^2) ));
    eitheta = exp(1i*theta);
    % Comment: exp(1i*(theta+pi)) == -exp(1i*theta)
    if norm(x - eitheta*x_prev) > norm(x + eitheta*x_prev)
        eitheta = -eitheta;
    end
    x0 = soln_prev * eitheta;
end

function soln_prev = passThrough(soln_prev,varargin)
    % Intentionally does nothing
end

function x0 = randomVector(soln_prev,varargin)
    x0 = rand(length(soln_prev),1);
end

function x0 = zeroVector(soln_prev,varargin)
    x0 = zeros(length(soln_prev),1);
end
