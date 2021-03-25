function solver = phiSystemsSolver(D_obj,varargin)
%   phiSystemsSolver:
%       This routine creates a solver object for the solving the following
%       pair of linear systems:
%       
%           Phi_m = (I - epsilon^2 D' * D)
%           Phi_p = (I - epsilon^2 D * D'),
%
%       which arise when approximating the complex-valued spectral value
%       set abscissa and/or radius.  The solver object automatically
%       determines the most efficient method for solving the systems based
%       on the properties of matrix D.
%
%       The available methods are:
%           0 - D is zero, nothing to solve!
%           1 - Use back slash
%           2 - Use two Cholesky factorizations of Phi_m and Phi_p
%               See phiSystemsFactorizeBoth
%           3 - Use one Cholesky factorization: smaller of Phi_m and Phi_p
%               See phiSystemsFactorizeSmaller
%           4 - Exploit low-rank structure D when given as D = Dl*Dr'
%               See phiSystemsFactorizeLowRank
%           5 - Solve the systems iteratively using pcg
%               See phiSystemsPCG.
%       
%       By default, methods 0, 3, 4, 5, will be chosen automatically when D
%       is respectively 0, an explicit matrix, an outer product, or a
%       function handle.  However, one may optionally force a method to be
%       used, provided the supplied D has a representation that is
%       compatible with one's choice of method.
%
%       When iterative solves are used, this solver object also provides a
%       facility to recycle information from previous solves in order to
%       help warm start subsequent related solves, with the hope that this
%       recycled information will help reduce the number of pcg iterations
%       needed on the next solves.
%       
%   INPUT:
%       D_obj
%           A matrixObject representation of matrix D.
% 
%       opts
%           Optional struct of options.  For more details, see
%           phiSystemsSolverOptions.
%   
%   OUTPUT:
%       solver
%           An "object", a struct containing the following functions for
%           solving the pair of linear systems: 
%
%       [b,c] = solver.solve(By,Cx,y,x)
%           Given right hand sides By and Cx, this returns:
%               b = Phi_m^(-1) By
%               c = Phi_p^(-1) Cx.
%           Input arguments y and x are also required when iterative solves
%           are employed, which used to prepare the initial vectors when
%           recycling is enabled.
%
%       .updateEpsilon(new_epsilon)  
%           Update the value of epsilon appearing in the linear systems.
%
%       .updateInitialVectors()
%           Update pcg's pair of initial vectors (one for each linear
%           system) with information from the most recently completely pair
%           of solves.  If the previous and next linear systems are
%           similar, doing this recycling will hopefully reduce the number
%           of pcg iterations needed to solve the next pair of linear
%           systems.  Note if iterative solves aren't being used, then this
%           function has no effect.
% 
%       totals = solver.getTotals()
%           Returns a struct containing the total number of times the pair
%           of linear systems has been solved (totals.phi_pair_solves) and,
%           if iterative solves are being used, the sum of the all pcg
%           iterations incurred over all of these solves (respectively
%           divided into totals.phi_m_iters and totals.phi_p_iters for the
%           two linear systems).
%
%   See also phiSystemsBackSlash, phiSystemsFactorizeBoth,
%   phiSystemsFactorizeLowRank, phiSystemsFactorizeSmaller, phiSystemsPCG,
%   and phiSystemsSolverOptions.
%
%
%   For more details, see [GGO13, Section 3], [Mit14, Section 4.1] and 
%   [MO16, Section 6].  
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   phiSystemsSolver.m introduced in ROSTAPACK Version 1.0
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

    % autodetect default method from D's type and get/process parameters
    [method,is_outer_product]   = getDefaultMethod(D_obj);
    opts                        = phiSystemsSolverOptions(varargin{:});
    
    if opts.force_method > 0
        assert( method > 0,                                             ...
                'phiSystemsSolver: D == 0, cannot force any method!'    );
        method = opts.force_method;
    end
          
    total_solves                = 0;    
    updated_solve_fn            = @solveEpsilonDZero;    
    update_initial_vectors_fn   = @NOP;
    get_totals_fn               = @getTotals;
    
    switch method
        case 0 
            solve_fn            = @solveEpsilonDZero;
            update_fn           = @NOP;
        case 1 
            [solve_fn,update_fn]    = phiSystemsBackSlash(D_obj);
        case 2 
            [solve_fn,update_fn]    = phiSystemsFactorizeBoth(D_obj);
        case 3 
            [solve_fn,update_fn]    = phiSystemsFactorizeSmaller(D_obj);
        case 4
            assert( is_outer_product,                                   ...
                [   'phiSystemsSolver: cannot force low-rank method as '...
                    'matrix D was not supplied as an outer product!'    ]);
            [solve_fn,update_fn]    = phiSystemsFactorizeLowRank(D_obj);
        case 5 
            get_totals_fn           = @getTotalsIterative;  
            [   solve_fn,                   ...
                update_fn,                  ...
                update_initial_vectors_fn,  ...
                get_iter_totals_fn          ] = phiSystemsPCG(D_obj,opts);
    end
    
    % This could be infinity if norm(D) is zero (which is fine).
    epsilon_ub = epsilonUpperBound(D_obj);
    
    solver = struct(                                                ...
                'solve',                @solve,                     ...
                'updateEpsilon',        @updateSolver,              ...
                'updateInitialVectors', update_initial_vectors_fn,  ... 
                'getTotals',            get_totals_fn               );
     
    function varargout = solve(varargin)
        total_solves            = total_solves + 1;
        [varargout{1:nargout}]  = updated_solve_fn(varargin{:}); 
    end
            
    function updateSolver(epsilon)
        assert( epsilon >= 0 && epsilon < epsilon_ub,               ...
                'phiSystemsSolver: epsilon must be in [0,%g).',     ...
                epsilon_ub                                          );
        if epsilon == 0 
            updated_solve_fn = @solveEpsilonDZero;
        else
            update_fn(epsilon);
            updated_solve_fn = solve_fn;
        end
    end

    function data = getTotals()
        data    = struct(   'phi_pair_solves',  total_solves    );
    end

    function data = getTotalsIterative()
        [total_m_iters, total_p_iters] = get_iter_totals_fn();
        data    = struct(   'phi_pair_solves',  total_solves,   ...
                            'phi_m_iters',      total_m_iters,  ...
                            'phi_p_iters',      total_p_iters   );
    end
end

function [method,is_outer_product] = getDefaultMethod(D_matrixObject)   
    is_outer_product = false;
    if D_matrixObject.isZero()
        method = 0;
        return
    end
    sparse_type = D_matrixObject.isSparse();
    if sparse_type == 3
        % matrix is only given as a function handle -> must use pcg
        method = 5;
    elseif sparse_type == 2
        % matrix is given as an outer product -> use low rank formulas
        is_outer_product = true;
        method = 4;
    else
        % matrix is given explicitly -> use min(m,p) form
        method = 3;
    end
end

function [b,c] = solveEpsilonDZero(By,Cx,varargin)
    % can't do anonymous function since b and c might be different lengths
    b = By;
    c = Cx;
end
