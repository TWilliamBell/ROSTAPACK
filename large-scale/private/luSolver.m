function solver = luSolver(use_permutations)
%   luSolver:
%       An "object" for applying the inverse of a given matrix to the left
%       or right side using LU factorizations, computing using regular LU
%       or with permutations.
%
%   INPUT:
%       use_permutations            [logical | {true}]
%           A logical indicating whether to use permutations when computing
%           LU factorizations.
%
%   OUTPUT:         
%       solver
%           An luSolver "object" with the following methods:
%
%           solver.factor(A)
%               Computes an internal LU factorization of a matrix A.
%               Matrix A can be dense or sparse.  Matrix A will be
%               automatically converted to a sparse matrix internally as
%               necessary, i.e. when use_permutations is true. 
%
%           [L,U,p,q] = solver.get()
%               Returns the current LU factorization.  Vectors p and q are
%               only nonempty if use_permutations was given as true.  If
%               solver.factor(A) has not been called, all output arguments
%               are returned as empty matrices.
%
%           tf = solver.hasInfNan()
%               Returns true if infs/nans are present in the LU
%               decomposition.
%
%           diff = solver.getNormDiff(A)
%               Computes the relative difference of a user-supplied matrix
%               A with the internally computed L*U, using the Frobenius
%               norm.  Returns nan if solver.factor() has not yet been
%               called.  Generally, one would provide the same matrix as
%               that given to solver.factor(A), for verifying that A was
%               indeed invertible and that LU decomposition is valid.
%           
%           X = solver.applyLeft(B)
%           X = solver.applyLeftCtrans(B)
%           X = solver.applyRight(B)
%           X = solver.applyRightCtrans(B)
%               Respectively computes X = A\B, A'\B, B/A, and B/A' using
%               the current LU factorization.  NOTE: solver.factor(A) must
%               have been already called at least once before any of these
%               functions can be called.
% 
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   luSolver.m introduced in ROSTAPACK Version 1.0
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

    L       = [];
    U       = [];
    p       = [];
    q       = [];
    diff    = [];
    
    if nargin < 1 || use_permutations
        factor_fn       = @factorVector;
        left_fn         = @applyLeftVector;
        left_trans_fn   = @applyLeftVectorCtrans;
        right_fn        = @applyRightVector;
        right_trans_fn  = @applyRightVectorCtrans;
    else
        factor_fn       = @factor;
        left_fn         = @applyLeft;
        left_trans_fn   = @applyLeftCtrans;
        right_fn        = @applyRight;
        right_trans_fn  = @applyRightCtrans;
    end
    
    solver = struct(    'factor',               factor_fn,      ...
                        'get',                  @get,           ...           
                        'getNormDiff',          @getNormDiff,   ...
                        'hasInfNan',            @hasInfNan,     ...
                        'applyLeft',            left_fn,        ...
                        'applyLeftCtrans',      left_trans_fn,   ...
                        'applyRight',           right_fn,       ...
                        'applyRightCtrans',     right_trans_fn   );

    function factor(A)
        diff        = [];
        [L,U]       = lu(A);
    end

    function factorVector(A)
        diff        = [];
        if ~issparse(A)
            A       = sparse(A);
        end
        [L,U,p,q]   = lu(A,'vector');
    end

    function [Lo,Uo,po,qo] = get()
        Lo          = L;
        Uo          = U;
        po          = p;
        qo          = q;  
    end

    function d = getNormDiff(A)
        if isempty(L)
            d       = nan;
            return
        end
        if isempty(diff)
            diff    = norm(A - L*U,'fro')/norm(A,'fro');
        end
        d           = diff;
    end

    function tf = hasInfNan()
        tf           = ~isFiniteValued(L) || ~isFiniteValued(U);
    end

    % No permutation vectors
    function X = applyLeft(B)
        X           = U \ (L \ B);
    end

    function X = applyRight(B)
        X           = (B / U) / L;
    end

    function X = applyLeftCtrans(B)
        X           = L' \ (U' \ B);
    end

    function X = applyRightCtrans(B)
        X           = (B / L') / U';
    end

    % Permutation vectors p,q
    function X = applyLeftVector(B)
        X(q,:)      = U \ (L \ B(p,:));
    end

    function X = applyRightVector(B)
        X(:,p)      = (B(:,q) / U) / L;
    end
   
    function X = applyLeftVectorCtrans(B)
        X(p,:)      = L' \ (U' \ B(q,:));
    end

    function X = applyRightVectorCtrans(B)
        X(:,q)      = (B(:,p) / L') / U';
    end

end