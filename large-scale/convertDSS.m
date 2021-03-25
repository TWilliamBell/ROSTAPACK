function [A_dss,B_dss] = convertDSS(E,A,B,use_permutations)
%   convertDSS:
%       This routine converts a descriptor system with nonsingular E:
% 
%           E*dx = A*x + B*u, 
%   
%       to a standard LTI:
%
%           dx = A_dss*x + B_dss*u,
%
%       where A_dss = E^{-1}*A and B_dss = E^{-1}*B.
%
%       Generally, E should be provided as a sparse matrix, so that a
%       (hopefully) sparse LU decomposition can be computed.  In this case,
%       A_dss and B_dss are not formed explicitly.  Instead, they are
%       returned as two cell arrays containing function handles and
%       metadata for respectively applying A_dss,A_dss' and B_dss,B_dss' to
%       vectors/matrices.  These two cell arrays can then be directly
%       passed to getStabRadBound as inputs for the A and B inputs.
%
%       Optionally, E can instead be provided as a full matrix, in which
%       case full matrix representations of matrices A and B will first be
%       formed (even if they were only supplied implicitly, e.g. as
%       function handles) and then, a dense LU factorization will be
%       computed in order to return matrices A_dss and B_dss, also
%       given explicitly as full matrices.
%
%       NOTE: when E is sparse, every matrix-vector product done with
%       either A_dss or B_dss will involve doing a backsolve with the LU
%       decomposition of E.  The additional cost of doing these backsolves
%       can be significant, particularly if the LU decomposition is not
%       sufficiently sparse.  
%
%   USAGE:
%       A_dss           = convertDSS(E);     
%       A_dss           = convertDSS(E,A);    
%       [A_dss,B_dss]   = convertDSS(E,A,B);
%       [A_dss,B_dss]   = convertDSS(E,A,B,use_permutations);
%
%   INPUT:
%       System matrix E         [required]
%           Nonsingular matrix E, dense or sparse, of the same dimension as
%           matrix A.  If the LU decomposition fails, i.e. it contains
%           infs/nans, an error will be thrown.  If the relative difference
%           between E and its LU matrix is too large, a message will
%           printed indicating that E may not actually be invertible and
%           thus the descriptor system cannot be transformed into a regular
%           LTI.
%
%       System matrix A         [optional | {[], which sets to identity}]
%           Matrix A can be specified in any of the following formats:
%           - as an explicit matrix, dense or sparse: A
%           - as an explicit matrix, dense or sparse, enclosed in single 
%             element cell array: {A}
%           - as an outer product A = U*V': {U,V}
%           - as a single function handle for multiplying A and A':
%                   {applyA,rows,cols,is_real}
%             where:
%                   applyA(x,false) returns A*x
%                   applyA(x,true) returns A'*x 
%                   rows is the number of rows of matrix A
%                   cols is the number of columns of matrix A
%                   is_real:logical, true if A only contains real entries
%           - as two separate function handles for multiplying A and A':
%                   {applyA,applyAh,rows,cols,is_real}
%             where:
%                   applyA(x,false) returns A*x
%                   applyAh(x,true) returns A'*x 
%                   rows is the number of rows of matrix A
%                   cols is the number of columns of matrix A
%                   is_real: logical, true if A only contains real entries
%
%           Note that the function handles should be able to multiply by A
%           and A' when x is either a vector or matrix.
% 
%       System matrix B          [optional | {[], which sets to identity}]
%           Matrix B can also be provided in any of above formats,
%           as well as via the following shortcut:
%
%       Logical use_permutations
%           By default, the LU decomposition will be computed using
%           permutations (i.e. lu(M,'vector') when M is sparse).  If you
%           wish to use the standard LAPACK routine, set use_permutations
%           to false.
%
%       Note: Although it is recommended to use the shortcut [] for
%       indicating that A or B should be the identity, you can of course
%       provide them explicitly as the identity.  However, in this case,
%       they should be sparse identities, i.e. via speye.  Otherwise,
%       computation and memory usage may be very inefficient if the
%       dimension of the system is large.  
%     
%   OUTPUT:
%       A_dss, B_dss
%           When E is a full matrix, A_dss and B_dss will be full matrices.
%
%           When E is a sparse matrix, A_dss will be a cell array of the
%           following form:
%
%               A_dss = {apply,applyHermitian,rows,cols,is_real}.
%           
%           where
%   
%               apply(x) 
%                   function handle to compute A_dss*x
%               applyHermitian(x)
%                   function handle to compute A_dss'*x
%               rows 
%                   the number of rows of A_dss
%               cols 
%                   the number of columns of A_dss
%               is_real
%                   logical, true if A_dss only has real-valued entries.
%
%           B_dss will also be returned in the same format.
%
%           Note that if A was provided as a function handle, then the
%           value of is_real for A_dss will be dependent on the caller's
%           declaration of whether or not A was real-valued (since it is
%           not practical to check this when the matrix is given
%           implicitly).  The same caveat applies for matrix B and B_dss.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   convertDSS.m introduced in ROSTAPACK Version 1.0   
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

    if nargin < 1
        error('convertDSS requires at least the first input, matrix E!');
    end
    
    assert(~isempty(E),'E must not be an empty matrix!');
    assert(isFiniteValued(E),'E must not contain infs or nans!');
    [n,E_cols]          = size(E);
    assertSquare(n,E_cols,'E');
    
    % If A or B is not provided, it is taken to be the identity.
    if nargin < 2 || isempty(A)
        A               = matrixObject({@(x) x, @(x) x, n, n, true});
    else
        A               = matrixObject(A);
        [A_rows,A_cols] = A.getSize();
        assertSquare(A_rows,A_cols,'A');
        assert(A_rows == n,'E and A must have the same dimensions!');
    end
    
    if nargin < 3 || isempty(B)
        B               = matrixObject({@(x) x, @(x) x, n, n, true});
    else
        B               = matrixObject(B);
        [B_rows,~]      = B.getSize();
        assert(B_rows == n,'E and B must have the same number of rows!');
    end
    
    if nargin < 4
        use_permutations     = true;
    else
        assert( isscalar(use_permutations) &&                       ...
                islogical(use_permutations),                        ...
                'use_permutations must be a single logical!'        );
    end
       
    lu_solver           = luSolver(use_permutations);
    lu_solver.factor(E);
    assert(     ~lu_solver.hasInfNan(),                             ...
                'LU decomposition of E contains infs/nans!'         );
    if lu_solver.getNormDiff(E) > 1e-8
        warning('convertDSS:singular',                              ...
                [   'the relative difference of E and its LU is '   ...
                    'greater than 1e-8!\nThis indicates that '      ...
                    'matrix E may not actually be invertible!'      ]);
    end
    applyEinv           = lu_solver.applyLeft;
    applyEinvHermitian  = lu_solver.applyLeftCtrans;
    build_fn            = ternOp(issparse(E),@buildSparse,@buildDense);

    A_dss               = build_fn(A);
    if nargout > 1  
        B_dss           = build_fn(B);
    end
    
    % Private functions 
 
    function M = buildDense(M)
        M               = M.formFull();
        M               = applyEinv(M);
    end
    
    function mat_args = buildSparse(M)
        [rows,cols]     = M.getSize();
        is_real         = M.isReal() && isRealValued(E);
        applyM          = @(x) applyEinv(M.apply(x));
        applyMh         = @(x) M.applyHermitian(applyEinvHermitian(x));
        mat_args        = {applyM,applyMh,rows,cols,is_real};
    end
    
end

function assertSquare(M_rows,M_cols,M_name)
    assert( M_rows == M_cols && M_rows > 0,                             ...
            'Matrix %s must be a nonempty square matrix!', M_name       );
end
