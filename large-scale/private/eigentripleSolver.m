function solver = eigentripleSolver(matrix,varargin)
%   eigentripleSolver:
%       An "object" for computing the largest real part or largest modulus
%       eigenvalues of matrix plus a low-rank perturbation A + U*V', along
%       with their associated right and left eigenvectors, satisfying RP or
%       RP(-z) compatibility.
%
%       eigentripleSolver supports both dense and sparse eigenvalue solves
%       (respectively via eig and eigs/eigsPlus).  Whether dense or sparse
%       eigenvalue solves or used is determined automatically based on the
%       properties of how matrix A was defined (dense/sparse matrix,
%       low-rank outer product, given as a function handle).
% 
%       For sparse eigensolves, eigenSolver also provides a facility for
%       recycling information from the last computation, to help reduce the
%       number of eigs iterations for a subsequent related computation.      
%
%   USAGE:
%       solver = eigentripleSolver(A);
%       solver = eigentripleSolver(A,sort_by_modulus);
%       solver = eigentripleSolver(A,sort_by_modulus,opts);
%
%   INPUT:
%       A                   [required]
%           The matrix must be passed as a matrixObject or in one of the
%           formats that matrixObject accepts.  Dense eigensolves are only
%           used if A is a dense matrix or is a matrixObject and
%           A.isSparse() returns false.
%       
%       sort_by_modulus     [optional: logical | {false}]
%           By default, eigenSolver will sort eigenvalues by largest real
%           part first.  When this input argument is provided and true, 
%           eigenSolver will instead sort eigenvalues by largest modulus
%           first.
%          
%       opts                [optional: struct of parameters]
%           An optional struct of settable parameters relevant for sparse
%           eigensolves or []. To see available parameters and their
%           descriptions, type: 
%           >> help eigenSolverOptions
%
%   OUTPUT:
%       solver
%           An "object", a struct containing the following functions for
%           computing eigenvalues and eigenvectors of a matrix:
%
%       [z,x] = solver.computeEigenpairsOfA(k)
%           Perform the computation to calculate the kth eigenpair of
%           matrix A, that is eigenvalue z and right eigenvector x.  If k
%           is not provided, then all eigenpairs computed are returned,
%           which is either all of them (dense) or those that converged
%           (for the sparse eigenvalue computations).
%
%       [x,y,absyhx] = solver.computeEigentripleOfA(z,x)
%           Perform the computation to calculate the corresponding left
%           eigenvector y for eigenvalue z and eigenvector x of A, while
%           scaling x and y so that they are RP or RP(-z) compatibile with
%           absyhx=abs(y'*x).
%
%       [z,X] = solver.computeEigenpairsOfAUV(U,V,k)
%           Analog of computeEigepairsOfA(k) for matrix A + U*V'.  
%           Argument k is likewise optional.
%
%       [x,y,absyhx] = solver.computeEigentripleOfAUV(z,x)
%           Analog of computeEigentripleOfA(k) for matrix A + U*V'.
%           Argument k is likewise optional.
%
%       tf = solver.isSparseSolver()
%           Returns true if sparse eigensolves are being used.
%
%       [z,x] = solver.computeRightEigenpair()
%           Perform the computation to calculate all/some eigenvalues of
%           matrix A and their right eigenvectors and returns the
%           desired/"selected" eigenvalue z and its right eigenvector x.
%           For more details, see solver.setWhichEigenvalue.
%
%       [x,y,absyhx] = solver.computeLeftEigenpair(z,x)
%           Given an eigenvalue z and its right eigenvector, this function
%           computes all/some eigenvalues of matrix A along with their
%           corresponding left eigenvectors.  This function then returns
%           the right eigenvector x and left eigenvector y, both scaled so
%           that they are RP or RP(-z) compatible, where absyhx=abs(y'*x).
%           Note that this function chooses y by finding the closest match
%           in the conjugates of the eigenvalue of A', such that absyhx is
%           not zero (if possible).
%      
%       tf = solver.isSparseSolver()
%           Returns true if sparse eigensolves are being used.
%
%       .updateInitialVectors()
%           Updates the pair of initial vectors to be used for the next
%           iterative sparse eigensolve computations to obtain both right
%           and left eigenvectors, based on the strategy selected by
%           opts.v0_recycling_type.
% 
%       .resetInitialVectorsToRandom()
%           Sets the pair of initial vectors to be used for the next
%           iterative sparse eigensolve computation to obtain both right
%           and left eigenvectors to random vectors.
% 
%       .setNumberRequested(k)
%           Sets the sparse eigensolvers to each request k eigenvalues to
%           obtain k right and k left eigenvectors.
%
%       .updateEigsOptions(opts)
%           Update the set of options used for eigs.  Valid fields for the
%           struct opts are:  .issym, .tol, .maxit, .p, and .disp.
%
%       tf = solver.isSortedByModulus()
%           Returns true if eigenvalues are sorted largest modulus first.
%           Returns false if eigenvalues are sorted by largest real part
%           first.
%   
%       totals = solver.getTotals()
%           Returns a struct containing:
%           .eig_solves_left    Total number of left eigenvector solves
%           .eig_solves_right   Total number of right eigenvector solves
%           
%           If eigs iteration counts are also available (eigsPlus is
%           being used), then totals will also contain the following two
%           fields, each set to the sum of all eigs iterations incurred
%           over all of their respective left/right eigenvector solves:
%           .eigs_iters_left 
%           .eigs_iters_right
%        
%   See also eigenSolverOptions and eigenSolver.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   eigentripleSolver.m introduced in ROSTAPACK Version 1.0
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
   
    try 
        [modulus,opts] = processOptions(varargin{:});
    catch err
        switch err.identifier
            case 'eigenSolver:invalidUserOption'
                err.throwAsCaller();
            otherwise
                error(invalidInputsMsg());
        end
    end

    rp_compatible_fn        = ternOp(modulus,@rpzCompatible,@rpCompatible);
    
    try 
        right_solver        = eigenSolver(matrix,modulus,false,opts);
        left_solver         = eigenSolver(matrix,modulus,true,opts);
    catch err
        err.throwAsCaller();
    end
    dim                     = right_solver.getDimension();
    
    solver = struct(                                                    ...
        'eigenpairsOfA',                @computeEigenpairsOfA,          ...
        'eigentripleOfA',               @computeEigentripleOfA,         ...
        'eigenpairsOfAUV',              @computeEigenpairsOfAUV,        ...
        'eigentripleOfAUV',             @computeEigentripleOfAUV,       ...
        'isSparseSolver',               @right_solver.isSparseSolver,   ...
        'updateInitialVectors',         @updateInitialVectors,          ...
        'resetInitialVectorsToRandom',  @resetInitialVectors,           ...
        'setNumberRequested',           @setNumberRequested,            ...
        'updateEigsOpts',               @updateEigsOptions,             ...
        'isSortedByModulus',            @() modulus,                    ...
        'getTotals',                    @getTotals                      );
  
    function [z,x] = computeEigenpairsOfA(varargin)
        right_solver.computeEigenvaluesOfA();
        [z,x] = getRightEigenpairs(varargin{:});
    end

    function [x,y,absyhx] = computeEigentripleOfA(z,x)
        left_solver.computeEigenvaluesOfA();
        [x,y,absyhx] = getEigentriple(z,x);
    end

    function [z,x] = computeEigenpairsOfAUV(U,V,varargin)
        right_solver.computeEigenvaluesOfAUV(U,V);
        [z,x] = getRightEigenpairs(varargin{:});
    end

    function [x,y,absyhx] = computeEigentripleOfAUV(U,V,z,x)
        left_solver.computeEigenvaluesOfAUV(U,V);
        [x,y,absyhx] = getEigentriple(z,x);
    end

    function [z,x] = getRightEigenpairs(kth_eval)
        if nargin < 1
            [z,x] = right_solver.getEigenpairsAll();
        else
            k = length(kth_eval);
            if k == 1
                [z,x] = right_solver.getEigenpair(kth_eval);
                return
            end
            z = zeros(k,1);
            x = zeros(dim,k);
            for j = 1:k
                [z(j), z(:,j)] = right_solver.getEigenpair(kth_eval(j));
            end
        end
    end

    function [x,y,absyhx] = getEigentriple(z,x)
        [~,y,diff] = left_solver.getEigenpairBestMatch(z,x);
        if diff > 1e-8*max(abs(z),1)
            warning('eigentripleSolver: eigenvalues differ by %g!', diff);
        end
        [x,y,absyhx] = rp_compatible_fn(x,y,z);
    end
    
    function setNumberRequested(k)
        right_solver.setNumberRequested(k);
        left_solver.setNumberRequested(k);
    end
  
    function updateInitialVectors()
        right_solver.updateInitialVector();
        left_solver.updateInitialVector();
    end

    function resetInitialVectors()
        right_solver.resetInitialVectorToRandom();
        left_solver.resetInitialVectorToRandom();
    end

    function updateEigsOptions(opts)
        right_solver.updateEigOptions(opts);
        left_solver.updateEigOptions(opts);
    end

    function totals = getTotals()
        left    = left_solver.getTotals();
        right   = right_solver.getTotals();
        
        if isfield(left,'iters') 
            iters   = {     'eigs_iters_left',  left.iters,     ...
                            'eigs_iters_right', right.iters      };
        else
            iters   = {};
        end
        
        totals  = struct(   'eig_solves_left',  left.solves,    ...
                            'eig_solves_right', right.solves,   ...
                            iters{:}                            );
    end
end

function [radius,opts] = processOptions(varargin)
    
    if nargin > 0 
        radius = varargin{1};
        assert(islogical(radius));
    else
        radius = false;
    end
  
    if nargin > 1
        opts = varargin{2};
        if isfield(opts,'verified') && opts.verified
            return
        end
    else
        opts = [];
    end
    opts = eigenSolverOptions(opts);
    opts.verified = true;
end

function m = invalidInputsMsg()
m = [                                                                   ...
'eigentripleSolver''s format is: \n'                                    ...
'  solver = eigentripleSolver(A)\n'                                     ...
'  solver = eigentripleSolver(A,sort_by_modulus)\n'                     ...
'  solver = eigentripleSolver(A,sort_by_modulus,opts)\n'                ];
end
