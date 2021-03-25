function solver = eigenSolver(matrix,varargin)
%   eigenSolver:
%       An "object" for computing the largest real part or largest modulus
%       eigenvalues of matrix plus a low-rank perturbation A + U*V', along
%       with either their associated right or left eigenvectors.
%
%       eigenSolver supports both dense and sparse eigenvalue solves
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
%       solver = eigenSolver(A);
%       solver = eigenSolver(A,sort_by_modulus);
%       solver = eigenSolver(A,sort_by_modulus,left_eigenvector);
%       solver = eigenSolver(A,sort_by_modulus,left_eigenvector,opts);
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
%       left_eigenvector    [optional: logical | {false}]
%           By default, eigenSolver will compute right eigenvectors of
%           matrix A.  When this input argument is provided and true,
%           eigensolver will instead compute left eigenvectors of matrix A.
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
%       .computeEigenvaluesOfA()
%           Perform the computation to calculate all/some eigenvalues of
%           matrix A.
%
%       .computeEigenvaluesOfAUV(U,V)
%           Perform the computation to calculate all/some eigenvalues of
%           matrix A + U*V'.  Note that this matrix is never formed 
%           explicitly when sparse eigenvalue solves are employed.  
%
%       z = solver.getEigenvalue(kth_eigenvalue)
%           This "selects" the kth eigenvalue from the computed list of
%           eigenvalues sorted by largest real part or modulus and then
%           returns this eigenvalue to the caller.  This function should
%           only ever be called after one of the solver.computeEigenvalues
%           functions has been called at least once.
%
%       x = solver.getEigenvector()
%           Return the eigenvector for the currently "selected" eigenvalue.
%           This function should only ever be called after one of the
%           solver.computeEigenvalues functions has been called at least
%           once.
%
%       [z,x] = solver.getEigenpair(kth_eigenvalue)
%           This "selects" the kth eigenvalue from the computed list of
%           eigenvalues sorted by largest real part or modulus and then
%           returns this eigenvalue z and its corresponding eigenvector x
%           to the caller.  This function should only ever be called after
%           one of the solver.computeEigenvalues functions has been called
%           at least once.
%
%       [z,X] = solver.getEigenpairsAll()
%           Returns all computed eigenvalues z with associated eigenvectors 
%           given by the columns of X.
% 
%       [z,x,eig_diff] = solver.getEigenpairBestMatch(z_target,x_target)
%           This "selects" the eigenvalue from the computed list of
%           eigenvalues that is closest to the target value z_target and
%           then returns this eigenvalue z and its corresponding
%           eigenvector x to the caller.  When vector x_target is also
%           provided, z is chosen to be the closest vlaue to z_target such
%           that x'*x_target is not zero; if this is not possible, z will
%           simply remain the closest eigenvalue to z_target.  eig_diff is
%           the distance between z and z_target.  This function should only
%           ever be called after one of the solver.computeEigenvalues
%           functions has been called at least once.
%
%       n = solver.getDimension()
%           Returns the dimension of A.
%       
%       tf = solver.isSparseSolver()
%           Returns true if sparse eigensolves are being used.
%
%       .updateInitialVector()
%           Updates the initial vector to be used for the next iterative
%           sparse eigensolve computation, based on the strategy selected
%           by opts.v0_recycling_type.  This function should only ever be
%           called after one of the solver.computeEigenvalues functions has
%           been called at least once.
% 
%       .resetInitialVectorToRandom()
%           Sets initial vector to be used for the next iterative sparse
%           eigensolve computation to a random vector.
% 
%       .setNumberRequested(k)
%           Sets the sparse eigensolver to request k eigenvalues per solve.
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
%       tf = solver.isLeftEigenvector()
%           Returns true if eigenSolver is computing left eigenvectors.
%           Returns false if eigenSolver is computing right eigenvectors.
%
%       tf = solver.isPerturbed()
%           Returns true if eigenSolver has computed eigenvalues of A plus
%           a perturbation U*V'.  Returns false if eigenSolver has computed
%           eigenvalues of just A.
%   
%       totals = solver.getTotals()
%           Returns a struct containing the total number of eigenvalue
%           solves computed (totals.solves), and, if eigs iterations
%           counts are also available (eigsPlus is being used), the sum of
%           the all eigs iterations incurred over all of these solves
%           (totals.iters).
%        
%   See also eigenSolverOptions and eigentripleSolver.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   eigenSolver.m introduced in ROSTAPACK Version 1.0
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
        [modulus,left_evector,opts] = processOptions(varargin{:});
    catch err
        switch err.identifier
            case 'eigenSolver:invalidUserOption'
                err.throwAsCaller();
            otherwise
                error(invalidInputsMsg());
        end
    end
   
    if isnumeric(matrix) || iscell(matrix)
        try 
            matrix = matrixObject(matrix);
        catch err
            err.throwAsCaller();
        end
    else       
        assert(isfield(matrix,'getSize'),invalidMatrixMsg());
        assert(isfield(matrix,'getMatrix'),invalidMatrixMsg());
        assert(isfield(matrix,'apply'),invalidMatrixMsg());
        assert(isfield(matrix,'applyHermitian'),invalidMatrixMsg());
        assert(isfield(matrix,'isSparse'),invalidMatrixMsg());
        assert(isfield(matrix,'isReal'),invalidMatrixMsg());
    end
    
    [dim,dim2]          = matrix.getSize();
    assert(dim == dim2,'eigenSolver: matrix must be square!');
    sparse_solver       = matrix.isSparse(); 
    
    eigs_maxit          = opts.maxit;
    eigs_tol            = opts.tol;
    k_requested         = min(opts.request_k_eigenvalues,dim);
    eigs_krylov_dim     = opts.krylov_dim;
    recycling_type      = opts.v0_recycling_type;
    use_default_eigs    = opts.use_default_eigs;
                              
    z                   = [];
    X                   = [];
    sort_indices        = [];
    selected_index      = [];
    iters               = 0;
    total_iters         = 0;
    total_solves        = 0;
    is_perturbed        = false;
       
    sort_fn             = ternOp(modulus, @sortByModulus, @sortByRealPart);      
    conj_left_eigenvalues_fn    = @NOP;
    if sparse_solver   
        [applyA,applyM,set_perturbation_fn] = setupMatvecFns(matrix);
        Areal                   = matrix.isReal();
       
        setNumberRequested(k_requested);
        if isempty(eigs_krylov_dim)
            eigs_krylov_dim     = max(2*k_requested,20);
        end  
        eigs_krylov_dim         = max(1,min(eigs_krylov_dim,dim));
    
        v0          = [];
        eigs_opts   = struct(   'maxit',    eigs_maxit,         ...
                                'tol',      eigs_tol,           ...
                                'p',        eigs_krylov_dim,    ...
                                'disp',     0,                  ...
                                'isreal',   Areal               );
        sort_code   = ternOp(modulus, 'LM', 'LR');
        
        if use_default_eigs
            eigsplus_fn             = @eigs;
            compute_eigenvalues_fn  = @computeEigenvaluesEigs;
            get_totals_fn           = @getTotals;
        else
            assert( exist('eigsPlus','file') > 0,       ...
                    'Error: eigsPlus is not installed!' );
            eigsplus_fn             = @eigsPlus;
            compute_eigenvalues_fn  = @computeEigenvaluesEigsPlus;
            get_totals_fn           = @getTotalsEigsPlus;
        end
        clear_perturbation_fn       = @clearPerturbationSparse;
       
        get_initial_vector_fn       = @getRandomVector; 
        reset_initial_vector_fn     = @setInitialVectorToRandom;
        switch recycling_type
            case 0 
                recycle_initial_vector_fn   = @NOP;
                update_initial_vector_fn    = @NOP;
                reset_initial_vector_fn     = @NOP;
            case 1 
                recycle_initial_vector_fn   = @recycleSelectedEigenvector;
                update_initial_vector_fn    = @setInitialVectorToRecycled;
            case 2
                recycle_initial_vector_fn   = @recycleAllEigenvectors;
                update_initial_vector_fn    = @setInitialVectorToRecycled;
        end
  
        number_requested_fn         = @setNumberRequested;
        update_eigs_options_fn      = @updateEigsOptions;
        
    else
        [A,M,set_perturbation_fn]   = setupMatrix(matrix);
        compute_eigenvalues_fn      = @computeEigenvaluesDense;
        clear_perturbation_fn       = @clearPerturbationDense;
        get_totals_fn               = @getTotals;
       
        update_initial_vector_fn    = @NOP;  
        reset_initial_vector_fn     = @NOP;
        number_requested_fn         = @NOP;
        update_eigs_options_fn      = @NOP;
    end  
        
    solver = struct(                                                    ...
            'computeEigenvaluesOfA',        @computeEigenvaluesOfA,     ...
            'computeEigenvaluesOfAUV',      @computeEigenvaluesOfAUV,   ...
            'getEigenvalue',                @getEigenvalue,             ...
            'getEigenvector',               @getEigenvector,            ...
            'getEigenpair',                 @getEigenpair,              ...
            'getEigenpairsAll',             @getEigenpairsAll,          ...
            'getEigenpairBestMatch',        @getEigenpairBestMatch,     ...
            'getDimension',                 @() dim,                    ...
            'isSparseSolver',               @() sparse_solver,          ...
            'updateInitialVector',          update_initial_vector_fn,   ...
            'resetInitialVectorToRandom',   reset_initial_vector_fn,    ...
            'setNumberRequested',           number_requested_fn,        ...    
            'updateEigsOptions',            update_eigs_options_fn,     ...
            'isSortedByModulus',            @() modulus,                ...
            'isLeftEigenvector',            @() left_evector,           ...
            'isPerturbed',                  @isPerturbed,               ...
            'getTotals',                    get_totals_fn               );
  
    % public functions    
    function computeEigenvaluesOfA()
        is_perturbed = false;
        clear_perturbation_fn();
        computeEigenvalues();
    end

    function computeEigenvaluesOfAUV(U,V)
        is_perturbed = true;
        set_perturbation_fn(U,V);
        computeEigenvalues()
    end

    function tf = isPerturbed()
        tf = is_perturbed;
    end

    function z_out = getEigenvalue(which)
        selected_index = getActualIndex(which);
        z_out = z(selected_index);
    end

    function x_out = getEigenvector()
        x_out = X(:,selected_index);
    end

    function [z_out, x_out] = getEigenpair(which)
        selected_index = getActualIndex(which);
        z_out = z(selected_index);
        x_out = X(:,selected_index);
    end

    function [z_out, X_out] = getEigenpairsAll()
        if isempty(sort_indices) % get sort order of z
            [~,sort_indices] = sort_fn(z);
        end
        z_out = z(sort_indices);
        X_out = X(:,sort_indices);
    end

    function [z_out, x_out, eig_diff] = getEigenpairBestMatch(z_target,x)
        diffs           = abs(z_target - z);
        [~, indx]       = sort(diffs);
        selected_index  = indx(1);
        
        % z_target may be a double eigenvalue so we must check whether
        % y'*x is zero for the best matching eigenvalue in the array of 
        % candidate eigenvalues z, when x is provided.
        % If y'*x is zero, we instead will take the closest value in z such
        % that y'*x is not zero.
        if nargin > 1
            for j = 1:length(diffs)
                test_index = indx(j);
                y = X(:,test_index);
                if y'*x ~= 0
                    selected_index = test_index;
                    break
                end
            end
        end
        
        eig_diff    = diffs(selected_index);
        z_out       = z(selected_index);
        x_out       = X(:,selected_index);
        % if z is empty, so is eig_diff, so make it infinity instead
        if isempty(eig_diff)
            eig_diff = inf;
        end
    end

    function data = getTotals()
        data    = struct(   'solves',   total_solves    );
    end

    function data = getTotalsEigsPlus()
        data    = struct(   'solves',   total_solves,   ...
                            'iters',    total_iters     );
    end

    % private functions
    function computeEigenvalues()
        total_solves    = total_solves + 1;
        
        try 
            compute_eigenvalues_fn();
        catch err
            ME = MException('getStabRadBound:eigenSolver',              ...
                            'eig(s) threw an error! See below for cause');
            ME = addCause(ME,err);
            ME.throw();
        end
        if min(size(z)) == 0 
            error(  'getStabRadBound:eigenSolver',                      ...
                    'eig(s) did not return any eigenvalues!'            );
        end
        if min(size(X)) == 0
            error(  'getStabRadBound:eigenSolver',                      ...
                    'eig(s) did not return any eigenvectors!'           );
        end
        
        conj_left_eigenvalues_fn();
        sort_indices    = [];
        selected_index  = [];
    end
    
    function clearPerturbationDense()
        M = A;
    end

    function setPerturbationDense(U,V)
        M = A + U*V';
    end

    function clearPerturbationSparse()
        applyM              = applyA;
        eigs_opts.isreal    = Areal;
    end

    function setPerturbationSparse(U,V)
        applyM              = @(x) applyA(x)  + U*(V'*x);
        eigs_opts.isreal    = Areal && isreal(U) && isreal(V);
    end

    function computeEigenvaluesDense()
        [X, Z]      = eig(M);
        z           = diag(Z);
    end

    function computeEigenvaluesEigsPlus()
        eigs_opts.v0    = get_initial_vector_fn(eigs_opts.isreal);
        [X, z, iters]   = eigsplus_fn(  applyM, dim, k_requested,       ...
                                        sort_code, eigs_opts            );
        total_iters     = total_iters + iters;
    end

    function computeEigenvaluesEigs()
        eigs_opts.v0 = get_initial_vector_fn(eigs_opts.isreal);
                
        [X, Z]  = eigs(applyM,dim,k_requested,sort_code,eigs_opts);
        z       = diag(Z);
        indx    = ~isnan(z);
        z       = z(indx);
        X       = X(:,indx);
    end

    function actual_index = getActualIndex(which)
        n_converged     = length(z);
        if n_converged < 1
            actual_index = [];
            return 
        end
        if isempty(sort_indices) % get sort order of z 
            [~,sort_indices] = sort_fn(z);
        end  
        index           = max(1,min(round(which),n_converged));
        actual_index    = sort_indices(index);
    end
    
    function setNumberRequested(k)
        k_requested = max(1,min(round(k),dim-2));
    end

    function setInitialVectorToRandom()
        get_initial_vector_fn = @getRandomVector;
    end

    function setInitialVectorToRecycled()
        if ~isempty(X)
            recycle_initial_vector_fn();
            get_initial_vector_fn = @getRecycledVector;
        end
    end

    function [v] = getRandomVector(real_matrix)
        v = randn(dim,1);
        if ~real_matrix
            v = v + 1i*randn(dim,1);
        end        
    end

    function v0_out = getRecycledVector(real_matrix)
        if isempty(v0)
            v0_out = getRandomVector(real_matrix);
            return
        end
        v0_out = v0;
        if real_matrix
           v0_out = real(v0_out);
        end
    end

    function recycleSelectedEigenvector()
        v0 = X(:,selected_index);
    end

    function recycleAllEigenvectors()
        v0 = sum(X,2);
    end

    function conjugateLeftEigenvalues()
        z = conj(z);
    end

    function [A,M,set_pert_fn] = setupMatrix(matrix)
        set_pert_fn = @setPerturbationDense;
        A           = matrix.getMatrix();
        if left_evector
            A                           = A';
            set_pert_fn                 = @(U,V) set_pert_fn(V,U);
            conj_left_eigenvalues_fn    = @conjugateLeftEigenvalues;  
        end
        M = A;
    end

    function [applyA,applyM,set_pert_fn] = setupMatvecFns(matrix)
        set_pert_fn = @setPerturbationSparse;
        if left_evector
            applyA                      = matrix.applyHermitian;
            set_pert_fn                 = @(U,V) set_pert_fn(V,U);
            conj_left_eigenvalues_fn    = @conjugateLeftEigenvalues; 
        else
            applyA  = matrix.apply;
        end
        applyM      = applyA;
    end

    function updateEigsOptions(user_opts)
        allowable_fields = {    'issym',    ...
                                'tol',      ...
                                'maxit',    ...
                                'p',        ...
                                'disp'      };
                            
        for j = 1:length(allowable_fields)
            field_name = allowable_fields{j};
            if isfield(user_opts,field_name)
                eigs_opts.(field_name) = user_opts.(field_name);
            end
        end
    end
end

function [radius,left_eigenvector,opts] = processOptions(varargin)
    
    if nargin > 0 
        radius = varargin{1};
        assert(islogical(radius));
    else
        radius = false;
    end
    
    if nargin > 1 
        left_eigenvector = varargin{2};
        assert(islogical(left_eigenvector));
    else
        left_eigenvector = false;
    end
    
    if nargin > 2
        opts = varargin{3};
        if isfield(opts,'verified') && opts.verified
            return
        end
    else
        opts = [];
    end
    opts = eigenSolverOptions(opts);
end

function m = invalidInputsMsg()
m = [                                                                   ...
'eigenSolver''s format is: \n'                                          ...
'  solver = eigenSolver(A)\n'                                           ...
'  solver = eigenSolver(A,sort_by_modulus)\n'                           ...
'  solver = eigenSolver(A,sort_by_modulus,left_eigenvector)\n'          ...
'  solver = eigenSolver(A,sort_by_modulus,left_eigenvector,opts)\n'     ];
end

function m = invalidMatrixMsg()
m = [                                                                   ...
'eigenSolver requires the matrix be supplied as a matrixObject or in '  ...
'one of the formats that matrixObject accepts itself as input.'         ];
end
