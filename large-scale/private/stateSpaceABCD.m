function sys = stateSpaceABCD(A,B,C,D,varargin)
%   stateSpaceABCD:
%       An "object" for representing an LTI dynamical system, in
%       state-space form, along with an appropriate eigenvalue solver
%       (dense or sparse) and an appropriate root function for determining
%       the stability of the perturbed system matrix M(Delta) for
%       perturbations Delta = epsilon*U*V', where U*V' is unit norm and
%       
%           M(Delta) = A + B Delta (I - D Delta)^{-1} C.
%
%   INPUT:
%       System matrix A         [required]
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
%       System matrices B,C,D   [required]
%           Matrices B,C,D can also be provided in any of above formats,
%           as well as via the following shortcuts:
%           B: [] indicates B is an appropriately-sized sparse identity.
%           C: [] indicates C is an appropriately-sized sparse identity.
%           D: [] indicates D is an appropriately-sized zero matrix.
%           If B and D are both [], m is automatically set to n.
%           If C and D are both [], p is automatically set to n.
% 
%       Note: when explicitly providing B=C=I and D=0, B and C should be
%       given as sparse identities (using speye or via function handles)
%       and D should be given as a sparse zero matrix.  Otherwise,
%       computation and memory usage may be very inefficient.
%
%       opts                    [optional: struct of parameters]
%           An optional struct of settable parameters or [].
%           The two main parameters control:
%           - whether the system is considered continuous or discrete time 
%           - whether dense or sparse eigenvalue solves will be used.
%   
%   OUTPUT:
%       system
%           An "object", a struct containing the following functions
%           relevant to the continuous-time or discrete-time, dense or
%           sparse, LTI system:
%   
%       [n,m,p] = system.getDimensions()
%           Returns the dimensions of the LTI, where:
%               - A is n by n 
%               - B is n by m
%               - C is p by n
%               - D is p by m
%       
%       M_obj = system.getMatrixObject(mat_name)
%           Returns a matrixObject representation of the desired matrix.
%           mat_name must be one of 'A', 'B', 'C', or 'D'.
%
%       tf = system.isReal()
%           Returns true if all the matrices only have real-valued entries.
%           If some/all of the matrices were given implicitly, this returns
%           true based on the user's declarations of whether or not these
%           implicitly-provided matrices only contain real-valued entries.
%
%       tf = system.isDiscrete()
%           Returns true if the system is discrete time.  Otherwise returns
%           false when the system is continuous time.
%
%       fn = system.getRootFunction()
%           Returns a function handle for computing the root function
%           determining stability.  This function handle takes a single
%           real-valued or complex-valued scalar z as its input.  The root
%           function is:
%           - continuous-time systems:  real(z)
%           - discrete-time systems:    abs(z) - 1.
%
%       solver = system.getEigSolver()
%           Returns an eigentripleSolver instance, preconfigured for this
%           system (dense or sparse eigenvalue solves and continuous-time
%           or discrete-time eigenvalue ordering).
% 
%       counts = system.getCounts()
%           This returns a struct with the following fields: A, Ah, B, Bh,
%           C, Ch, D, and Dh, where each contains the number of times each
%           matrix operator has been applied.  Ah means A'.  If
%           opts.count_multiplies was set to false (or not provided), all
%           the field values will be nan.
%           
% 
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   stateSpaceABCD.m introduced in ROSTAPACK Version 1.0
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
        opts            = processOptions(varargin{:});
    catch err
        err.throwAsCaller();
    end
    discrete_time       = opts.discrete_time;
    sparse_mode         = opts.sparse_mode;
    eig_solver_opts     = opts.eig_solver_opts;
    count_multiplies    = opts.count_multiplies;
    makeMatrix          = @(M) matrixObject(M,count_multiplies);
   
    try 
        current = 'A';
        A           = makeMatrix(A);
        [n,Ac]      = A.getSize();
        assert(n == Ac,'Matrix A must be square!');
        [p,m]       = deal(n);
        
        current = 'B';
        if ~isempty(B)
            B       = makeMatrix(B);
            [~,m]   = B.getSize();
        end
        
        current = 'C';
        if ~isempty(C)
            C       = makeMatrix(C);
            [p,~]   = C.getSize();
        end
        
        current = 'D';
        if ~isempty(D)
            D       = makeMatrix(D);
            [p,m]   = D.getSize();
        end
    catch err
        processError(current);
        err.throw();
    end
    
    if isempty(B)
        B = makeMatrix(getNonSquareIdentity(n,m));
    end
    if isempty(C)
        C = makeMatrix(getNonSquareIdentity(p,n));
    end
    if isempty(D)
        D = makeMatrix(sparse(p,m));
    end
          
    [Br,Bc] = B.getSize();
    [Cr,Cc] = C.getSize();
    [Dr,Dc] = D.getSize();
    
    assert(n == Br,'A and B must have the same number of rows!');
    assert(Cr == Dr,'C and D must have the same number of rows!');
    assert(Ac == Cc,'A and C must have the same number of columns!');
    assert(Bc == Dc,'B and D must have the same number of columns!');
    
    ABCD_real   = A.isReal() && B.isReal() && C.isReal() && D.isReal();  
    root_fn     = ternOp(   discrete_time,              ...
                            @rootFunctionModulus,       ...
                            @rootFunctionRealPart       );
    
    if sparse_mode == 0 && A.isSparse()     
        % the user wants to force that dense eig computations are used
        A.formFull();
    elseif sparse_mode == 1 && ~A.isSparse()
        % the user wants to force that sparse eigs computations are used
        A.removeFull();
    end
    
    % Set up the eigentriple solver      
    try 
        eig_solver  = eigentripleSolver(A,discrete_time,eig_solver_opts);
    catch err
        % invalid eig_solver_opts parameter set by user is most likely
        err.throwAsCaller();
    end
    
    sys = struct(   'getDimensions',    @getDimensions,         ...
                    'getMatrixObject',  @getMatrixObject,       ...
                    'isReal',           @() ABCD_real,          ...
                    'isDiscreteTime',   @() discrete_time,      ...
                    'getRootFunction',  @() root_fn,            ...
                    'getEigSolver',     @() eig_solver,         ...
                    'getCounts',        @getCounts              );
    
    function [n,m,p] = getDimensions()
        [n,~]   = A.getSize();
        [p,m]   = D.getSize();
    end
    
    function M_data = getMatrixObject(mat_name)
        mat_name = upper(mat_name);
        switch mat_name
            case {'A','B','C','D'}
                M_data = eval(mat_name);
            otherwise
                error('System matrix must be one of A,B,C,D.');
        end
    end
    
    function c = getCounts()
        if count_multiplies
            c = makeCountStruct(    A.getCounts(),  B.getCounts(),      ...
                                    C.getCounts(),  D.getCounts()       );
        else
            c = makeCountStruct([nan nan],[nan nan],[nan nan],[nan nan]);
        end
    end

    function processError(name)
        id          = 'stateSpaceABCD:invalidInput';
        msg_lines   = strsplit(err.message,'\n');
        msg         = [ msg_lines{1} sprintf(': Matrix %s\n',name)      ...
                        msg_lines{2} '\n'                               ...
                        msg_lines{3}                                    ];
        err         = MException(id,msg);
    end
end

function s = makeCountStruct(A,B,C,D)
    s = struct( 'A',A(1),   'Ah',A(2),  'B',B(1),   'Bh',B(2),          ...
                'C',C(1),   'Ch',C(2),  'D',D(1),   'Dh',D(2)           ); 
end

function opts = processOptions(opts)
    if nargin > 0
        if isfield(opts,'verified') && opts.verified
            return
        end
    else
        opts = [];
    end
    opts = stateSpaceABCDOptions(opts);
end

function mat = getNonSquareIdentity(m,n)
    k       = min(m,n);
    l       = max(m,n);
    k_fn    = @(x) x(1:k,:);
    l_fn    = @(x) [x(1:k,:); zeros(l-k,size(x,2))];
    if m < n
        [I_fn,It_fn] = deal(k_fn,l_fn);
    else
        [I_fn,It_fn] = deal(l_fn,k_fn);
    end
    mat     = {I_fn,It_fn,m,n,true};
end