function obj = matrixObject(M,count)
%   matrixObject:
%       This object encapsulates and abstracts matrices so that operations
%       involving a given matrix M can be done via matrixObject's
%       interface, regardless of how the matrix was originally supplied.
%       Thus, code that uses a matrixObject does not need specific code
%       paths to handle different matrix representations.
%
%       Note that matrixObject requires at least MATLAB R2016a when
%       computing the norm of sparse or implicitly-given matrices, since it
%       uses the recently added ability to provide svds a matrix via a
%       function handle representation.
%
%   INPUT:
%       M           [required matrix]
%           The following formats are supported: 
%           - as an explicit matrix, dense or sparse: M
%           - as an explicit matrix, dense or sparse, enclosed in single
%             element cell array: {M}
%           - as an outer product M = U*V': {U,V}
%           - as a single function handle for multiplying M and M':
%                   {applyM,rows,cols,is_real}
%             where:
%                   applyM(x,false) returns M*x
%                   applyM(x,true) returns M'*x 
%                   rows is the number of rows of matrix M
%                   cols is the number of columns of matrix M
%                   is_real:logical, true if M only contains real entries
%           - as two separate function handles for multiplying M and M':
%                   {applyM,applyMh,rows,cols,is_real}
%             where:
%                   applyM(x,false) returns M*x
%                   applyMh(x,true) returns M'*x 
%                   rows is the number of rows of matrix M
%                   cols is the number of columns of matrix M
%                   is_real: logical, true if M only contains real entries
%
%           Note that the function handles should be able to multiply by M
%           and M' when x is either a vector or matrix.
%
%       count       [optional logical | {false}]
%           When this is provided and true, the matrixObject will also keep
%           a count of how many times multiplications have each been done
%           with M and M'.  Note that enabling this adds some overhead.  
% 
%   OUTPUT:
%       obj         [matrixObject representation of matrix M]
%           An "object", a struct containing the following functions:
%       
%           dim = obj.getSize()
%           [m,n] = obj.getSize()
%           [m,n,s] = obj.getSize()
%               Returns the dimensions of the m by n matrix:
%               - dim = [m n]
%               - [m,n,s], where s is the common dimension if the matrix 
%                 was given as an outer product.  If not, s is set to [].
%
%           M = obj.getMatrix()
%               Returns an explicit representation of the matrix, either
%               full or sparse.  If an explicit representation is not 
%               currently available, this will return [].  See
%               obj.formFull() and obj.removeFull().
%
%           [U,V] = obj.getOuterProduct()
%               If the matrix was given as an outer product U*V', this will
%               returns those two factors.  If not, U and V will each be
%               returned as []. 
%
%           value = obj.getNorm()                
%               This returns the two norm of the matrix.  If the matrix is
%               full, norm() is used to perform the computation.
%               Otherwise, svds is used and R2016a or later is then
%               required.  Note that the norm computation is only ever
%               performed so subsequent calls to this function will just
%               return the saved value.
%
%           b = obj.apply(x)      
%               Returns b = M*x, where x may be a vector or matrix
%
%           b = obj.applyHermitian(x)
%               Returns b = M'*x, where x may be a vector or matrix
%
%           b = obj.innerProduct(v,u)
%               Return b = v'*M*u, where v and u may be vectors or
%               matrices, provided that they are dimensionally compatible.
%               The order of evaluation is chosen to minimize the number of
%               flops.
%
%           t = obj.isSparse()
%               Returns zero if a full matrix representation is currently
%               available or if not, a positive number representing the
%               type of sparse representation that is available:
%                   0 - full matrix available
%                   1 - sparse matrix is available
%                   2 - outer product representation is available
%                   3 - only function handles
%               See obj.formFull() and obj.removeFull().
%               
%           tf = obj.isReal()
%               Returns true if the matrix only contains real-valued
%               entries.  If the matrix was supplied as an outer product
%               U*V', this will return true if both U and V are real, even
%               though it is possible that U*V' may be real even if U
%               and/or V aren't.  Furthermore, if the matrix was supplied
%               implicitly via a function handle, this will simply return
%               the user's declaration regarding whether or not the matrix
%               is real-valued.
%
%           tf = obj.isZero()
%               Returns true if the matrix only has zero entries.  However,
%               if the matrix was supplied implicitly via a function
%               handle, it is assumed the matrix is not zero adn this will
%               always return false.
%
%           M = obj.formFull()
%               Creates a full (dense) representation of the matrix, if it
%               does already exist and returns it to the caller.  Note that
%               this will often be a very expensive operation.  After this
%               routine is called, the full representation will be
%               obtainable via obj.getMatrix() and obj.isSparse() will
%               return 0.
%
%           .removeFull()
%               If a full (dense) representation has been created of a
%               sparse or implicit form, by calling obj.formFull(), this
%               will routine will delete that full representation and
%               revert the matrixObject's type back to its original type
%               value.  However, if the original matrix was already a full
%               matrix, then nothing is actually deleted internally.
%               Instead, the matrixObject will simply act as if only a
%               function handle representation now is available and
%               obj.isSparse() will return 3 and calling obj.formForm() 
%               will return the matrixObject back to its original full
%               matrix type state.
%               
%          counts = obj.getCounts()
%               Returns a 1 by 2 array containing the number of
%               multiplications applied for M and M', respectively stored
%               in counts(1) and counts(2).  If input arg count was either
%               not provided or was false, this will return [0 0].
%            
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   matrixObject.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  matrixObject.m                                                       |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is originally from URTM.                                   |
% |                                                                       |
% |  URTM is free software: you can redistribute it and/or modify         |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  URTM is distributed in the hope that it will be useful,              |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
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
    
    if nargin < 2
        count = false;
    end
    
    try 
        assert(nargin > 0);
        if ~iscell(M)
            M = {M};
        end
        assert(numel(M) == length(M));
        obj = buildMatrixObject(M,count);
    catch err
        msg = err.message;
        switch err.identifier
            case 'matrixObject:explicitMatrixError'
                type    = 'explicit matrix';
            case 'matrixObject:outerProductMatrixError'
                type    = outerProductStr();
            case 'matrixObject:functionMatrixError'
                type    = functionHandleStr();
            otherwise
                type    = 'unrecognized format';
                msg     = sprintf(invalidMatrixFormatMsg());
        end                                            
        error('matrixObject:invalidInput',                              ...
              'Invalid input\nMatrix format: %s\nIssue: %s',type,msg   );
    end   
end

function obj = buildMatrixObject(Mc,count)
    
    M           = [];
    M_left      = [];
    M_right     = [];
    common      = [];
    
    switch numel(Mc)
        case 1
            [   M,M_fn,Mh_fn,                   ...
                rows,cols,is_real,is_zero       ] = matrixExplicit(Mc);
            type = double(issparse(M));
        case 2
            [   M_left,M_right,M_fn,Mh_fn,      ...
                rows,cols,common,is_real,is_zero] = matrixOuterProduct(Mc);
            type = 2;
        case {4,5}
            [   M_fn,Mh_fn,                     ...
                rows,cols,is_real,is_zero       ] = matrixImplicit(Mc);
            type = 3;
        otherwise
            error('Invalid number of inputs!');
    end
    
    type_original               = type;                        
    cached_norm                 = [];
    [inner_fn,get_counts_fn]    = setupApplyFns();
            
    obj = struct(   'getSize',          @getSize,           ...
                    'getMatrix',        @getMatrix,         ...
                    'getOuterProduct',  @getOuterProduct,   ...
                    'getNorm',          @computeNorm,       ...          
                    'apply',            M_fn,               ...
                    'applyHermitian',   Mh_fn,              ...
                    'innerProduct',     inner_fn,           ...
                    'isSparse',         @getSparsityType,   ...
                    'isReal',           @() is_real,        ...
                    'isZero',           @() is_zero,        ...
                    'formFull',         @formFull,          ...
                    'removeFull',       @removeFull,        ...
                    'getCounts',        get_counts_fn       );        
    
    % Functions that the user may call
    
    function t = getSparsityType()
        t = type;
    end

    function A = getMatrix()
        A = M;
    end

    function [U,V] = getOuterProduct()
        U = M_left;
        V = M_right;
    end
            
    function [m,n,s] = getSize()
        if nargout < 2
            m = [rows cols];
        else
            m = rows;
            n = cols;
            s = common;
        end
    end

    function value = computeNorm()
        if ~isempty(cached_norm)
            value = cached_norm;
            return
        end
        
        if is_zero
            value = 0;
        else
            if type == 0 
                value = norm(M);
            else
                % I believe svds support for function handles was only
                % added in R2016a so this won't work on earlier releases
                value = max(svds(@applyForSvds,[rows,cols],6));
            end
        end
        cached_norm = value;
    end

    function M_out = formFull()
        if type == 0
            M_out = M;
            return
        end
        
        % Allow the norm to be recomputed by the dense method 
        cached_norm = [];
       
        switch type
            case 1
                form_fn = @() full(M);
                msg     = 'from a sparse matrix';
            case 2
                form_fn = @() M_left*M_right';
                msg     = 'from outer product UV''';
            case 3
                form_fn = @() M_fn(speye(cols));
                msg     = 'from applyM function handle';
        end
         
        type = 0;
        if cols > 100 && rows > 100
            warning(                                                    ...
                'matrixObject:formFull',                                ...
                'Forming a full matrix %s may be very expensive!',msg   );
        end
        M = form_fn();
        M_out = M;
    end

    function removeFull()
        if type == type_original
            % If the original type is a full matrix, simulate a sparse
            % matrix by setting M to [] so that only the function handles
            % for multiplying by M and M' remain (thus it becomes type 3).
            % This object will then act like a sparse matrix.
            if type == 0
                M       = [];
                type    = 3;
            end
            return
        end
        
        type = type_original;
        if type == 1
            M = sparse(M);
        else
            M = [];
        end
    end

    % Setup helper function (to be only called once internally)
    function [inner_fn,get_counts_fn] = setupApplyFns()
        if count
            [M_fn, M_count_fn]      = functionCounter(M_fn);
            [Mh_fn, Mh_count_fn]    = functionCounter(Mh_fn);
            get_counts_fn           = @getCounts;
        else
            get_counts_fn           = @() [0 0];
        end
        
        function counts = getCounts()
            counts  = [ M_count_fn() Mh_count_fn() ];
        end
        
        % Set the inner product function; do the variant which is cheaper 
        if cols < rows
            inner_fn = @(v,u) Mh_fn(v)'*u;
        else
            inner_fn = @(v,u) v'*M_fn(u);
        end
    end
    
    % We need a specific multiply function for svds
    function b = applyForSvds(x,tflag)
        if strcmp(tflag,'notransp')
            b = M_fn(x);
        else
            b = Mh_fn(x);
        end
    end

end

function [M,M_fn,Mh_fn,rows,cols,is_real,is_zero] = matrixExplicit(M)

    M           = M{1};
    [rows,cols] = size(M);
    
    assertMatrix(M,[],'explicitMatrixError');
    
    M_fn        = @(x) M*x;
    Mh_fn       = @(x) M'*x;
    is_real     = isreal(M);
    is_zero     = nnz(M) == 0;
end

function [U,V,M_fn,Mh_fn,r,c,sU,is_real,is_zero] = matrixOuterProduct(M)
    
    U       = M{1};
    V       = M{2};           
    % r = rows, c = cols, s = shared dimension
    [r,sU]  = size(U);
    [c,sV]  = size(V);
    
    assertMatrix(U,'U','outerProductMatrixError');
    assertMatrix(V,'V','outerProductMatrixError');
  
    assert( sU == sV && sU > 0,                                         ...
            'matrixObject:outerProductMatrixError',                     ...
            'Matrices U and V are not dimensionally compatible!'        );

    M_fn    = @(x) U*(V'*x);
    Mh_fn   = @(x) V*(U'*x);
    is_real = isreal(U) && isreal(V);
    is_zero = min(nnz(U),nnz(V)) == 0;
end

function [M_fn,Mh_fn,rows,cols,is_real,is_zero] = matrixImplicit(M)

    apply_fn    = M{1};
    rows        = M{end-2};
    cols        = M{end-1};
    is_real     = M{end};
    is_zero     = false;
    
    assert( isa(apply_fn,'function_handle'),                            ...
            'matrixObject:functionMatrixError',                         ...
            'applyM must be a function handle!'                         );
    
    if numel(M) < 5
        M_fn    = @(x) apply_fn(x,true);
        Mh_fn   = @(x) apply_fn(x,false);
    else
        apply_h_fn = M{2};
        assert( isa(apply_h_fn,'function_handle'),                      ...
                'matrixObject:functionMatrixError',                     ...
                'applyM must be a function handle!'                     );   
        M_fn    = apply_fn;
        Mh_fn   = apply_h_fn;
    end
        
    assert( isAnInteger(rows) && rows > 0,                              ...
            'matrixObject:functionMatrixError',                         ...
            '# of rows must be a positive integer!'                     );
    assert( isAnInteger(cols) && cols > 0,                              ...
            'matrixObject:functionMatrixError',                         ...
            '# of cols must be a positive integer!'                     );
    assert( islogical(is_real) && numel(is_real) == 1,                  ...
            'matrixObject:functionMatrixError',                         ...
            'is_real must be a logical!'                                );
end

function assertMatrix(M,matrix_name,err_name)

    if nargin < 2 || isempty(matrix_name)
        name        = 'Matrix';
    else
        name        = sprintf('Matrix %s',matrix_name);
    end
    
    if nargin < 3
        err_name    = 'invalidInput';
    end
    
    assertFn = @(tf,msg) assert(tf,['matrixObject:' err_name],msg,name);
  
    assertFn(ismatrix(M),'%s must be a matrix!');
    assertFn(~isempty(M),'%s must not be empty!');
    assertFn(isFiniteValued(M),'%s must not contain infs or nans!');
  
end

function s = outerProductStr()
    s = 'outer product UV'' -> {U,V}';
end

function s = singleFunctionHandleStr()
    s = 'single function handle -> {applyM,rows,cols,is_real}';
end

function s = twoFunctionHandleStr()
    s = 'two function handles -> {applyM,applyMh,rows,cols,is_real}';
end

function s = invalidMatrixFormatMsg()
s = [                                                                   ...
'A matrixObject can be specified in the following ways:\n'              ...
'- as an explicit matrix M or a single element cell array {M}\n'        ...
'- as an ' outerProductStr() '\n'                                       ...
'- as a ' singleFunctionHandleStr() '\n'                                ...
'- as ' twoFunctionHandleStr() '.'                                      ];
end
