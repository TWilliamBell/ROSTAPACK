function assertSystem(A,B,C,D,form_full)
%   assertSystem:
%       This is a utility to check whether or not the user-provided system
%       (A,B,C,D) has been correctly specified when some of the matrices
%       are given implicitly.  If all the matrices are given explicitly,
%       this routine simply verifies that the dimensions are all correct.
%
%       This method not only checks simple properties, which can (and are)
%       checked at runtime by the other methods in ROSTAPACK, but also
%       other necessary properties which cannot be checked so
%       easily/efficiently.  If some of the matrices are given implicitly,
%       either as outer products or via function handles, then the user
%       must declare the dimensions of these matrices and whether or not
%       they are real-valued.  Furthermore, in the latter case, they must
%       provide a way to apply a matrix and its conjugate transpose.  For
%       these implicitly given matrices, this utility will attempt to check
%       whether or not these user declarations are correct.
%       
%       By default, this method will attempt to perform multiplications
%       with the implicit matrices and their conjugate transposes to check:
%
%           - if there are dimensional incompatibilities
%
%           - if real-valued declared matrices actually return
%             complex-valued products when applied to real-valued inputs
%
%           - if matrices return products with infs or nans when applied to
%             finite-valued inputs
%
%           - if the function handles do indeed return the correct output
%             (requires form_full to be set to true).
%
%       These checks rely on doing matrix multiplications on random vectors
%       and matrices and, as such, it is unlikely that these checks will
%       miss errors.  Of course, this routine can be run multiple times to
%       increase the confidence of these tests.
%       
%       Finally, this method can also optionally construct full versions of
%       all the matrices, so that explicit assertions can be done without
%       any doubt.  However, this is not practical for large-scale systems,
%       due to the memory cost of storing the matrices is full/dense form.
%
%   USAGE:
%       assertSystem(A,B,C,D)
%       assertSystem(A,B,C,D,form_full)
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
%           B: [] indicates B is an appropriately-sized sparse identity
%           C: [] indicates C is an appropriately-sized sparse identity
%           D: [] indicates D is an appropriately-sized zero matrix
%           If B and D are both [], m is automatically set to n.
%           If C and D are both [], p is automatically set to n.
%
%       form_full               [optional | {false}]
%           Set this to true to form full versions of all the matrices, so
%           that all explicit checks can be done.  Note that this is likely
%           impractical for large-scale matrices.  
%
%   THROWS:
%       Throws errors if it detects that the system has been incorrectly
%       specified by the user.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   assertSystem.m introduced in ROSTAPACK Version 1.0
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

    if nargin < 5
        form_full = false;
    end
    if nargin < 4
        D = [];
    end
    if nargin < 3
        C = [];
    end
    if nargin < 2
        B = [];
    end
    
    try 
        sys = stateSpaceABCD(A,B,C,D);
        checkMatrix('A');
        if ~isempty(B)
            checkMatrix('B');
        end
        if ~isempty(C)
            checkMatrix('C');
        end
        if ~isempty(D)
            checkMatrix('D');
        end
    catch err
        err.throwAsCaller();
    end

    function checkMatrix(Mname)
        
        Mobj    = sys.getMatrixObject(Mname);
        M       = Mobj.getMatrix();
        [m,n]   = Mobj.getSize();
        
        % If the matrix was given explicitly, there is nothing to check.
        if ~isempty(M)
            return
        end
        
        % When a matrix is given as an outer product, we only need to
        % confirm the declaration that it is truly real-valued. 
        [U,~]   = Mobj.getOuterProduct();
        if ~isempty(U) && Mobj.isReal()
            assertRealValued();
            return
        end
        
        % Otherwise, if the matrix was given implicitly via function
        % handles, we need to check whether or not:
        % - the declared dimensions are correct
        % - it is truly real valued, when declared as such
        % - it contains any infs or nans.
        
        % Check whether the function handles work for vectors and matrices
        % compatible with the declared dimensions of the matrix.
        assertDimensions()
        
        % Check whether the matrix is truly real-valued if declared as such
        if Mobj.isReal()
            assertRealValued();
        end
        
        % Ensure that the implicitly supplied matrix has no infs/nans
        assertFiniteValued();

        if form_full
            [diff,diff_h] = measureFunctionHandleDiffs();
            fprintf('The %s handle has error %g.\n',Mname,diff);
            fprintf('The %s'' handle has error %g.\n',Mname,diff_h);
        end
        
        % DONE WITH CHECKS
        
        % Nested helper functions to assert properties.
        
        function assertRealValued()
            
            if form_full
                assert(isRealValued(Mobj.formFull()),notRealMsg());
                return
            end
            
            x       = randn(n,1);
            assert(isRealValued(Mobj.apply(x)),notRealMsg(false));
            y       = randn(m,1);
            assert(isRealValued(Mobj.applyHermitian(y)),notRealMsg(true));
            
            function m = notRealMsg(hermitian)
                s1 = sprintf('Matrix %s was reported as real but ',Mname);
                if nargin < 1
                    s2 = 'it actually';
                else
                    op  = ternOp(hermitian,'''*','*');
                    s2  = sprintf('%s%sx',Mname,op);
                end
                m = sprintf('%s%s contains complex-valued entries!',s1,s2);
            end
        end
        
        function assertDimensions()
            % We can't form the full matrix via applying the matrix to an
            % identity because if the declared dimensions are incorrect,
            % the multiply will fail.  So we must assert the dimensions.
            
            x = randn(n,2);
            try
                Mobj.apply(x(:,1));
            catch
                error(badFnHandleMsg(false,n));
            end
            try
                Mobj.apply(x);
            catch
                error(badFnHandleMsg(false,n,2));
            end
            
            y = randn(m,2);
            try
                Mobj.applyHermitian(y(:,1));
            catch
                error(badFnHandleMsg(true,m));
            end
            try
                Mobj.applyHermitian(y);
            catch
                error(badFnHandleMsg(true,m,2));
            end
            
            function m = badFnHandleMsg(hermitian,m,n)
                s1 = sprintf('Matrix %s error:', Mname);
                op = ternOp(hermitian,'''*','*');
                s2 = sprintf('%s%sx function handle failed for',Mname,op);
                if nargin < 3
                    s3 = sprintf('a column vector of length of %d!',m);
                else
                    s3 = sprintf('a %d by %d matrix!',m,n);
                end
                m = [s1 ' ' s2 ' ' s3];
            end
        end
        
        function assertFiniteValued()
       
            if form_full
                checkFiniteValued(Mobj.formFull());
                return
            end
           
            x       = randn(n,1);
            Mx      = Mobj.apply(x);
            checkFiniteValued(Mx,false);
            
            y       = randn(m,1);
            Mhy     = Mobj.applyHermitian(y);
            checkFiniteValued(Mhy,true);
            
            function checkFiniteValued(X,hermitian) 
                if nargin < 2
                    msg_fn = @(nans) notFiniteMsg(nans);
                else
                    msg_fn = @(nans) notFiniteMsg(nans,hermitian);
                end
                
                if isFiniteValued(X)
                    return
                end
                % Otherwise, M has infs or nans
                assert(~any(isinf(X(:))),msg_fn(false));
                assert(~any(isnan(X(:))),msg_fn(true));
            end
            
            function m = notFiniteMsg(nans,hermitian)
                type = ['contains ' ternOp(nans,'nans','infs')];
                
                if nargin < 2
                    m = sprintf('Matrix %s %s!',Mname,type);
                    return
                end
                
                s1 = sprintf('Matrix %s error:', Mname);
                op = ternOp(hermitian,'''*','*');
                s2 = sprintf('the result of %s%sx %s!',Mname,op,type);
                m = [s1 ' ' s2];
            end
        end
        
        function [diff,diff_h] = measureFunctionHandleDiffs()
            M       = Mobj.formFull();
            diff    = norm(M-Mobj.apply(speye(n)),'fro');          
            diff_h  = norm(M'-Mobj.applyHermitian(speye(m)),'fro'); 
        end
    end
end