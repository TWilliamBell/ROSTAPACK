function [fn0,fn1,fn2] = makeSigmaFn(A,B,C,D,E,polar,solve_type)
%   makeSigmaFn:
%       For a given spectral value set, specified by matrices A,B,C,D,E,
%       this method returns function handles to compute f(x,y), the norm of
%       the transfer function at a point given by x and y (Cartesian or
%       polar coordinates), and its gradient and Hessian with respect to x
%       and/or y.
%
%   INPUT:
%       A,B,C,D,E                   [matrices]
%           System matrices.  A must be provided while B,C,D,E can be given
%           explicitly or as [] for their shortcuts.
% 
%       polar                       [logical]
%           False:  x and y are Cartesian coordinates, z = x + 1i*y
%           True:   x and y are polar coordiantes, z = x*exp(1i*y)
%
%       solve_type                  [0,1,2]
%           Sets how inverses of zE-A should be applied:
%           0: via upper triangular Hessenberg factorization 
%           1: via LU
%           2: via LU with permutations 
% 
%           Option 0 is intended for small-scale problems, where it is
%           possible to do an upper triangular Hessenberg factorization of
%           E and A, which has O(n^3) work.  The benefit of this is that
%           all subsequent norm of transfer function evaluations can be
%           done in O(n^2) work, for any point in the complex plane.
%
%           Options 1 and 2 may be used with small or large systems.  With
%           dense matrices, each evaluation of the norm of the transfer
%           function requires O(n^3) work, since the dominant cost is doing
%           the LU decomposition.  For large sparse matrices, the work is
%           proportional to the fill-in factor in the LU matrices, which
%           hopefully remains O(n).  The LU approaches are preferrable when
%           A and E are large and sparse or if the norm of the transfer
%           function will only be evaluated very few times.
%   
%   OUTPUT:
%       fn0,fn1,fn2                 [function handles]
%           where [f,g,H,u,v] = fnX(x,y,type) 
%               f:  the norm of the transfer function 
%               g:  the first partial derivative(s) 
%               H:  the second partial derivative(s) 
%               u:  corresponding left singular vector 
%               v:  corresponding right singular vector 
%           at the point specified by x and y.  
%           The partial_type input selects what partials are returned:
%           fn0:
%               1:          g is nan, H = nan
%               2:          g is nan, H = nan
%               otherwise:  g is nan(2,1), H is nan(2)
%           fn1:
%               1:          g is f_dx, H = nan
%               2:          g is f_dy, H = nan
%               otherwise:  g is the gradient, H is nan(2)
%           fn2:
%               1:          g is f_dx, H = f_dx2
%               2:          g is f_dy, H = f_dy2
%               3:          g is the gradient, H = [f_dx2; f_dy2]
%               otherwise:  g is the gradient, H is the Hessian
%           
%           If f cannot be evaluated (e.g. at/near a pole of the transfer
%           function, f is returned as realmax while g and H are nans, 
%           consistently dimensioned. 
%
%   See also makeSigmaMinFn.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   makeSigmaFn.m introduced in ROSTAPACK Version 2.0
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
       
    fn0 = @computeFunction;
    fn1 = @computeFirstOrder;
    fn2 = @computeSecondOrder;
          
    [n,m,p]     = systemDimensions(A,B,C,D);
    if isempty(B) || (isEye(B) && ~issparse(B))
        B       = speye(n,m);
    end
    if isempty(C) || (isEye(C) && ~issparse(C))
        C       = speye(p,n);
    end
    if isempty(D) || isZero(D)
        D       = 0;    % Can be a scalar since it is only added
    end
    
    [newZ,prependZinv,appendZinv,prependE,appendE]  = setupSolver();
    [getG,getTdx,getGdxGdx2]                        = getFunctionsBC();
    
    if nargin > 5 && polar
        toComplexScalar = @partialsPolar;
        getSdy2Fn       = @getSdy2Theta;
        getSdxdyFn      = @getSdxdyTheta;
    else
        toComplexScalar = @partialsCartesian;
        getSdy2Fn       = @getSdx2;
        getSdxdyFn      = @getSdxdy;
    end
    
    % Shared variables
    u       = [];
    v       = [];
    ZinvB   = [];
    CZinv   = [];
    ZinvBv  = [];
    uhCZinv = [];
           
    % PRIVATE FUNCTIONS   
    
    % COMPUTE FUNCTION VALUE ONLY
    function [f,g,H,uo,vo] = computeFunction(x,y,partial)
                   
        z       = toComplexScalar(x,y);
        newZ(z);
        G       = getG();
        if isFiniteValued(G)
            [U,S,V] = svd(G);
            % svd can sometimes return -0 so make sure to take the absolute
            % value, to avoid 1/f being -inf instead of +inf.  This is
            % almost certain to not be an issue for the largest singular
            % value but better safe than sorry.
            f   = abs(S(1));
            u   = U(:,1);
            v   = V(:,1);
        else
            f   = realmax;
        end   
        if nargin > 2 && (partial == 1 || partial == 2)
            g   = nan;
            H   = nan;
        else
            g   = nan(2,1);
            H   = nan(2);
        end
        uo      = u;
        vo      = v;
    end
    
    % COMPUTE GRADIENT ONLY   
    function [f,g,H,uo,vo] = computeFirstOrder(x,y,partial)     
   
        [z,z_dx,z_dy] = toComplexScalar(x,y);
        newZ(z);
        G           = getG();
        if isFiniteValued(G)
            [U,S,V] = svd(G);
            % svd can sometimes return -0 so make sure to take the absolute
            % value, to avoid 1/f being -inf instead of +inf.  This is
            % almost certain to not be an issue for the largest singular
            % value but better safe than sorry.
            f       = abs(S(1));
            u       = U(:,1);
            v       = V(:,1);
            t_dx    = full(getTdx());
        else
            f       = realmax;
            t_dx    = 0;
        end
              
        % The derivative of the singular value with respect to x (the
        % rectangular coordinate) can be converted to the derivative with
        % respect to y, its radius, or its angle by multiplying it by the
        % corresponding derivatives of z_dx or z_dy.
        if nargin < 3
            partial = 0;
        end
        switch partial
            case 1
                g   = z_dx * t_dx;
                H   = nan;
            case 2
                g   = z_dy * t_dx;
                H   = nan;
            otherwise
                g   = [z_dx; z_dy] * t_dx;
                H   = nan(2);
        end
        g       = real(g);
        uo      = u;
        vo      = v;
    end

    % COMPUTE FIRST AND SECOND DERIVATIVES
    function [f,g,H,uo,vo] = computeSecondOrder(x,y,partial)
        
        [z,z_dx,z_dy,~,z_dy2,z_dxdy] = toComplexScalar(x,y);
        newZ(z);
        G           = getG();
        if isFiniteValued(G)
            [U,S,V] = svd(G);
            % svd can sometimes return -0 so make sure to take the absolute
            % value, to avoid 1/f being -inf instead of +inf.  This is
            % almost certain to not be an issue for the largest singular
            % value but better safe than sorry.  Need to do all of the
            % singular values for the second partial derivatives.
            S       = abs(S);
            f       = S(1);
            u       = U(:,1);
            v       = V(:,1);
        else
            if nargin < 3
                partial = 0;
            end
            switch partial
                case {1,2}
                    g   = 0;
                    H   = 0;
                case 3
                    g   = zeros(2,1);
                    H   = g;
                otherwise
                    g   = zeros(2,1);
                    H   = zeros(2);
            end   
            f       = realmax;
            return
        end
       
        [Gp,Gpp]    = getGdxGdx2();
        t_dx        = u'*(Gp*v);   
        [UV,dUV]    = svdToEigs(U,S,V);        
        partialFn   = secondPartialSingularValue(G,UV,dUV);
        
        % The first and second derivatives Gp and Gpp of the matrix G with
        % respect to x (the cartesian coordinate) can be converted to the
        % corresponding first and second derivatives with respect to y
        % (either the x coordinate or its radius) by multiplying Gp and Gpp
        % by the corresponding derivatives of z_dx or z_dy.  For theta, Gp
        % can be converted identically but converting Gpp also requires a
        % second term involving the converted Gp multiplied by the second
        % derivative of z with respect to theta.  The second derivatives
        % with respect to x and y can be similarly converted. Note that
        % z_dx2 is zero in either cartesian or polar form.
        if nargin < 3
            partial = 0;
        end
        switch partial
            case 1
                g       = z_dx * t_dx;              
                H       = getSdx2(partialFn,Gp,Gpp,z_dx);
            case 2
                g       = z_dy * t_dx;
                H       = getSdy2Fn(partialFn,Gp,Gpp,z_dy,z_dy2);
            case 3
                g       = [z_dx; z_dy] * t_dx;
                f_dx2   = getSdx2(partialFn,Gp,Gpp,z_dx);
                f_dy2   = getSdy2Fn(partialFn,Gp,Gpp,z_dy,z_dy2);
                H       = [ f_dx2; f_dy2 ];                
            otherwise
                g       = [z_dx; z_dy] * t_dx;
                f_dx2   = getSdx2(partialFn,Gp,Gpp,z_dx);
                f_dy2   = getSdy2Fn(partialFn,Gp,Gpp,z_dy,z_dy2);
                f_dxdy  = getSdxdyFn(partialFn,Gp,Gpp,z_dx,z_dy,z_dxdy);
                H       = [f_dx2 f_dxdy; f_dxdy f_dy2];
        end
        g       = real(g);
        H       = real(H);
        uo      = u;
        vo      = v;
    end
       
    % ABSTRACTIONS OF THE BASIC OPERATIONS 
        
    function X = prependZinvE(Y)
        X       = prependZinv(prependE(Y));
    end

    function X = appendEZinv(Y)
        X       = appendZinv(appendE(Y));
    end
    
    % COMPUTE TRANSFER MATRIX
    
    function G = getGViaB()
        ZinvB   = prependZinv(B);
        G       = full((C*ZinvB) + D);
    end

    function G = getGViaC()
        CZinv   = appendZinv(C);
        G       = full((CZinv*B) + D);
    end

    % FOR GRADIENT ONLY 

    function t_dx = getTdxViaB()
        ZinvBv  = ZinvB*v;
        uhCZinv = appendZinv(u'*C); 
        t_dx    = full(-uhCZinv * prependE(ZinvBv));
    end

    function t_dx = getTdxViaC()
        ZinvBv  = prependZinv(B*v);
        uhCZinv = u'*CZinv;
        t_dx    = full(-uhCZinv * prependE(ZinvBv));
    end
    
    % FOR GRADIENT AND HESSIAN
    
    function [Gdx,Gdx2] = getGdxGdx2ViaB()
        ZinvEZinvB  = prependZinvE(ZinvB);
        Gdx         = full(-C*ZinvEZinvB);
        Gdx2        = full(C*prependZinvE(ZinvEZinvB));
    end

    function [Gdx,Gdx2] = getGdxGdx2ViaC()
        CZinvEZinv  = appendEZinv(CZinv);
        Gdx         = full(-CZinvEZinv*B);
        Gdx2        = full(appendEZinv(CZinvEZinv)*B);
    end

    % SETUP FUNCTIONS

    function [newZ,prependZinv,appendZinv,prependE,appendE] = setupSolver()
        E_ident = isempty(E) || isEye(E);
        
        switch solve_type
            case 0
                if E_ident
                    [solv,L,R]  = hessSolver(A);
                else
                    [solv,L,R]  = hessSolver(A,E);
                    E               = R*E*L;
                end   
                C               = C*L;
                B               = R*B;
                newZ            = @solv.factor;
                prependZinv     = @solv.applyLeftUH;
                appendZinv      = @solv.applyRightUH;
            case {1,2}
                if E_ident
                    I           = speye(n);     
                    getZ        = @(z) z*I - A; 
                else
                    getZ        = @(z) z*E - A;
                end
                solv            = luSolver(solve_type == 2);
                newZ            = @(z) solv.factor(getZ(z));
                prependZinv     = @solv.applyLeft;
                appendZinv      = @solv.applyRight;
        end
       
        if E_ident
            prependE    = @(X) X;
            appendE     = @(X) X;
        else
            prependE    = @(X) E*X;
            appendE     = @(X) X*E;
        end
    end

    function [getG,getTdx,getGdxGdx2] = getFunctionsBC()
        if m <= p
            getG        = @getGViaB;
            getTdx      = @getTdxViaB;
            getGdxGdx2  = @getGdxGdx2ViaB;
        else
            getG        = @getGViaC;
            getTdx      = @getTdxViaC;
            getGdxGdx2  = @getGdxGdx2ViaC;
        end
    end
end

function s_dx2 = getSdx2(partialFn,Gp,Gpp,z_dx,~)
    % z_dx2 term is zero when x is x, y, or radius.   
    Gdx     = z_dx * Gp;
    Gdx2    = (2*z_dx^2)*Gpp;
    s_dx2   = partialFn(1,Gdx,Gdx2);
end

function s_dy2 = getSdy2Theta(partialFn,Gp,Gpp,z_dy,z_dy2)
    % z_dy2 term is not zero when y is theta.
    Gdy     = z_dy * Gp;
    Gdy2    = (2*z_dy^2)*Gpp + z_dy2*Gp;
    s_dy2   = partialFn(1,Gdy,Gdy2);
end

function s_dxdy = getSdxdy(partialFn,Gp,Gpp,z_dx,z_dy,~)
    % z_dxdy term is zero for cartesian form
    Gdx     = z_dx * Gp;
    Gdy     = z_dy * Gp;
    Gdxdy   = (2*z_dx*z_dy)*Gpp;
    s_dxdy  = partialFn(1,Gdx,Gdy,Gdxdy);
end

function s_dxdy = getSdxdyTheta(partialFn,Gp,Gpp,z_dx,z_dy,z_dxdy)
    % z_dxdy term is not zero for polar form
    Gdx     = z_dx * Gp;
    Gdy     = z_dy * Gp;
    Gdxdy   = (2*z_dx*z_dy)*Gpp + z_dxdy*Gp;
    s_dxdy  = partialFn(1,Gdx,Gdy,Gdxdy);        
end