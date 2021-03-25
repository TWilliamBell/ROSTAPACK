function [fn0,fn1,fn2] = makeSigmaMinFn(A,E,polar,varargin)
%   makeSigmaMinFn:
%       For a given spectral value set with only A and E, i.e. B=C=I, D=0,
%       and n=m=p, this method returns function handles to compute f(x,y),
%       the norm of the transfer function at a point given by x and y
%       (Cartesian or polar coordinates), and its gradient and Hessian with
%       respect to x and/or y.

%       However, since B=C=I and D=0, this routine does this by instead
%       computing the reciprocal of the smallest singular value of zE - A,
%       which is equivalent to computing the norm of the transfer function.
%       As this alternative formula does not involve an inverse, it should
%       be more efficient.  The returned function also optionally computes
%       the first and second partial derivatives of f(x,y) with respect to
%       x and/or y.
%
%       NOTE: When derivatives are requested and n > 25, svd will call
%       LAPACK's GESDD (divide and conquer) routine.  Unfortunately, for
%       poorly-scaled matrices, GESDD can sometimes be much less accurate
%       than GESVD when the smallest singular values are desired.  For such
%       problems, it may be necessary to instead use makeSigmaFn.
%
%   INPUT:
%       A,E                         [matrices]
%           System matrices.  A must be provided while E can be given
%           explicitly or as [] to denote the identity matrix
% 
%       polar                       [logical]
%           False:  x and y are Cartesian coordinates, z = x + 1i*y
%           True:   x and y are polar coordiantes, z = x*exp(1i*y)
%   
%   OUTPUT:
%       fn0,fn1,fn2                 [function handles]
%           where [f,g,H,u,v] = fnX(x,y,partial_type) 
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
%   See also makeSigmaFn.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   makeSigmaMinFn.m introduced in ROSTAPACK Version 2.0
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
       
    n           = length(A);
    if isempty(E) || isEye(E)
        I       = speye(n);
        getZ    = @(z) z*I - A;
        E       = I;
    else
        getZ    = @(z) z*E - A;
    end
    
    if nargin > 2 && polar
        toComplexScalar = @partialsPolar;
    else
        toComplexScalar = @partialsCartesian;
    end
 
    % PRIVATE FUNCTIONS   
    
    % COMPUTE FUNCTION VALUE ONLY
    function [f,g,H,uo,vo] = computeFunction(x,y,partial)
                   
        z       = toComplexScalar(x,y);
        [U,S,V] = svd(getZ(z));
        % svd can sometimes return -0 for exactly zero singular values so
        % make sure to take the absolute value, to avoid 1/S(end)
        % inadvertently being -inf instead of +inf.
        f       = 1/abs(S(end));
        if nargin > 2 && (partial == 1 || partial == 2)
            g   = nan;
            H   = nan;
        else
            g   = nan(2,1);
            H   = nan(2);
        end
        % We want to return the left and right singular vectors for the
        % largest singular value of (zE - A)^{-1} so U and V are swapped!
        uo      = V(:,end);
        vo      = U(:,end);
    end

    % COMPUTE GRADIENT ONLY   
    function [f,g,H,uo,vo] = computeFirstOrder(x,y,partial)     
   
        % For consistency with makeSigmaFn, we will first just compute the
        % minimum singular value and the first partial derivative(s).  Then
        % we will convert these all to be for the reciprocal of the minimum
        % singular value.
        [z,z_dx,z_dy]   = toComplexScalar(x,y);
        [U,S,V]         = svd(getZ(z));
        % svd can sometimes return -0 so make sure to take the absolute
        % value, to avoid 1/f being -inf instead of +inf.
        f               = abs(S(end));
        u               = U(:,end);
        v               = V(:,end);
        t_dx            = u'*E*v;
           
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
                H   = nan(2,1);
        end
        g       = real(g);
        % Do the conversion from smallest singular value to its reciprocal.
        [f,g]   = toReciprocal1st(f,g);
        % Remember U and V are swapped for (zE - A)^{-1}!
        uo      = V(:,end);
        vo      = U(:,end);
    end

    % COMPUTE FIRST AND SECOND DERIVATIVES
    function [f,g,H,uo,vo] = computeSecondOrder(x,y,partial)
        
        % For consistency with makeSigmaFn, we will first just compute the
        % minimum singular value and the first and second partial
        % derivative(s).  Then we will convert these all to be for the
        % reciprocal of the minimum singular value.
        [z,z_dx,z_dy,z_dx2,z_dy2,z_dxdy] = toComplexScalar(x,y);

        Z           = getZ(z);
        [U,S,V]     = svd(Z);
        % svd can sometimes return -0 so make sure to take the absolute
        % value, to avoid 1/f being -inf instead of +inf.  Need to do all
        % of the singular values for the second partial derivatives.
        S           = abs(S);
        f           = S(end);
        u           = U(:,end);
        v           = V(:,end);
        t_dx        = u'*E*v;
       
        [UV,dUV]    = svdToEigs(U,S,V);
        partialFn   = secondPartialSingularValue(Z,UV,dUV);
        
        % The first and second derivatives Gp and Gpp of the matrix G with
        % respect to x (the rectangular coordinate) can be converted to the
        % corresponding first and second derivatives with respect to y
        % (either the rectangular coordinate or its radius) by multiplying
        % Gp and Gpp by the corresponding derivatives of z_dx or z_dy.  For
        % theta, Gp can be converted identically but converting Gpp also
        % requires a second term involving the converted Gp multiplied by
        % the second derivative of z with respect to theta.  The second
        % derivatives with respect to x and y can be similarly converted.  
        % Note that z_dx2 is zero in either rectangular or polar form.
        if nargin < 3
            partial = 0;
        end
        switch partial
            case 1
                g       = z_dx * t_dx;  
                H       = partialFn(n,z_dx*E,z_dx2*E);
            case 2
                g       = z_dy * t_dx;
                H       = partialFn(n,z_dy*E,z_dy2*E);
            case 3
                g       = [z_dx; z_dy] * t_dx;
                f_dx2   = partialFn(n,z_dx*E,z_dx2*E);
                f_dy2   = partialFn(n,z_dy*E,z_dy2*E);
                H       = [ f_dx2; f_dy2 ];
            otherwise
                zE_dx   = z_dx*E;
                zE_dy   = z_dy*E;
                g       = [z_dx; z_dy] * t_dx;
                f_dx2   = partialFn(n,zE_dx,z_dx2*E);
                f_dy2   = partialFn(n,zE_dy,z_dy2*E);
                f_dxdy  = partialFn(n,zE_dx,zE_dy,z_dxdy*E);
                H       = [f_dx2 f_dxdy; f_dxdy f_dy2];
        end
        g       = real(g);
        H       = real(H);
        % Do the conversion from smallest singular value to its reciprocal.
        [f,g,H] = toReciprocal2nd(f,g,H);
        % Remember U and V are swapped for (zE - A)^{-1}!
        uo      = V(:,end);
        vo      = U(:,end);
    end

end

function [f,g] = toReciprocal1st(s,grad_s)
    % Given function value s and its gradient, this returns 1/s and the
    % gradient of 1/s.
    f   = 1/s;
    g   = -grad_s/(s^2);
end

function [f,g,H] = toReciprocal2nd(s,grad_s,hess_s)
    % Given function value s, its first partial derivatives, and second
    % partial derivatives, this returns 1/s and corresponding first and
    % second partial derivatives for 1/s.
    s2  = s^2;
    s3  = s2*s;
    f   = 1/s;
    g   = -grad_s/s2;
    if length(grad_s) == 1
        % Single variable problem
        sx  = grad_s;
        sxx = hess_s;

        H   = (2*sx^2 - sxx*s)/s3;
    elseif numel(hess_s) == 2
        % Two variable problem but only given dx2 and dy2 so H is a vector
        % and does not include the dxdy term      
        sx  = grad_s(1);
        sy  = grad_s(2);
        sxx = hess_s(1);
        syy = hess_s(2);

        fxx = (2*sx^2 - sxx*s)/s3;
        fyy = (2*sy^2 - syy*s)/s3;
        H   = [fxx; fyy];
    else
        % Two variable problem with complete Hessian so return gradient
        % and Hessian of 1/s   
        sx  = grad_s(1);
        sy  = grad_s(2);
        sxx = hess_s(1,1);
        syy = hess_s(2,2);
        sxy = hess_s(1,2);

        fxx = (2*sx^2 - sxx*s)/s3;
        fyy = (2*sy^2 - syy*s)/s3;
        fxy = (2*sx*sy - sxy*s)/s3;
        H   = [fxx fxy; fxy fyy];
    end
end